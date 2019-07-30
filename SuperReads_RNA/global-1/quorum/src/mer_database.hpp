/* Quorum
 * Copyright (C) 2012  Genome group at University of Maryland.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __QUORUM_MER_DATABASE_HPP__
#define __QUORUM_MER_DATABASE_HPP__

#include <fstream>

#include <jellyfish/file_header.hpp>
#include <jellyfish/large_hash_array.hpp>
#include <jellyfish/atomic_bits_array.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/rectangular_binary_matrix.hpp>
#include <jellyfish/err.hpp>
#include <jellyfish/locks_pthread.hpp>
#include <jellyfish/atomic_field.hpp>

#include <src/verbose_log.hpp>

namespace err = jellyfish::err;

using jellyfish::mer_dna;
typedef jellyfish::large_hash::array<mer_dna> mer_array;
typedef jellyfish::large_hash::array_raw<mer_dna> mer_array_raw;
typedef jellyfish::atomic_bits_array<uint64_t> val_array;
typedef jellyfish::atomic_bits_array_raw<uint64_t> val_array_raw;

class database_header : public jellyfish::file_header {
public:
  database_header() : jellyfish::file_header() { }
  database_header(std::istream& is) : jellyfish::file_header(is) { }

  void bits(uint32_t b) { root_["bits"] = (Json::UInt)b; }
  uint32_t bits() const { return root_["bits"].asUInt(); }

  size_t value_bytes() const { return root_["value_bytes"].asLargestUInt(); }
  void value_bytes(size_t bytes) { root_["value_bytes"] = (Json::UInt64)bytes; }

  size_t key_bytes() const { return root_["key_bytes"].asLargestUInt(); }
  void key_bytes(size_t bytes) { root_["key_bytes"] = (Json::UInt64)bytes; }

  void set_format() {
    this->format("binary/quorum_db");
  }
  bool check_format() const {
    return "binary/quorum_db" == this->format();
  }
};

class hash_with_quality {
  mer_array*                         keys_;
  mer_array*                         new_keys_;
  val_array*                         vals_;
  val_array*                         new_vals_;
  const uint64_t                     max_val_;
  jellyfish::locks::pthread::barrier size_barrier_;
  jflib::atomic_field<uint16_t>       done_threads_;
  jflib::atomic_field<uint16_t>        size_thid_;
  const uint16_t                     nb_threads_;
  enum status { OK, DONE, FULL };

public:
  hash_with_quality(size_t size, uint16_t key_len, int bits, uint16_t nb_threads, uint16_t reprobe_limit = 126) :
    keys_(new mer_array(size, key_len, 0, reprobe_limit)),
    new_keys_(0),
    vals_(new val_array(bits + 1, keys_->size())),
    new_vals_(0),
    max_val_((uint64_t)-1 >> (sizeof(uint64_t) * 8 - bits)),
    size_barrier_(nb_threads),
    done_threads_(0), size_thid_(0),
    nb_threads_(nb_threads)
  { }

  ~hash_with_quality() {
    delete keys_;
    delete vals_;
  }

  bool add(const mer_dna& key, unsigned int quality) {
    bool is_new;
    size_t id;
    while(__builtin_expect(!keys_->set(key, &is_new, &id), 0)) {
      if(handle_full_ary() == FULL)
        return false;
    }

    auto     entry = (*vals_)[id];
    uint64_t nval = entry.get();;
    do {
      if((nval & 1) < quality)
        nval = 3;
      else if((nval >> 1) == max_val_ || (nval & 1) > quality)
        return true;
      else
        nval += 2;
    } while(__builtin_expect(!entry.set(nval), 0));
    return true;
  }

  void write(std::ostream& os, database_header* header = 0) const {
    if(header) {
      header->set_format();
      header->update_from_ary(*keys_);
      header->bits(vals_->bits() - 1);
      header->key_bytes(keys_->size_bytes());
      header->value_bytes(vals_->size_bytes());
      header->write(os);
    }
    keys_->write(os);
    vals_->write(os);
  }

  void done() {
    done_threads_ += 1;
    while(handle_full_ary() == OK);
  }

  mer_array& keys() { return *keys_; }
  val_array& vals() { return *vals_; }

private:
  status handle_full_ary() {
    bool serial_thread = size_barrier_.wait();
    if(done_threads_ >= nb_threads_) // All done?
      return DONE;

    if(serial_thread) {
      new_keys_ = 0;
      new_vals_ = 0;
      try {
        new_keys_ = new mer_array(keys_->size() * 2, keys_->key_len(), keys_->val_len(),
                                  keys_->max_reprobe(), keys_->reprobes());
        new_vals_ = new val_array(vals_->bits(), new_keys_->size());
      } catch(...) {
        delete new_keys_;
        delete new_vals_;
        new_keys_ = 0;
        new_vals_ = 0;
      }
      size_thid_ = 0;
    }
    size_barrier_.wait();
    mer_array* new_keys = *(mer_array* volatile *)&new_keys_;
    val_array* new_vals = *(val_array* volatile *)&new_vals_;
    if(!new_keys || !new_vals) {
      size_barrier_.wait();
      return FULL;
    }

    uint16_t thid = (size_thid_ += 1) - 1;
    auto it = keys_->eager_slice(thid, nb_threads_);

    bool   is_new;
    size_t id;
    while(it.next()) {
      new_keys->set(it.key(), &is_new, &id);
      auto entry = (*new_vals)[id];
      entry.get();
      uint64_t ov = (*vals_)[it.id()].get();
      entry.set(ov);
    }
    size_barrier_.wait();
    if(serial_thread) {
      delete keys_;
      delete vals_;
      keys_ = new_keys;
      vals_ = new_vals;
    }
    size_barrier_.wait();

    return OK;
  }
};

class suck_in_file {
public:
  suck_in_file(const char* path) : base_(0) { read_in(path); }
  suck_in_file(int fd) : base_(0) { read_in(fd); }
  ~suck_in_file() { }

  char* base() const { return base_; }
  define_error_class(ErrorReading);

protected:
  void read_in(int fd, const char* path = "<unknown>") {
    delete[] base_;
    struct stat buf;
    if(fstat(fd, &buf) < 0)
      throw ErrorReading(err::msg() << "Can't stat file '" << path << "'" << err::no);
    base_ = new (std::nothrow) char[buf.st_size];
    if(!base_)
      throw ErrorReading(err::msg() << "Not enough memory to read in file '" << path << "'" << err::no);
    ssize_t sucked = 0;
    while(sucked < buf.st_size) {
      ssize_t s = read(fd, base_ + sucked, buf.st_size - sucked);
      if(s == -1)
        throw ErrorReading(err::msg() << "Failed to read in file '" << path << "'");
      sucked += s;
    }
  }
  void read_in(const char* path) {
    int fd = open(path, O_RDONLY);
    if(fd < 0)
      throw ErrorReading(err::msg() << "Can't open file '" << path << "'" << err::no);
    read_in(fd, path);
    close(fd);
  }

private:
  char* base_;
};

class map_or_read_file {
  std::unique_ptr<const jellyfish::mapped_file> mapped;
  std::unique_ptr<const suck_in_file>           sucked;

public:
  map_or_read_file(const char* filename, bool no_map) {
    if(no_map) {
      sucked.reset(new suck_in_file(filename));
    } else {
      mapped.reset(new jellyfish::mapped_file(filename));
      vlog << "Mer database bogus checksum: " << (int)mapped->load();
    }
  }

  char* base() {
    if(mapped)
      return mapped->base();
    else
      return sucked->base();
  }
};


class database_query {
  const database_header        header_;
  map_or_read_file             file_;
  const mer_array_raw          keys_;
  const val_array_raw          vals_;

  static database_header parse_header(const char* filename) {
    std::ifstream file(filename);
    if(!file.good())
      throw std::runtime_error(err::msg() << "Can't open '" << filename << "' for reading");
    database_header res;
    if(!res.read(file))
      throw std::runtime_error(err::msg() << "Can't parse header of file '" << filename << "'");
    if(!res.check_format())
      std::runtime_error(err::msg() << "Wrong type '" << res.format() << "' for file '" << filename << "'");
    return res;
  }

public:
  database_query(const char* filename, bool map = false) :
  header_(parse_header(filename)),
  file_(filename, map),
  keys_(file_.base() + header_.offset(), header_.key_bytes(),
        header_.size(), header_.key_len(), header_.val_len(),
        header_.max_reprobe(), header_.matrix()),
  vals_(file_.base() + header_.offset() + header_.key_bytes(), header_.value_bytes(),
        header_.bits() + 1, header_.size())
  { }

  const database_header& header() const { return header_; }
  const mer_array_raw& keys() const { return keys_; }
  const val_array_raw& vals() const { return vals_; }

  std::pair<uint64_t, int> operator[](const mer_dna& m) const {
    std::pair<uint64_t, int> res = { 0, 0 };
    size_t                   id  = 0;
    if(keys_.get_key_id(m, &id)) {
      uint64_t v = vals_[id];
      res.first  = v >> 1;
      res.second = v & 0x1;
    }
    return res;
  }

  // Get value of m in the high quality database
  uint64_t get_val(const mer_dna& m) const {
    auto v = operator[](m);
    return v.second ? v.first : 0;
  }

  // Get all alternatives at the best level
  template<typename mer_type>
  int get_best_alternatives(mer_type& m, uint64_t counts[4], int& ucode, int& level) const {
    mer_dna tmp_mer;
    int     count = 0;
    memset(counts, '\0', sizeof(uint64_t) * 4);
    level = 0;
    int ori_code = m.code(0);

    for(int i = 0; i < 4; ++i) {
      m.replace(0, i);
      auto v = operator[](m.canonical());
      if(v.first > 0) {
        if(v.second >= level) {
          if(v.second > level && count > 0) {
            for(int j = 0; j < i; ++j)
              counts[j] = 0;
            count = 0;
          }
          counts[i] = v.first;
          ucode     = i;
          level     = v.second;
          ++count;
        }
      }
    }
    m.replace(0, ori_code); // Reset m to original value
    return count;
  }

  class const_iterator :
    public std::iterator<std::forward_iterator_tag, std::pair<const mer_dna*, std::pair<uint64_t, int> > > {
    mer_array_raw::const_iterator it_;
    const val_array_raw&          vals_;
    value_type                    content_;
  public:
    const_iterator(const mer_array_raw::const_iterator it, const val_array_raw& vals) :
      it_(it), vals_(vals)
    { }

    bool operator==(const const_iterator& rhs) const { return it_ == rhs.it_; }
    bool operator!=(const const_iterator& rhs) const { return it_ != rhs.it_; }
    const_iterator& operator++() { ++it_; return *this; }
    const_iterator operator++(int) {
      const_iterator res(*this);
      ++*this;
      return res;
    }

    const value_type& operator*() {
      content_.first         = &it_.key();
      uint64_t v             = vals_[it_.id()];
      content_.second.first  = v >> 1;
      content_.second.second = v & 0x1;
      return content_;
    }
    const value_type* operator->() { return &this->operator*(); }
  };

  const_iterator begin() const { return const_iterator(keys_.begin(), vals_); }
  const_iterator end() const { return const_iterator(keys_.end(), vals_); }
};

#endif /* __QUORUM_MER_DATABASE_HPP__ */
