/* SuperRead pipeline
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


#include <iostream>
#include <memory>
#include <algorithm>
#include <set>

#include <jellyfish/jellyfish.hpp>
#include <src/create_k_unitigs_common.hpp>

#include <gzip_stream.hpp>
#include <jflib/multiplexed_io.hpp>
#include <src/create_k_unitigs_large_k2_cmdline.hpp>

cmdline_parse args;

struct mer_counter {
  mer_hash hash_;

  mer_counter(size_t size, uint16_t key_len, uint16_t val_len, uint16_t nb_threads) :
    hash_(size, key_len, val_len, nb_threads)
  { }

  struct element {
    mer_hash&      hash_;
    const mer_dna& m_;

    element(mer_hash& hash, const mer_dna& m) : hash_(hash), m_(m) { }
    element& operator++() {
      hash_.add(m_, 1);
      return *this;
    }
  };

  struct const_element {
    const mer_hash& hash_;
    const mer_dna&  m_;

    const_element(const mer_hash& hash, const mer_dna& m) : hash_(hash), m_(m) { }
    operator uint64_t() const {
      uint64_t v = 0;
      hash_.ary()->get_val_for_key(m_, &v);
      return v;
    }
  };

  element operator[](const mer_dna& m) { return element(hash_, m); }
  const_element operator[](const mer_dna& m) const { return const_element(hash_, m); }
  void done() {
    hash_.done();
  }
};
typedef populate_mer_set<mer_counter, sequence_parser, read_mers> mer_populate;

struct mer_set {
  mer_hash::array set_;

  mer_set(size_t size, uint16_t key_len) :
    set_(size, key_len, 0, 126)
  { }
  std::pair<int, bool> insert(const mer_dna& m) {
    size_t id;
    bool   is_new;
    set_.set(m, &is_new, &id);
    return std::make_pair(0, is_new);
  }
};

struct used_type {
  template<typename T>
  std::pair<int, bool> insert(const T&) const { return std::make_pair(0, true); }
};

// Here, the parser is basically the hash itself, and a stream is an
// iterator into the hash.
struct mer_parser {
  const mer_hash::array& ary_;
  const int              nb_;

  mer_parser(const mer_hash::array& ary, int nb) : ary_(ary), nb_(nb) { }
  friend class mer_stream;
};
class mer_stream {
  mer_hash::eager_iterator it_;
  bool                         more_;
public:
  mer_stream(const mer_parser& p, int id) :
    it_(p.ary_.eager_slice(id, p.nb_)),
    more_(it_.next())
  { }
  mer_stream& operator++() {
    if(more_)
      more_ = it_.next();
    return *this;
  }
  const mer_dna& operator*() const { return it_.key(); }
  const mer_dna* operator->() const { return &it_.key(); }
  operator bool() const { return more_; }
};

typedef create_k_unitig<mer_counter, used_type, mer_set,
                        mer_parser, mer_stream,
                        cmdline_parse> unitiger_type;

std::ostream* open_output() {
  if(!args.output_given)
    return new std::ostream(std::cout.rdbuf());
  if(args.gzip_flag)
    return new gzipstream(args.output_arg);
  return new std::ofstream(args.output_arg);
}

int main(int argc, char *argv[])
{
  args.parse(argc, argv);
  mer_dna::k(args.mer_arg);
  if(!args.min_len_given)
    args.min_len_arg = args.mer_arg + 1;

  std::unique_ptr<mer_counter> kmers;
  {
    kmers.reset(new mer_counter(args.nb_mers_arg, 2 * args.mer_arg, args.val_arg, args.threads_arg));
    stream_manager  manager(args.input_arg.begin(), args.input_arg.end());
    sequence_parser parser(mer_dna::k(), manager.nb_streams(), 3 * args.threads_arg, 4096, manager);
    mer_populate    populate(args.mer_arg, *kmers, parser);
    populate.exec_join(args.threads_arg);
  }

  {
    std::unique_ptr<std::ostream> output_ostream(open_output());
    used_type                     used;
    mer_set                       ends(kmers->hash_.size(), kmers->hash_.key_len());
    mer_parser                    parser(*kmers->hash_.ary(), args.threads_arg);
    unitiger_type unitiger(*kmers, used, ends, args.threads_arg // 1
                           , parser,
                           *output_ostream, args);
    unitiger.exec_join(args.threads_arg // 1
                       );
  }

  return 0;
}
