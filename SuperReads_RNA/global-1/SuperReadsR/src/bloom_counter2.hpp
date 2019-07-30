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


#ifndef __BLOOM_COUNTER2_HPP__
#define __BLOOM_COUNTER2_HPP__

#include <math.h>
#include <limits.h>
#include <reallocators.hpp>
#include <jellyfish/divisor.hpp>
#include <jellyfish/atomic_field.hpp>
#include <jellyfish/compare_and_swap.hpp>
#include <src/bloom_hash.hpp>

using jflib::divisor64;

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <gcc_builtins.hpp>

#define LOG2    0.6931471805599453
#define LOG2_SQ 0.4804530139182014
/* Bloom counter with 3 values: 0, 1 or 2. It is thread safe and lock free.
 */
template<typename Key, typename HashPair = hash_pair<Key>, typename R=remaper<unsigned char> >
class bloom_counter2 {
  // The number of bits in the structure, previously known as m_, is
  // know stored as d_.d()
  const divisor64       d_;
  const unsigned long   k_;
  unsigned char * const data_;
  HashPair              hash_fns_;

  struct prefetch_info {
    size_t         boff;
    unsigned char* pos;
  };

  static size_t nb_bytes(size_t l) {
    return l / 5 + (l % 5 != 0);
  }

public:
  typedef Key key_type;

  bloom_counter2(double fp, size_t n) :
    d_(n * (size_t)lrint(-log(fp) / LOG2_SQ)),
    k_(lrint(-log(fp) / LOG2)),
    data_(R::realloc(0, 0, nb_bytes(d_.d()))),
    hash_fns_() { }

  bloom_counter2(size_t m, unsigned long k) :
    d_(m), k_(k),
    data_(R::realloc(0, 0, nb_bytes(d_.d()))),
    hash_fns_() { }

  bloom_counter2(size_t m, unsigned long k, unsigned char* ptr) :
    d_(m), k_(k), data_(ptr), hash_fns_() { }

  ~bloom_counter2() {
    R::realloc(data_, nb_bytes(d_.d()), 0);
  }

  void write_bits(std::ostream& out) {
    out.write((char*)data_, m() / 5 + (m() % 5 != 0));
  }

  // Number of hash functions
  unsigned long k() const { return k_; }
  // Size of bit vector
  size_t m() const { return d_.d() ; }

  // Insert key k. Returns previous value of k
  unsigned int insert(const Key &k) {
    uint64_t hashes[2];
    hash_fns_(k, hashes);
    return insert(hashes);
  }

  // Insert key with given hashes
  unsigned int insert(const uint64_t* hashes) {
    // Prefetch memory locations
    static_assert(std::is_pod<prefetch_info>::value, "Prefetch info is a POD");
    prefetch_info pinfo[k()];
    const size_t base = d_.remainder(hashes[0]);
    const size_t inc  = d_.remainder(hashes[1]);
    for(unsigned long i = 0; i < k_; ++i) {
      const size_t p   = d_.remainder(base + i * inc);
      const size_t off = p / 5;
      pinfo[i].boff    = p % 5;
      pinfo[i].pos     = data_ + off;
      prefetch_write_no(pinfo[i].pos);
    }

    // Insert element
    unsigned char res = 2;
    for(unsigned long i = 0; i < k_; ++i) {
      size_t        boff = pinfo[i].boff;
      unsigned char v    = jflib::a_load(pinfo[i].pos);
      
      while(true) {
        unsigned char w = v;
        switch(boff) {
        case 0:          break;
        case 1: w /= 3;  break;
        case 2: w /= 9;  break;
        case 3: w /= 27; break;
        case 4: w /= 81; break;
        }
        w = w % 3;
        if(w == 2) break;
        unsigned char nv = v;

        switch(boff) {
        case 0: nv += 1;  break;
        case 1: nv += 3;  break;
        case 2: nv += 9;  break;
        case 3: nv += 27; break;
        case 4: nv += 81; break;
        }
        if(jflib::cas(pinfo[i].pos, v, nv, &v)) {
          if(w < res)
            res = w;
          break;
        }
      }
    }
    return res;
  }

  unsigned int check(const Key &k) const {
    uint64_t hashes[2];
    hash_fns_(k, hashes);
    return check(hashes);
  }

  unsigned int check(uint64_t *hashes) const {
    // Prefetch memory locations
    static_assert(std::is_pod<prefetch_info>::value, "Prefetch info is a POD");
    prefetch_info pinfo[k()];
    const size_t base = d_.remainder(hashes[0]);
    const size_t inc  = d_.remainder(hashes[1]);
    for(unsigned long i = 0; i < k_; ++i) {
      const size_t p   = d_.remainder(base + i * inc);
      const size_t off = p / 5;
      pinfo[i].boff    = p % 5;
      pinfo[i].pos     = data_ + off;
      prefetch_read_no(pinfo[i].pos);
    }

    // Check element
    unsigned char res = 2;
    for(unsigned long i = 0; i < k_; ++i) {
      size_t        boff = pinfo[i].boff;
      unsigned char w    = jflib::a_load(pinfo[i].pos);

      switch(boff) {
      case 0:          break;
      case 1: w /= 3;  break;
      case 2: w /= 9;  break;
      case 3: w /= 27; break;
      case 4: w /= 81; break;
      }
      w = w % 3;
      if(w < res)
        res = w;
    }
    return res;
  }

  // Limited std::map interface compatibility
  class element_proxy {
    bloom_counter2& bc_;
    const Key&      k_;

  public:
    element_proxy(bloom_counter2& bc, const Key& k) : bc_(bc), k_(k) { }

    unsigned int operator++() { 
      unsigned int res = bc_.insert(k_);
      return res == 0 ? 1 : 2;
    }

    unsigned int operator++(int) { return bc_.insert(k_); }
    unsigned int operator*() const { return bc_.check(k_); }
    operator unsigned int() const { return bc_.check(k_); }
  };

  class const_element_proxy {
    const bloom_counter2& bc_;
    const Key&            k_;

  public:
    const_element_proxy(const bloom_counter2& bc, const Key& k) : bc_(bc), k_(k) { }

    unsigned int operator*() const { return bc_.check(k_); }
    operator unsigned int() const { return bc_.check(k_); }
  };

  element_proxy operator[](const Key& k) { return element_proxy(*this, k); }
  const_element_proxy operator[](const Key& k) const { return const_element_proxy(*this, k); }
};

#endif // __BLOOM_COUNTER2_HPP__
