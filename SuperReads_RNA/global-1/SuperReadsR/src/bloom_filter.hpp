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


#ifndef __BLOOM_FILTER_HPP__
#define __BLOOM_FILTER_HPP__

#include <math.h>
#include <vector>
#include <algorithm>
#include <jellyfish/divisor.hpp>
#include <reallocators.hpp>
#include <src/bloom_hash.hpp>
#include <gcc_builtins.hpp>

using jflib::divisor64;

/* Bloom filter using Kirsh & Mitzenmacher double hashing. I.e., only
   two hash functions are computed and the k functions have values (0
   <= i < k):

   hash_0 + i * hash_1
 */

/* Memory operator for serial access.
 */
template<typename T>
struct serial_access {
  static T fetch_and_or(T* ptr, T x) {
    T res = *ptr;
    *ptr |= x;
    return res;
  }

  static T fetch(T* ptr) {
    return *ptr;
  }
};

/* Memory operator for muti-threaded access.
 */
template<typename T>
struct mt_access {
  static T fetch_and_or(T* ptr, T x) { return __sync_fetch_and_or(ptr, x); }
  static T fetch(T* ptr) { return *(volatile T*)ptr; }
};

#define LOG2    0.6931471805599453
#define LOG2_SQ 0.4804530139182014
template<typename Key, typename HashPair = hash_pair<Key>, typename M=serial_access<unsigned int>, typename R=remaper<unsigned int>>
class bloom_filter {
  typedef typename R::element_type  element_type;
  typedef typename R::element_type* element_pointer;
  // The number of bits in the structure, previously known as m_, is
  // know stored as d_.d()
  const divisor64     d_;
  const unsigned long k_;       // Number of hashes
  element_pointer     data_;
  const HashPair      hash_fns_;
  M                   mem_access_;
  
  typedef std::vector<bool>::reference bit;

  struct prefetch_info {
    size_t          boff;
    element_pointer pos;
  };

  static const size_t elt_size = sizeof(element_type) * 8;
  static size_t nb_elements(size_t m) { return m / elt_size + (m % elt_size != 0); }

public:
  // BF with false positive rate of fp and estimated number of entries
  // of n.
  bloom_filter(double fp, size_t n) : 
    d_(n * (size_t)lrint(-log(fp) / LOG2_SQ)),
    k_(lrint(-log(fp) / LOG2)),
    data_(R::realloc(0, 0, nb_elements(d_.d()))),
    hash_fns_(), mem_access_()
  { }

  bloom_filter(size_t m, unsigned long k) :
    d_(m), k_(k),
    data_(R::realloc(0, 0, nb_elements(d_.d()))),
    hash_fns_(), mem_access_() { }

  ~bloom_filter() {
    R::realloc(data_, nb_elements(d_.d()), 0);
  }

  // Number of hash functions
  unsigned long k() const { return k_; }
  // Size of bit vector
  size_t m() const { return d_.d(); }

  // Std::set compatible insert method. There is no iterator for a
  // bloom filter, so return true for the first element of the
  // pair. The second element is true if the element was inserted and
  // was not in the set before.
  std::pair<bool, bool> insert(const Key &k) {
    return std::make_pair(true, !add(k));
  }

  // Insert a key. Return true if the element was already present.
  bool add(const Key& k) {
    uint64_t hashes[2];
    hash_fns_(k, hashes);
    return add(hashes);
  }
    
  // Insert key with given hashes
  bool add(const uint64_t *hashes) {
    // Prefetch memory locations
    static_assert(std::is_pod<prefetch_info>::value, "Prefetch info is a POD");
    prefetch_info pinfo[k()];
    const size_t base    = d_.remainder(hashes[0]);
    const size_t inc     = d_.remainder(hashes[1]);
    for(unsigned long i = 0; i < k_; ++i) {
      const size_t pos   = d_.remainder(base + i * inc);
      const size_t elt_i = pos / elt_size;
      pinfo[i].boff      = pos % elt_size;
      pinfo[i].pos       = data_ + elt_i;
      prefetch_write_no(pinfo[i].pos);
    }
    
    // Check if element present
    bool present = true;
    for(unsigned long i = 0; i < k_; ++i) {
      const element_type mask = ((element_type)1) << pinfo[i].boff;
      const element_type prev = mem_access_.fetch_and_or(pinfo[i].pos, mask);
      present                 = present && (prev & mask);
    }

    return present;
  }

  // Compute hashes of key k
  void hash(const Key &k, uint64_t *hashes) const { hash_fns_(k, hashes); }

  // True if k is a member of the set
  bool is_member(const Key &k) const {
    uint64_t hashes[2];
    hash_fns_(k, hashes);
    return is_member(hashes);
  }

  bool is_member(const uint64_t *hashes) const {
    // Prefetch memory locations
    prefetch_info pinfo[k()];
    const size_t base    = d_.remainder(hashes[0]);
    const size_t inc     = d_.remainder(hashes[1]);
    for(unsigned long i = 0; i < k_; ++i) {
      const size_t pos   = d_.remainder(base + i * inc);
      const size_t elt_i = pos / elt_size;
      pinfo[i].boff      = pos % elt_size;
      pinfo[i].pos       = data_ + elt_i;
      prefetch_read_no(pinfo[i].pos);
    }

    for(unsigned long i = 0; i < k_; ++i)
      if(!(mem_access_.fetch(pinfo[i].pos) & ((element_type)1 << pinfo[i].boff)))
        return false;
    return true;
  }
};

#endif // __BLOOM_FILTER_HPP__
