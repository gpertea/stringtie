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


#ifndef __BLOOM_HASH_HPP__
#define __BLOOM_HASH_HPP__

#include <src/MurmurHash3.h>

/* Hash pairs
 */
template<typename Key>
class hash_pair { };

template <>
class hash_pair<const char*> {
public:
  void operator()(const char* const key, uint64_t *hashes) const {
    MurmurHash3_x64_128(key, strlen(key), 0, hashes);
  }
};

#endif // __BLOOM_HASH_HPP__
