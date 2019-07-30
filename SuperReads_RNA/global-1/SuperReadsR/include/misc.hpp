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


#ifndef __MISC_HPP__
#define __MISC_HPP__

//#define _GNU_SOURCE
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>

// IO manipulator for substrings
class substr {

public:
  const char * const _str;
  size_t const       _len;
  substr(const char *str, size_t len) :
    _str(str), _len(len) {}
  substr(const char *str, const char *end) :
    _str(str), _len(end > str ? end - str : 0) {}
};
inline bool is_base(const char c) {
  switch(c) {
  case 'A': case 'C': case 'G': case 'T':
  case 'a': case 'c': case 'g': case 't':
    return true;

  default:
    return false;
  }
}
inline bool not_base(const char c) { return !is_base(c); }
inline std::ostream & operator<<(std::ostream &os, const substr &ss) {
  assert(std::count_if(ss._str, ss._str+ss._len, not_base) == 0);
  return os.write(ss._str, ss._len);
}

template<typename T>
int getFldsFromLine(char *line, T &res, const char* sep = " \t\n") {
  char *saveptr;
  res.clear();

  char *tok = strtok_r(line, sep, &saveptr);
  while(tok) {
    res.push_back(tok);
    tok = strtok_r(0, sep, &saveptr);
  }
  return res.size();
}

template<typename T>
int appendFldsFromLine(char *line, T &res, const char* sep = " \t\n") {
  char *saveptr;
  int numFlds = 0;

  char *tok = strtok_r(line, sep, &saveptr);
  while(tok) {
    res.push_back(tok);
    ++numFlds;
    tok = strtok_r(0, sep, &saveptr);
  }
  return numFlds;
}

#include <signal.h>

#define BREAKPOINT raise(SIGINT);

template<long int n>
struct ConstFloorLog2 {
  static const int val = ConstFloorLog2<n / 2>::val + 1;
};
template<>
struct ConstFloorLog2<1> {
  static const int val = 0;
};
// Return length random bits
inline uint64_t random_bits(int length) {
  uint64_t res = 0;
  for(int i = 0; i < length; i += ConstFloorLog2<RAND_MAX>::val) {
    res ^= (uint64_t)random() << i;
  }
  return res & ((uint64_t)-1 >> (8 * sizeof(uint64_t) - length));
}


#endif
