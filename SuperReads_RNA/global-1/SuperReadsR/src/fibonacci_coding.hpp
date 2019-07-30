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


#ifndef __FIBONACCI_CODING_HPP__
#define __FIBONACCI_CODING_HPP__

#include <errno.h>
#include <stdint.h>
#include <iostream>
#include <algorithm>

struct fibonacci {
  static const uint64_t *fibs;

  template<typename T>
  static int encode(T x, T &res) {
    const uint64_t * f = fibs + 8 * sizeof(T) - 1;
    x++;
    if(x <= (T)0 || x >= *f) {
      errno = EINVAL; // Number is negative or too large
      return -1;
    }
    
    f             = std::upper_bound(fibs, f, x) - 1;
    int high_bit  = f - fibs;
    res           = (T)1 << high_bit;
    x            -= *f;
    while(x > 0) {
      f    = std::upper_bound(fibs, f, x) - 1;
      res |= (T)1 << (f - fibs);
      x   -= *f;
    }
    
    return high_bit + 1;
  }

  template<typename T>
  static T decode(T x) {
    T res = 0;
    for(const uint64_t *f = fibs; x; x >>= 1, ++f)
      if(x & 0x1)
        res += *f;
    return res - 1;
  } 
};

#endif
