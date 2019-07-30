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


#ifndef __GCC_BUILTINS_HPP__
#define __GCC_BUILTINS_HPP__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// TODO: check in configure the existence of hardware version
template<typename T>
int ctz(T x);

template<>
int ctz<unsigned int>(unsigned int x);
template<>
int ctz<unsigned long>(unsigned long x);
template<>
int ctz<unsigned long long>(unsigned long long x);

#ifdef HAVE_BUILTIN_PREFETCH
inline void prefetch_read_no(const void* addr) { __builtin_prefetch(addr, 0, 0); }
inline void prefetch_read_low(const void* addr) { __builtin_prefetch(addr, 0, 1); }
inline void prefetch_read_med(const void* addr) { __builtin_prefetch(addr, 0, 2); }
inline void prefetch_read_high(const void* addr) { __builtin_prefetch(addr, 0, 3); }
inline void prefetch_write_no(const void* addr) { __builtin_prefetch(addr, 1, 0); }
inline void prefetch_write_low(const void* addr) { __builtin_prefetch(addr, 1, 1); }
inline void prefetch_write_med(const void* addr) { __builtin_prefetch(addr, 1, 2); }
inline void prefetch_write_high(const void* addr) { __builtin_prefetch(addr, 1, 3); }
#else // HAVE_BUILTIN_PREFETCH
inline void prefetch_read_no(const void* addr) { }
inline void prefetch_read_low(const void* addr) { }
inline void prefetch_read_med(const void* addr) { }
inline void prefetch_read_high(const void* addr) { }
inline void prefetch_write_no(const void* addr) { }
inline void prefetch_write_low(const void* addr) { }
inline void prefetch_write_med(const void* addr) { }
inline void prefetch_write_high(const void* addr) { }
#endif // HAVE_BUILTIN_PREFETCH


#endif /* __GCC_BUILTINS_HPP__ */
