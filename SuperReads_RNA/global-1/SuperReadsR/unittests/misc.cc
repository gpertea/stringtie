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


#include <unittests/misc.hpp>

void pdo(unsigned int n, void* (*f)(void*), void *data) {
  pthread_t tids[n];
  for(unsigned int i = 0; i < n; ++i) {
    int e = pthread_create(&tids[i], 0, f, data);
    if(e)
      throw std::runtime_error(strerror(e));
  }
  for(unsigned int i = 0; i < n; ++i) {
    int e = pthread_join(tids[i], 0);
    if(e)
      throw std::runtime_error(strerror(e));
  }
}
