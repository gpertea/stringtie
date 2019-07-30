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


/*  This file is part of k_unitig.

    k_unitig is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    k_unitig is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with k_unitig.  If not, see <http://www.gnu.org/licenses/>.
*/


/* Will eventually be replaced by direct indexing array withing
 * Jellyfish distribution. Supports only one operation: set_to_max
 */

#ifndef __ALIGNED_SIMPLE_ARRAY_HPP__
#define __ALIGNED_SIMPLE_ARRAY_HPP__

namespace simple_array {
  template <typename word, typename atomic_t, typename mem_block_t>
  class array {
  public:
    array(size_t _size) :
      size(_size), 
      mem_block(sizeof(word) * size),
      data((word *)mem_block.get_ptr()) {}

    inline word cas(const size_t &id, const word &old_value, const word &new_value) {
      return atomic.cas(data + id, old_value, new_value);
    }

    /**
     * Set entry at position id to val if the entry is empty (i.e. ==
     * 0). Return the value that was at id before the call.
     */
    word set_if_empty(const size_t id, const word &val) {
      return atomic.cas(data + id, (word)0, val);
    }

    /**
     * Set entry at position id to the max of the value at position
     * id and val. It returns the position at id before the call.
     */
    word set_to_max(const size_t id, const word &val) {
      word ov = data[id];
      while(val > ov) {
        word nov = cas(id, ov, val);
        if(nov == ov)
          break;
        ov = nov;
      }
      return ov;
    }
  private:
    size_t         size;
    mem_block_t    mem_block;
    volatile word *data;
    atomic_t       atomic;
  };
}

#endif
