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


#include <assert.h>
#include <limits.h>
#include <algorithm>
#include <vector>
#include <charb.hpp>
#include <jellyfish/err.hpp>
#include <src/fibonacci_coding.hpp>

namespace err = jellyfish::err;

// Return 1 plus the index of the first 1-bit in x.
#ifdef __builtin_ffsl
int ffsl(uint64_t x) { return __builtin_ffsl(x); }
#else
int ffsl(uint64_t x, int res = 1) {
  if(x == 0) return 0;
  if(x & 0x1) return res;
  return ffsl(x >> 1, res + 1);
}
#endif

class sr_name {
  const char * const data;
  sr_name() : data(0) { }

public:
  sr_name(const char *str) : data(from_str(str)) { }
  sr_name(const sr_name &rhs) : data(copy(rhs.data)) { }
  ~sr_name() { free((void*)data); }
  class iterator;
  iterator begin() const { return iterator(data); }
  iterator end() const { return iterator(); }

  const char* raw() const { return data; }
  void to_str(charb& str) { to_str(data, str); }
  
public:
  const char* copy(const char *str) {
    size_t len = strlen(str) + 1;
    char *res = (char*)malloc(len);
    if(!res)
      throw std::runtime_error(err::msg() << "Failed to allocate " << len << " bytes to copy sr_name");
    strcpy(res, str);
    return res;
  }
  // Iterator element. By definition begin has offset = 0 and end has
  // offset = (uint64_t)-1
  struct iterator_elt {
    uint64_t id;
    char     ori; // Either 'R' or 'F'
    uint64_t offset;
  };
  class iterator : public std::iterator<std::input_iterator_tag, iterator_elt> {
    iterator_elt elt;
    const char * data;
    uint64_t     chunk;
    uint32_t     chunk_used;
  public:
    iterator(const char* data_) : data(data_), chunk(0), chunk_used(0) {
      if(decode_fib(&elt.id)) {
        if(decode_ori(&elt.ori)) {
          elt.offset = 0;
          return;
        }
      }
      elt.offset = (uint64_t)-1;
    }

    friend class sr_name;
  public:
    iterator() : data(0), chunk(0), chunk_used(0) {
      elt.offset = -1;
    }
    // Default copy constructor
    iterator& operator=(const iterator rhs) {
      elt        = rhs.elt;
      data       = rhs.data;
      chunk      = rhs.chunk;
      chunk_used = rhs.chunk_used;
      return *this;
    }
    iterator& operator++() {
      if(!decode_fib(&elt.offset))
        goto failed;
      if(!decode_fib(&elt.id))
        goto failed;
      if(!decode_ori(&elt.ori))
        goto failed;
      return *this;

    failed:
      elt.offset = (uint64_t)-1;
      return *this;
    }
    iterator operator++(int) {
      iterator res(*this);
      ++*this;
      return res;
    }
    bool operator==(const iterator &rhs) const {
      if(rhs.elt.offset == (uint64_t)-1)
        return elt.offset == (uint64_t)-1;
      return data == rhs.data && chunk_used == rhs.chunk_used;
    }
    bool operator!=(const iterator &rhs) const { return !operator==(rhs); }
    const iterator_elt& operator*() const { return elt; }
    const iterator_elt* operator->() const { return &elt; }

  private:
    bool decode_fib(uint64_t *res) {
      while(true) {
        int fib_len  = ffsl(chunk & (chunk << 1));
        if(fib_len) {
          *res         = fibonacci::decode(chunk & mask(fib_len - 1, 0));
          chunk      >>= fib_len;
          chunk_used  -= fib_len;
          return true;
        }
        
        unsigned char c = *data++; 
        if(c == '\0')
          return false;
        c = ~c;
        uint64_t new_chunk = c;
        chunk       |= new_chunk << chunk_used;
        chunk_used  += CHAR_BIT;
      }
    }
    bool decode_ori(char *res) {
      if(chunk_used == 0) {
        unsigned char c = *data++;
        if(c == '\0')
          return false;
        c          = ~c;
        chunk      = (uint64_t)c;
        chunk_used = CHAR_BIT;
      }
      *res         = chunk & 0x1 ? 'R' : 'F';
      chunk      >>= 1;
      chunk_used  -= 1;
      return true;
    }
  };

  enum etype { FLAG, NUMBER };
  struct entry {
    etype    type;
    uint32_t nb_bits;
    uint64_t data;
  };
  typedef std::vector<entry> bin_array;
  static const char* from_str(const char* str) {
    bin_array fib_buffer;
    return from_str(str, fib_buffer);
  }
  static const char* from_str(const char* str, bin_array &fib_buffer) {
    bin_array ary;
    charb     tokens(str);
    char*     saveptr;
    char*     tok       = strtok_r(tokens, "_", &saveptr);
    size_t    total_len = 0;

    while(tok) {
      char* last_char = tok + strlen(tok) - 1;
      char ori = 0;
      switch(*last_char) {
      case 'R': ori = -1; break;
      case 'F': ori = 1; break;
      }
      if(ori)
        *last_char = '\0';

      entry e;
      e.type     = NUMBER;
      e.nb_bits  = fibonacci::encode((uint64_t)atol(tok), e.data);
      total_len += e.nb_bits + 1;
      ary.push_back(e);
      if(ori) {
        e.type    = FLAG;
        e.nb_bits = 1;
        e.data    = ori < 0 ? 0x1 : 0x0;
        ary.push_back(e);
        ++total_len;
      }

      tok = strtok_r(0, "_", &saveptr);
    }
    return from_bin_array(total_len, ary.begin(), ary.end());
  }

  static inline uint64_t mask(uint32_t len, uint32_t shift) {
    return (((uint64_t)1 << len) - 1) << shift;
  }

  template<typename Iterator>
  static const char *from_bin_array(size_t len, Iterator begin, Iterator end) {
    size_t    res_len    = len / CHAR_BIT + (len % CHAR_BIT != 0) + 1;
    char*     res        = (char*)malloc(res_len);
    uint64_t* chunk_ptr  = (uint64_t*)res;
    uint64_t  chunk      = 0;
    uint32_t  chunk_left = 64;

    if(!res)
      throw std::runtime_error(err::msg() << "Failed to allocate " << res_len << " bytes for sr_name");
    
    for(Iterator it = begin; it != end; ++it) {
      uint64_t append = 1;
      if(it->type == FLAG) {
        append = it->data;
      } else {
        chunk |= it->data << (64 - chunk_left);
        if(chunk_left <= it->nb_bits) {
          *chunk_ptr++ = ~chunk;
          chunk        = (it->data & mask(it->nb_bits - chunk_left, chunk_left)) >> chunk_left;
          chunk_left   = 64 - (it->nb_bits - chunk_left);
        } else {
          chunk_left -= it->nb_bits;
        }
        if(!chunk_left) {
          *chunk_ptr++ = ~chunk;
          chunk      = 0;
          chunk_left = 64;
        }
      }

      // Append last bit: orientation or end marker
      chunk |= append << (64 - chunk_left);
      chunk_left -= 1;
      if(!chunk_left) {
        *chunk_ptr++ = ~chunk;
        chunk      = 0;
        chunk_left = 64;
      }
    }

    // Deal with last chunk, which may be less than a 64 bit
    // word. Also append '\0' at the end.
    uint32_t nb_chars = (64 - chunk_left) / CHAR_BIT + ((64 - chunk_left) % CHAR_BIT != 0);
    if(nb_chars == 8) {
      *chunk_ptr++ = ~chunk;
      *(char*)chunk_ptr = '\0';
    } else {
      chunk |= (uint64_t)0xff << (CHAR_BIT * nb_chars);
      chunk = ~chunk;
      memcpy((char*)chunk_ptr, (char*)&chunk, nb_chars + 1);
    }
    
    return res;
  }

  static void to_str(const char* name, charb& str) {
    charb one_entry;
    iterator it = iterator(name);
    iterator end;
    if(it == end)
      return;
    sprintf(str, "%ld%c", it->id, it->ori);
    for(++it; it != end; ++it) {
      sprintf(one_entry, "_%ld_%ld%c", it->offset, it->id, it->ori);
      strcat(str, one_entry);
    }
  }
};
