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


#ifndef __CHARBUF_HPP__
#define __CHARBUF_HPP__

#include <charb.hpp>

/** Streambuf backed by a charb. Similar to std::stringbuf (used by
 * std::ostringstream) but is based on a charb rather than a
 * string. The associated stream (charstream) has an overloaded
 * operator<< which allows to write content since the last rewind.
 */

class charbuf : public std::streambuf {
  charb buf;
public:
  charbuf(size_t s) : std::streambuf(), buf(s) { 
    setp(buf.base(), buf.base() + buf.capacity());
  }
  virtual ~charbuf() { }

  virtual std::streamsize xsputn(const char *s, std::streamsize n) {
    std::streamsize left = epptr() - pptr();
    if(left < n) {
      int off = pptr() - pbase();
      buf.enlarge();
      setp(buf.base(), buf.base() + buf.capacity());
      pbump(off);
    }
    memcpy(pptr(), s, n);
    pbump(n);
    return n;
  }

  virtual int overflow(int c) {
    size_t csize = pptr() - pbase();
    buf.enlarge();
    setp(buf.base(), buf.base() + buf.capacity());
    pbump(csize);
    if(c != EOF) {
      *pptr() = c;
      pbump(1);
    }
    return !EOF;
  }

  size_t size() const { return pptr() - pbase(); }
  char *str() const { return pbase(); }
  void rewind() {
    setp(buf.base(), buf.base() + buf.capacity()); // Any other way to reset pptr?
  }
};

class charstream : public std::ostream {
  charbuf *buf;
public:
  charstream(size_t s = 1024) : std::ostream(new charbuf(s)), buf(0) { 
    buf = (charbuf*)std::ostream::rdbuf();
  }
  virtual ~charstream() {
    std::ostream::rdbuf(0);
    delete buf;
    buf = 0;
  }
  size_t size() const { return buf->size(); }
  const char *str() const { return buf->str(); }
  void rewind() { buf->rewind(); }
};

std::ostream &operator<<(std::ostream &os, const charstream &cs) {
  return os.write(cs.str(), cs.size());
}

#endif /* __CHARBUF_HPP__ */
