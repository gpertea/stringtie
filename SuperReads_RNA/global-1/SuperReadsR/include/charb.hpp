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


#ifndef __CHARB_HPP__
#define __CHARB_HPP__

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <stdexcept>
#include <string>
#include <cstring>
#include <istream>
#include <assert.h>

#include <exp_buffer.hpp>

/** An almost drop in replacement for `char*` with safe memory
    management.

    Many str* and IO function of the C standard library are overloaded
    to work with `charb`. A `charb` is cast automatically to a `char*`
    or `const char*` when necessary. Hence, the function `strlen(const
    char*)` can be called on a `charb`.

   The following calls are know not to be overloaded: strfmon,
   strftime, strxfrm.
*/

// For the friend business to work
template<typename R> class basic_charb;
template<typename R> char *fgets(basic_charb<R>& b, FILE* stream, char* cptr);
template<typename R> int vsprintf(basic_charb<R>& b, char*, const char* format, va_list ap);
template<typename R> char *strcat(basic_charb<R>& b, const char* src);
template<typename R> char *strcat(basic_charb<R>& b, const basic_charb<R>& src);
template<typename R> std::istream& getline(std::istream& is, basic_charb<R>& b, char delim, char *cptr);

/** Basic base charb class. This implement a 0
 * terminated growable array of char with overloaded functions for
 * many/most of the str* and IO function of the C standard library. It
 * is intended to replace string manipulation with minimal
 * modification to a C program. A `charb` is cast automatically to a
 * `char*` or `const char*` when necessary. Hence, the function
 * `strlen(const char*)` can be called on a `charb`.

 *
 * The base position of a charb may change due to reallocation. Assume
 * that any call to a standard library function that write into a
 * charb invalidates any pointers into that charb.
 *
 * See charb.hpp for more information on the overloaded functions.
 */
template<typename R>
class basic_charb : public ExpBuffer<char, R> {
  typedef ExpBuffer<char, R> super;
public:
  basic_charb() : super(1) {
    *super::ptr_ = '\0';
  }
  explicit basic_charb(size_t s) : super(s + 1) {
    *super::ptr_ = '\0';
  }
  basic_charb(const basic_charb &rhs) : super(rhs.base_, rhs.len() + 1) {
    *--super::ptr_ = '\0';
  }
  basic_charb(basic_charb&& rhs) : super(std::move(rhs)) { }
  basic_charb(const char *str) : super(str, strlen(str) + 1) {
    *--super::ptr_ = '\0';
  }
  basic_charb(const char *str, size_t str_len) : super(str, str_len + 1) {
    *--super::ptr_ = '\0';
  }
  basic_charb(const std::string &s) : super(s.c_str(), s.size() + 1) {
    *--super::ptr_ = '\0';
  }
  virtual ~basic_charb() { }

  basic_charb &operator=(const char *rhs) {
    size_t rhs_len = strlen(rhs);
    super::reserve(rhs_len + 1);
    strcpy(super::base_, rhs);
    super::ptr_ = super::base_ + rhs_len;
    return *this;
  }

  basic_charb& operator=(basic_charb&& rhs) {
    this->swap(rhs);
    return *this;
  }

  basic_charb& operator=(const basic_charb& rhs) {
    size_t rhs_len = rhs.len();
    super::reserve(rhs_len + 1);
    strcpy(super::base_, rhs.base_);
    super::ptr_ = super::base_ + rhs_len;
    return *this;
  }
  /** Length of string. The length of the string is updated by the
   * many methods (`fgets`, `chomp`, etc.). The length returned can in
   * some cases be different that what is returned by `strlen(b)`, for
   * example in the case that the charb was modified using other
   * methods (e.g. `strtok`) or directly by using pointers.
   */
  size_t len() const { return super::size(); }
  /** Set the length to be 0. */
  void clear() {
    super::clear();
    if(super::ptr_)
      *super::ptr_ = '\0';
  }
  /** Remove space characters (as defined by `isspace()`) from the end. */
  void chomp() {
    while(super::ptr_ > super::base_ && isspace(*(super::ptr_ - 1)))
      --super::ptr_;
    *super::ptr_ = '\0';
  }
  /** Truncate the string to the given length. The len() method will
      return the length passed. */
  void truncate(size_t s) {
    super::reserve(s + 1);
    super::ptr_  = super::base_ + s;
    *super::ptr_ = '\0';
  }

  friend char *fgets <> (basic_charb<R> &b, FILE *stream, char *cptr);
  friend int vsprintf <> (basic_charb<R> &b, char* start, const char *format, va_list ap);
  friend char *strcat <> (basic_charb<R> &b, const char *src);
  friend char *strcat <> (basic_charb<R> &b, const basic_charb<R>& src);
  friend std::istream& getline <> (std::istream& is, basic_charb<R> &b, char delim, char *cptr);
};

template<typename T>
struct c_string {
  typedef basic_charb<reallocator<T> > malloc;
  typedef basic_charb<remaper<T> > remap;
};
//typedef c_string<char>::malloc charb;
/** Charb specialized for `char` */
typedef basic_charb<reallocator<char> > charb;

/** Input of line for char buffer. Expand the size
    of the buffer if the line does not fit.

    @param b The charb to write to.
    @param stream The input stream.
    @param cptr The pointer inside the charb (not checked) to write to.
    @return On success, a pointer to the beginning of the charb. NULL on error or end of file while no characters have been read.
 */
template<typename R>
char *fgets(basic_charb<R> &b, FILE *stream, char *cptr) {
  long  pos   = ftell(stream);
  long  npos  = pos;
  char *start = cptr;

  if(b.capacity() <= 1) {
    b.reserve(20);
    start = cptr = b.base();
  }

  while(true) {
    char *res = fgets(cptr, b.capacity() - (cptr - b.base_), stream);
    if(!res)
      break;
    size_t char_read;
    if(pos == -1) {
      char_read = strlen(res);
    } else {
      npos = ftell(stream);
      char_read = npos - pos;
      pos = npos;
    }
    cptr      += char_read;
    if(cptr == b.base_)
      asm("int3;");
    if(cptr < b.end_ - 1 || *(cptr - 1) == '\n')
      break;
    size_t off  = cptr  - b.base_;
    size_t soff = start - b.base_;
    b.enlarge();
    cptr  = b.base_ + off;
    start = b.base_ + soff;
  }

  if(cptr == start)
    return 0;
  assert(cptr != NULL);
  b.ptr_ = cptr;
  return start;
}

/** Input of line for char buffer. The previous content of the charb is overwritten.

    @param b The charb to write to
    @param stream The input stream
    @return On success, a pointer to the beginning of the charb. NULL on error or end of file while no characters have been read.
*/
template<typename R>
char *fgets(basic_charb<R> &b, FILE *stream) { return fgets(b, stream, b.base()); }

/** Input of line from stdin. Equivalent to `fgets(b, stdin)`.

    @param b The charb to write to
*/
template<typename R>
char *gets(basic_charb<R> &b) { return fgets(b, stdin); }

/** Similar to `fgets`, but append to the `charb`. I.e. write starts
    after `length()` characters.

    @param b The charb to append to
    @param stream The input stream
*/
template<typename R>
char *fgets_append(basic_charb<R> &b, FILE *stream) { return fgets(b, stream, b.ptr()); }

/** Backward compatible `fgets` for char buffer.

    @param b The charb to write to
    @param size Ignored. Present for backward compatibility
    @param stream The input stream
    @return
 */
template<typename T, typename R>
char *fgets(basic_charb<R> &b, T size, FILE *stream) { return fgets(b, stream); }

/** Input one line into a char buffer from a standard stream. This is
    equivalent to fgets, except for the return value that conforms to
    getline(3).

    @param b The charb to write to
    @param stream The input stream
    @param start Pointer into b to write to (Caution: no check)
    @return On success, the number of characters read, including the newline character, but not includeing the terminating null byte. On failure, return -1 (including end-of-file condition).
 */
template<typename R>
ssize_t getline(basic_charb<R> &b, FILE *stream, char* start) {
  const char* ret = fgets(b, stream, b.base());
  if(!ret)
    return -1;
  return b.size();
}

/** Input one line into a char buffer from a standard stream. This is
    equivalent to fgets, except for the return value that conforms to
    getline(3).

    @param b The charb to write to
    @param stream The input stream
    @return On success, the number of characters read, including the newline character, but not includeing the terminating null byte. On failure, return -1 (including end-of-file condition).
 */
template<typename R>
ssize_t getline(basic_charb<R> &b, FILE *stream) {
  const char* ret = fgets(b, stream, b.base());
  if(!ret)
    return -1;
  return b.size();
}

/** Append one line to a char buffer from a standard stream. This is
    equivalent to fgets_append, except for the return value that conforms to
    getline(3).

    @param b The charb to append to
    @param stream The input stream
    @return On success, the number of characters read, including the newline character, but not includeing the terminating null byte. On failure, return -1 (including end-of-file condition).
 */
template<typename R>
ssize_t getline_append(basic_charb<R> &b, FILE* stream) {
  size_t  osize = b.size();
  ssize_t ret   = getline(b, stream, b.ptr());
  if(ret < 0)
    return ret;
  return ret - osize;
 }

/** Input one line into a char buffer from a C++ stream. Similar to
    the global function getline(istream& is, string& str). Unlike the
    C version of fgets and getline, the delimiter character is not
    included in the output string.

    @param is The input stream
    @param b The charb to write to
    @param delim The end of line delimiter
    @param cptr Pointer into the charb to write to (Caution: no check made)
    @return The input stream
 */
template<typename R>
std::istream& getline(std::istream &is, basic_charb<R> &b, char delim, char *cptr) {
  if(b.capacity() <= 1) {
    b.reserve(20);
    cptr = b.base();
  }

  while(true) {
    is.getline(cptr, b.capacity() - (cptr - b.base_), delim);
    if(is.bad()) // Bad, we quit
      break;
    cptr += strlen(cptr);
    if(is.eof())  // Eof, we quit
      break;
    if(!is.fail()) // Found delim, done
      break;
    // Need to extend buffer
    ssize_t off = cptr - b.base_;
    b.enlarge();
    cptr = b.base_ + off;
    is.clear();
  }
  b.ptr_ = cptr;
  return is;
}

/** Input one line into a char buffer from a C++ stream. Similar to
    the global function getline(istream& is, string& str). Unlike the
    C version of fgets and getline, the delimiter character is not
    included in the output string.

    @param is The input stream
    @param b The charb to write to
    @return The input stream
 */

template<typename R>
std::istream& getline(std::istream& is, basic_charb<R>& b) { return getline(is, b, '\n', b.base()); }

/** Input one line into a char buffer from a C++ stream. Similar to
    the global function getline(istream& is, string& str). Unlike the
    C version of fgets and getline, the delimiter character is not
    included in the output string.

    @param is The input stream
    @param b The charb to write to
    @param delim The end of line delimiter
    @return The input stream
 */
template<typename R>
std::istream& getline(std::istream& is, basic_charb<R>& b, char delim) { return getline(is, b, delim, b.base_); }

/** Append one line into a char buffer from a C++ stream. Similar to
    the global function getline(istream& is, string& str) except the
    line is appended instead of overwriting. Unlike the C version of
    fgets and getline, the delimiter character is not included in the
    output string.

    @param is The input stream
    @param b The charb to write to
    @return The input stream
 */
template<typename R>
std::istream& getline_append(std::istream& is, basic_charb<R>& b) { return getline(is, b, '\n', b.ptr()); }

/** Append one line into a char buffer from a C++ stream. Similar to
    the global function getline(istream& is, string& str) except the
    line is appended instead of overwriting. Unlike the C version of
    fgets and getline, the delimiter character is not included in the
    output string.

    @param is The input stream
    @param b The charb to write to
    @param delim The end of line delimiter
    @return The input stream
 */
template<typename R>
std::istream& getline_append(std::istream& is, basic_charb<R>& b, char delim) { return getline(is, b, delim, b.ptr()); }


/** Formatted output conversion. The charb grows as needed.

    @param b The charb to write to
    @param format The format string
    @return The number of characters written or a negative number on error
 */
template<typename R>
int sprintf(basic_charb<R> &b, const char* format, ...) __attribute__ ((format (printf, 2, 3)));
template<typename R>
int sprintf(basic_charb<R> &b, const char *format, ...) {
  va_list ap;

  va_start(ap, format);
  int res = vsprintf(b, (char*)b, format, ap);
  va_end(ap);

  return res;
}

/** Formatted output conversion. For backward compatibility (the
    length is ignored). The charb grows as needed.

    @param b The charb to write to
    @param size Ignored
    @param format The format string
    @return The number of characters written or a negative number on error
 */
template<typename T, typename R>
int snprintf(basic_charb<R> &b, T size, const char *format, ...) __attribute__ ((format (printf, 3, 4)));
template<typename T, typename R>
int snprintf(basic_charb<R> &b, T size, const char *format, ...) {
  va_list ap;
  va_start(ap, format);
  int res = vsprintf(b, (char*)b, format, ap);
  va_end(ap);

  return res;
}

/** Formatted output conversion. For backward compatibility (the
    length is ignored). The charb grows as needed.

    @param b The charb to write to
    @param size Ignored
    @param format The format string
    @param ap The variable argument list
    @return The number of characters written or a negative number on error
 */
template<typename T, typename R>
inline int vsnprintf(basic_charb<R> &b, T size, const char *format, va_list ap) __attribute__ ((format (printf, 3, 0)));
template<typename T, typename R>
inline int vsnprintf(basic_charb<R> &b, T size, const char *format, va_list ap) {
  return vsprintf(b, (char*)b, format, ap);
}

/** Formatted output conversion. The charb grows as needed.

    @param b The charb to write to
    @param format The format string
    @param ap The variable argument list
    @return The number of characters written or a negative number on error
 */
template<typename R>
inline int vsprintf(basic_charb<R>& b, const char* format, va_list _ap) __attribute__ ((format (printf, 2, 0)));
template<typename R>
inline int vsprintf(basic_charb<R>& b, const char* format, va_list _ap) {
  return vsprintf(b, (char*)b, format, _ap);
}

/** Formatted output conversion. The text is appended to the charb
    instead of overwriting. The charb grows as needed.

    @param b The charb to write to
    @param format The format string
    @return The number of characters written or a negative number on error
 */
template<typename R>
int sprintf_append(basic_charb<R>& b, const char* format, ...) __attribute__ ((format (printf, 2, 3)));
template<typename R>
int sprintf_append(basic_charb<R>& b, const char* format, ...) {
  va_list ap;

  va_start(ap, format);
  int res = vsprintf(b, b.ptr(), format, ap);
  va_end(ap);

  return res;
}

/** Formatted output conversion. The text is appended to the charb
    instead of overwriting. The charb grows as needed.

    @param b The charb to write to
    @param format The format string
    @param ap The variable argument list
    @return The number of characters written or a negative number on error
 */
template<typename R>
inline int vsprintf_append(basic_charb<R>& b, const char* format, va_list _ap) __attribute__ ((format (printf, 2, 0)));
template<typename R>
inline int vsprintf_append(basic_charb<R>& b, const char* format, va_list _ap) {
  return vsprintf(b, b.ptr(), format, _ap);
}

/** Formatted output conversion. The charb grows as needed.

    @param b The charb to write to
    @param start Where to write into the charb (Caution: no check made)
    @param format The format string
    @param ap The variable argument list
    @return The number of characters written or a negative number on error
 */
template<typename R>
int vsprintf(basic_charb<R> &b, char* start, const char* format, va_list _ap) __attribute__ ((format (printf, 3, 0)));
template<typename R>
int vsprintf(basic_charb<R> &b, char* start, const char* format, va_list _ap) {
  int res = 0;
  while(true) {
    va_list ap;
    va_copy(ap, _ap);
    size_t remain = b.end_ - start;
    res = vsnprintf(start, remain, format, ap);
    va_end(ap);
    if(res < 0)
      return res;
    if((size_t)res < remain)
      break;
    size_t offset = start - b.base_;
    b.reserve(offset + res + 1);
    start = b.base_ + offset;
  }
  b.ptr_ = start + res;
  return res;
}

/** Concatenate two strings. The charb grows as needed.

    @param b The charb to append to
    @param src The string to append
    @return A pointer to the beginning of b
 */
template<typename R>
char *strcat(basic_charb<R> &b, const char *src) {
  size_t b_len = (char*)b ? b.len() : 0;
  size_t src_len = strlen(src);
  b.reserve(b_len + src_len + 1);
  strncpy((char*)b + b_len, src, src_len + 1);
  b.ptr_ = b.base_ + b_len + src_len;
  return b;
}

/** Concatenate two charb strings. The first charb grows as needed.

    @param b The charb to append to
    @param src The string to append
    @return A pointer to the beginning of b
 */
template<typename R>
char *strcat(basic_charb<R> &b, const basic_charb<R>& src) {
  size_t b_len = (char*)b ? b.len() : 0;
  size_t src_len = (char*)src ? src.len() : 0;
  if(src_len != 0) {
    b.reserve(b_len + src_len + 1);
    strncpy((char*)b + b_len, src, src_len + 1);
    b.ptr_ = b.base_ + b_len + src_len;
  }
  return b;
}

/** Concatenate two strings. For backward compatibility (the size is
    ignored). The charb grows as needed.

    @param b The charb to append to
    @param size Ignored
    @param src The string to append
    @return A pointer to the beginning of b
 */
template<typename T, typename R>
char *strncat(basic_charb<R> &b, T size, const char *src) {
  return strcat(b, src);
}

/** Copy a string. Copy a string into a charb, which grows as needed.

    @param b The charb to write to
    @param src The input string
    @return A pointer to the beginning of b
 */
template<typename R>
char *strcpy(basic_charb<R> &b, const char *src) {
  b = src;
  return b;
}

/** Copy a string. Copy a string into a charb, which grows as
    needed. For backward compatibility (the length is ignored).

    @param b The charb to write to
    @param src The input string
    @param n Ignored
    @return A pointer to the beginning of b
 */
template<typename T, typename R>
char *strncpy(basic_charb<R> &b, const char *src, T n) {
  b = src;
  return b;
}

/** Return string describing error number. The behavior is slightly
  different than the original C version. In the XSI-compliant mode, 0
  is always returned and the charb is grown to accomodate the error
  message. In the GNU-specific version, the message is always copied
  into the charb and (char*)b is returned (i.e. it never returns a
  pointer to an immutable static string). See man strerror_r(3).

  @param errnum Error number
  @param b Charb to write to
  @return Always return 0
 */
#if (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && ! _GNU_SOURCE
template<typename R>
int strerror_r(int errnum, basic_charb<R> &b) {
  while(true) {
    int res = strerror_r(errnum, (char*)b, b.capacity());
    if(res == 0)
      break;
    if(errno != ERANGE)
      break;
    b.enlarge();
  }
  return res
}
/** Return string describing error number. Identical to the two
 argument version (the buflen is ignored).

  @param errnum Error number
  @param b Charb to write to
  @param buflen Ignored
  @return Always return 0
 */
template<typename T, typename R>
int strerror_r(int errnum, basic_charb<R> &b, T buflen) {
  return strerror_r(errnum, b);
}
#else
/** Return string describing error number. The behavior is slightly
  different than the original C version. In the XSI-compliant mode, 0
  is always returned and the charb is grown to accomodate the error
  message. In the GNU-specific version, the message is always copied
  into the charb and (char*)b is returned (i.e. it never returns a
  pointer to an immutable static string). See man strerror_r(3).

  @param errnum Error number
  @param b Charb to write to
  @return A pointer to the beginning of b
 */
template<typename R>
char *strerror_r(int errnum, basic_charb<R> &b) {
  char *res = strerror_r(errnum, (char*)b, b.capacity());
  if(res != (char*)b)
    strcpy(b, res);
  return b;
}

/** Return string describing error number. Identical to the two
 argument version (the buflen is ignored).

  @param errnum Error number
  @param b Charb to write to
  @param buflen Ignored
  @return A pointer to the beginning of b
 */
template<typename T, typename R>
char *strerror_r(int errnum, basic_charb<R> &b, T buflen) {
  return strerror_r(errnum, b);
}
#endif
#endif /* __CHARB_HPP__ */
