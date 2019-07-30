/* Quorum
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

#ifndef _VERBOSE_LOG_HPP_
#define _VERBOSE_LOG_HPP_

#include <ctime>
#include <iostream>
#include <iomanip>

template<typename T>
class basic_verbose_log : public std::basic_ostringstream<T> {
  typedef std::ostringstream super;

public:
  static bool verbose;

  basic_verbose_log() : super() {
    time_t t = time(0);
    struct tm d;
    localtime_r(&t, &d);
    //    auto save = this->flags();
    const T f = this->fill('0');
    *this << "["
          << std::setw(4) << (1900 + d.tm_year) << "/"
          << std::setw(2) << (d.tm_mon + 1) << "/"
          << std::setw(2) << d.tm_mday << " "
          << std::setw(2) << d.tm_hour << ":"
          << std::setw(2) << d.tm_min << ":"
          << std::setw(2) << d.tm_sec
          << "] ";
    this->fill(f);
  }
  ~basic_verbose_log() {
    if(verbose)
      std::cerr << super::str() << std::endl;
  }

  std::ostream& operator<<(const char* str) {
    return *static_cast<std::ostream*>(this) << str;
  }
};
template <typename T>
bool basic_verbose_log<T>::verbose = false;
typedef basic_verbose_log<char> verbose_log;

#define vlog if(1) verbose_log()

#endif /* _VERBOSE_LOG_HPP_ */
