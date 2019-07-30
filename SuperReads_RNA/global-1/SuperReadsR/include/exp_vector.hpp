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


#ifndef __EXP_VECTOR_HPP__
#define __EXP_VECTOR_HPP__

#include <memory>
#include <vector>

template <typename T, class Allocator = std::allocator<T> >
class exp_vector : public std::vector<T, Allocator> {
  typedef std::vector<T, Allocator> super;
  T def_val;
public:
  typedef typename super::reference              reference;
  typedef typename super::const_reference        const_reference;
  typedef typename super::iterator               iterator;
  typedef typename super::const_iterator         const_iterator;
  typedef typename super::size_type              size_type;
  typedef typename super::difference_type        difference_type;
  typedef typename super::value_type             value_type;
  typedef typename super::allocator_type         allocator_type;
  typedef typename super::pointer                pointer;
  typedef typename super::const_pointer          const_pointer;
  typedef typename super::reverse_iterator       reverse_iterator;
  typedef typename super::const_reverse_iterator const_reverse_iterator;

  explicit exp_vector(const Allocator& a = Allocator()) : super(a), def_val(T()) { }
  explicit exp_vector (size_type n, const T& value = T(), const Allocator& a = Allocator()) :
    super(n, value, a), def_val(value) { }
  exp_vector(size_type n, T&& value, const Allocator& a = Allocator()) :
    super(n, value, a), def_val(std::move(value)) { }
  explicit exp_vector(const T& value, const Allocator& a = Allocator()) :
    super(a), def_val(value) { }
  template <class InputIterator>
  exp_vector ( InputIterator first, InputIterator last, const Allocator& a = Allocator() ) :
    super(first, last, a), def_val(T()) { }
  exp_vector ( const std::vector<T,Allocator>& x ) : super(x) { }
  exp_vector(const exp_vector& rhs) : super(rhs), def_val(rhs.def_val) { }
  exp_vector(exp_vector&& rhs) : super(std::move(rhs)), def_val(std::move(rhs.dev_fal)) { }

  virtual ~exp_vector() { }

  exp_vector& operator=(const exp_vector& rhs) {
    super::operator=(rhs);
    def_val = rhs.def_val;
    return *this;
  }
  exp_vector& operator=(exp_vector&& rhs) {
    super::operator=(std::move(rhs));
    def_val = std::move(rhs.def_val);
    return *this;
  }

  reference default_value() { return def_val; }
  const_reverse_iterator default_value() const { return def_val; }
  void default_value(const T& df) { def_val = df; }

  reference operator[] ( size_type n ) {
    if(n >= super::size())
      super::resize(n + 1, def_val);
    return super::operator[](n);
  }
  const_reference operator[] ( size_type n ) const {
    if(n >= super::size())
      super::resize(n + 1, def_val);
    return super::operator[](n);
  }
};

#endif /* __EXP_VECTOR_HPP__ */
