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


#ifndef __SKIP_LIST_MAP_HPP__
#define __SKIP_LIST_MAP_HPP__

#include <skip_list_set.hpp>
#include <utility>

template<typename Key, typename T, class Compare = std::less<Key>,
         int p_ = 4, typename Random = xor_random>
class skip_list_map : 
  public skip_list_set<std::pair<const Key,T>,first_comp<std::pair<const Key,T>,Compare>,p_,Random>
{
 public:
  typedef Key                    key_type;
  typedef T                      mapped_type;
  typedef std::pair<const Key,T> value_type;
 private:
  typedef skip_list_set<value_type,first_comp<value_type,Compare>,p_,Random> super;
 public:
  typedef Compare                         key_compare;
  typedef first_comp<value_type,Compare>  value_compare;
  typedef value_type&                     reference;
  typedef const value_type&               const_reference;
  typedef typename super::iterator        iterator;
  typedef typename super::const_iterator  const_iterator;
  typedef typename super::size_type       size_type;
  typedef typename super::difference_type difference_type;
  typedef typename super::pointer         pointer;
  typedef typename super::const_pointer   const_pointer;
  
  explicit skip_list_map(const Compare& comp = Compare(),
                         const Random& rand = Random()) : 
    super(value_compare(comp), rand) { }
  explicit skip_list_map(int max_height,
                         const Compare& comp = Compare(),
                         const Random& rand = Random()) :
    super(max_height, value_compare(comp), rand)
  { }

  template<typename InputIterator>
  skip_list_map(int max_height, InputIterator first, InputIterator last,
                const Compare& comp = Compare(),
                const Random& rand = Random()) :
    super(max_height, first, last, comp, rand) { }
  virtual ~skip_list_map() { }
  
  // Some methods have different interface between the map and set
  // data structure, in particular with respect to const-ness. Deal
  // with those here.
  std::pair<iterator, iterator> equal_range(const key_type& x) {
    return super::equal_range(x);
  }
  std::pair<const_iterator, const_iterator> equal_range(const key_type& x) const {
    return super::equal_range(x);
  }

  iterator find(const key_type& x) { 
    return super::find(x);
  }
  const_iterator find(const key_type& x) const { 
    return super::find(x);
  }

  iterator lower_bound(const key_type& x) { return super::lower_bound(x); }
  const_iterator lower_bound(const key_type& x) const { return super::lower_bound(x); }

  iterator upper_bound(const key_type& x) { return super::upper_bound(x); }
  const_iterator upper_bound(const key_type& x) const { return super::upper_bound(x); }

  // Extra subscript operator
  mapped_type& operator[](const key_type& x) {
    return (*((this->insert(std::make_pair(x,T()))).first)).second;
  }
  
  // Differences in the comparators
  key_compare key_comp() const { return super::value_comp.comp; }
  value_compare value_comp() const { return super::value_comp; }
};

#endif // __SKIP_LIST_MAP_HPP__
