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


#ifndef __SKIP_LIST_SET_HPP
#define __SKIP_LIST_SET_HPP

/* An ordered set class based on the skip list data structure, instead
   of the more usual balanced binary tree. The performance is similar
   to a balanced tree, i.e. log(n) operations for insertion and
   search, but in average instead of worst case.  

   The interface is almost a drop-in replacement for std::set, with
   the exception that the iterators forward iterators instead of
   bidirectional iterators. And consequently there is no reverse
   operators (no rbegin() and rend()).

   As it stands, there is no allocator parameter.
 */

#include <iterator>
#include <algorithm>
#include <cstring>
#include <stdint.h>
#include <skip_list_common.hpp>

// Set based on a skip list
template <typename Key, typename Compare = std::less<Key>, int p_ = 4, typename Random = xor_random>
class skip_list_set {
  struct node {
    Key   k;
    int   height;
    node* tower[];
  };
  node**                    heads_;
  int                       max_height_;
  int                       cur_height_; // 0 < cur_height <= max_height_
  size_t                    size_;
  Compare                   comp_;
  random_height<Random, p_> rh_;

public:
  typedef Key                                   key_type;
  typedef Key                                   value_type;
  typedef Compare                               key_compare;
  typedef Compare                               value_compare;
  typedef Key&                                  reference;
  typedef const Key&                            const_reference;
  typedef size_t                                size_type;
  typedef ssize_t                               difference_type;
  typedef Key*                                  pointer;
  typedef const Key*                            const_pointer;
  static const int p = p_;

public:
  class node_iterator {
  protected:
    node* item;
    node_iterator(node* item_) : item(item_) { }
    node_iterator(const node_iterator& rhs) : item(rhs.item) { }
    void next() { item = item->tower[0]; }
  public:
    bool operator==(const node_iterator& rhs) const { return item == rhs.item; }
    bool operator!=(const node_iterator& rhs) const { return item != rhs.item; }
  };

  class iterator : 
    public std::iterator<std::forward_iterator_tag, key_type>,
    public node_iterator {
    friend class skip_list_set;
    iterator(node* item_) : node_iterator(item_) { }
  public:
    iterator() : node_iterator(0) { }
    iterator(const node_iterator& rhs) : node_iterator(rhs) { }

    iterator& operator=(iterator rhs) {
      std::swap(node_iterator::item, rhs.item);
      return *this;
    }    
    reference operator*() { return node_iterator::item->k; }
    pointer operator->() { return &node_iterator::item->k; }
    iterator& operator++() {
      node_iterator::next();
      return *this;
    }
    iterator operator++(int) {
      iterator c(*this);
      node_iterator::next();
      return c;
    }
  };
  class const_iterator :
    public std::iterator<std::forward_iterator_tag, key_type>,
    public node_iterator {
    friend class skip_list_set;
    const_iterator(node* item_) : node_iterator(item_) { }
  public:
    const_iterator() : node_iterator(0) { }
    const_iterator(const const_iterator& rhs) : node_iterator(rhs) { }
    const_iterator(const iterator& rhs) : node_iterator(rhs) { }

    const_iterator& operator=(node_iterator rhs) {
      swap(node_iterator::item, rhs.item);
      return *this;
    }    
    const_reference operator*() { return node_iterator::item->k; }
    const_pointer operator->() { return &node_iterator::item->k; }
    const_iterator& operator++() {
      node_iterator::next();
      return *this;
    }
    const_iterator operator++(int) {
      const_iterator c(*this);
      node_iterator::next();
      return c;
    }
  };

  
  explicit skip_list_set(const Compare& comp = Compare(), 
                         const Random& rand = Random()) :
    heads_(new node*[10]), max_height_(10), cur_height_(1), size_(0),
    comp_(comp), rh_(rand)
  {
    memset(heads_, '\0', sizeof(node*) * max_height_);    
  }
                         
  explicit skip_list_set(int max_height,
                         const Compare& comp = Compare(), 
                         const Random& rand = Random()) :
    heads_(new node*[max_height]),
    max_height_(max_height), cur_height_(1), size_(0),
    comp_(comp), rh_(rand)
  {
    memset(heads_, '\0', sizeof(node*) * max_height_);
  }
  template<class InputIterator>
  skip_list_set(int max_height, InputIterator first, InputIterator last,
                const Compare& comp = Compare(), 
                const Random& rand = Random()) :
    heads_(new node*[max_height]),
    max_height_(max_height), cur_height_(1), size_(0),
    comp_(comp), rh_(rand)
  {
    memset(heads_, '\0', sizeof(node*) * max_height_);    
    insert(first, last);
  }    
  skip_list_set(const skip_list_set& rhs) :
    heads_(new node*[rhs.max_height_]),
    max_height_(rhs.max_height_), cur_height_(1), size_(0)
  {
    memset(heads_, '\0', sizeof(node*) * max_height_);    
    insert(rhs.begin(), rhs.end());
  }
  virtual ~skip_list_set() {
    clear();
    delete [] heads_;
  }

  skip_list_set& operator=(skip_list_set rhs) {
    swap(rhs);
    return *this;
  }

  void swap(skip_list_set& rhs) {
    std::swap(heads_, rhs.heads_);
    std::swap(max_height_, rhs.max_height_);
    std::swap(cur_height_, rhs.cur_height_);
    std::swap(size_, rhs.size_);
    std::swap(comp_, rhs.comp_);
    std::swap(rh_, rhs.rh_);
  }

  size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  iterator begin() { return iterator(heads_[0]); }
  const_iterator begin() const { return const_iterator(heads_[0]); }
  iterator end() { return iterator(); }
  const_iterator end() const { return const_iterator(); }

  void clear() {
    node* cnode = heads_[0];
    while(cnode) {
      node* nnode = cnode->tower[0];
      delete cnode;
      cnode = nnode;
    }
    memset(heads_, '\0', sizeof(node*) * max_height_);
    size_ = 0;
  }

  template<typename T>
  size_type count(const T& x) const {
    return find_node(x) ? 1 : 0;
  }

  template<typename T>
  std::pair<iterator, iterator> equal_range(const T& x) const {
    node **path[max_height_];
    node *n = find_node_path(x, path);
    if(n)
      return std::make_pair(iterator(n), ++iterator(n));
    return std::make_pair(iterator(*path[0]), iterator(*path[0]));
  }

  size_type erase(const key_type& x) {
    node **path[max_height_];
    node *n = find_node_path(x, path);
    if(!n)
      return 0;
    for(int i = 0; i < n->height; ++i)
      *path[i] = n->tower[i];
    delete n;
    --size_;
    return 1;
  }

  template<typename T> // T typically key_type
  iterator find(const T& x) const {
    return iterator(find_node(x));
  }

  key_compare key_comp() const { return comp_; }
  value_compare value_comp() const { return comp_; }

  template<typename T> // T typically key_type
  iterator lower_bound(const T& x) const {
    node** path[max_height_];
    find_node_path(x, path);
    return iterator(*path[0]);
  }

  template<typename T>
  iterator upper_bound(const T& x) const {
    node** path[max_height_];
    node*  n = find_node_path(x, path);
    if(n)
      return ++iterator(n);
    return iterator(*path[0]);
  }

  size_type max_size() const {
    size_type res = 1;
    size_type pp  = p;
    int       n   = max_height_;

    while(n) {
      if(n & 0x1)
        res *= pp;
      pp *= pp;
      n >>= 1;
    }
    return res;
  }

  std::pair<iterator, bool> insert(const value_type& x) {
    node** path[max_height_];
    node*  n = find_node_path(x, path);
    if(n)
      return std::make_pair(iterator(n), false);
    n = new_node(x);
    const int height = std::min(n->height, cur_height_);
    int i;
    for(i = 0; i < height; ++i) {
      n->tower[i] = *path[i];
      *path[i]    = n;
    }
    if(n->height > cur_height_) {
      for(i = cur_height_; i < n->height; ++i) {
        n->tower[i] = 0; 
        heads_[i]   = n;
      }
      cur_height_ = n->height;
    }
    ++size_;
    return std::make_pair(iterator(n), true);
  }

  // Unlike std::set, this provides no advantage at this point.
  iterator insert(iterator position, const value_type& x) {
    return insert(x).first;
  }

  template<class InputIterator>
  void insert(InputIterator first, InputIterator last) {
    for( ; first != last; ++first)
      insert(*first);
  }
  
protected:
  // Find the node path. I.e. the addresses of the pointers that leads
  // to the largest element less than x. If a node with a value
  // equivalent to x is found, a pointer to the node is returned and
  // path has an undefined value. Otherwise, 0 is returned and the
  // path array is initialized.
  //
  // T is whatever type that can be compared with the elements in the set
  template<typename T>
  node* find_node_path(const T& x, node*** path) const {
    int i = cur_height_ - 1;
    for( ; i >= 0 && !(heads_[i] && comp_(heads_[i]->k, x)); --i)
      path[i] = &heads_[i];
    if(i >= 0) {
      for(node* cnode = heads_[i]; i >= 0; --i) {
        node* nnode = cnode->tower[i];
        while(nnode && comp_(nnode->k, x)) {
          cnode = nnode;
          nnode = nnode->tower[i];
        }
        path[i] = &(cnode->tower[i]);
      }
    }
    // Check if we found a node equal to x. If so, return it.
    node* lnode = *path[0];
    return lnode && !comp_(x, lnode->k) ? lnode : 0;
  }

  // Find a node. Equivalent to find_node_path, but the path is not
  // recorded. Maybe a tad faster.
  template<typename T>
  node* find_node(const T& x) const {
    int i = cur_height_ - 1;
    for( ; i >= 0 && !(heads_[i] && comp_(heads_[i]->k, x)); --i) ;
    if(i < 0) {
      if(heads_[0] && !comp_(x, heads_[0]->k))
        return heads_[0];
      return 0;
    }

    for(node* cnode = heads_[i]; i >= 0; --i) {
      node* nnode = cnode->tower[i];
      while(nnode && comp_(nnode->k, x)) {
        cnode = nnode;
        nnode = nnode->tower[i];
      }
      if(nnode && !comp_(x, nnode->k))
        return nnode;
    }
    return 0;
  }

  // Allocate a new node. Does raw memory allocation of a node with
  // enough space for the tower. Then in place copy construction of
  // the key from x.
  node* new_node(const value_type& x) {
    //    int height  = std::min(max_height_, std::min(cur_height_ + 1, rh_()));
    int height = std::min(max_height_, rh_());
    node* res   = (node*)operator new(sizeof(node) + height * sizeof(node*));
    res->height = height;
    new ((void*)&res->k) value_type(x);
    return res;
  }

  // // Debugging routines.
  // void print_lists(std::ostream &os = std::cout) const {
  //   for(int i = max_height_ - 1; i >= 0; --i) {
  //     node* cnode = heads_[i];
  //     os << i << ": " << (void*)&heads_[i];
  //     do {
  //       os << " " << (void*)cnode;
  //       if(cnode) {
  //         os << ":" << cnode->k;
  //         cnode = cnode->tower[i];
  //       }
  //     } while(cnode);
  //     os << "\n";
  //   }
  //   os << std::flush;
  // }

  // void check_path(const value_type& x, node*** path) const {
  //   for(int i = cur_height_ - 1; i >= 0; --i) {
  //     std::cout << "check " << i << ": " << (void*)path[i] << ":"
  //               << (path[i] ? (void*)*path[i] : (void*)0) << std::endl;
  //   }
  // }  
};

// Equality operators between iterators
// template <typename Key, typename Compare, int p_, typename Random>
// bool operator==(skip_list_set<Key, Compare, p_, Random::iterator,
//                 skip_list_set<Key, Compare, p_, Random::const_iterator) {
//   return 
// }

#endif /* __SKIP_LIST_SET_HPP */
