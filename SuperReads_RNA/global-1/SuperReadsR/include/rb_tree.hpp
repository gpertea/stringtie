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


/************************************************************************/
/*      rb_tree.h                                                       */
/*      ---------                                                       */
/*      an implementation of a red-black balanced binary search         */
/*      tree. This class implements an STL-compatible container.        */
/*      The implementation is directly ported from CLR.                 */
/* original author: Nir Idisis                                          */
/************************************************************************/

#pragma once
#ifndef __RB_TREE_H
#define __RB_TREE_H

#include <functional>
#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>
#include <rb_node.hpp>
#include <binary_search_tree_iterator.hpp>

using std::less;
using std::vector;

namespace red_black {
  /* new/delete based node allocator */
  template<typename T>
  class new_delete_allocator {
  public:
    typedef node<T, void *>    value_type;
    typedef void *             pointer;
    typedef value_type &       reference;
    typedef const void *       const_pointer;
    typedef const value_type & const_reference;
    typedef size_t             size_type;
    typedef ptrdiff_t          difference_type;

    new_delete_allocator() { }
    ~new_delete_allocator() { }

    pointer alloc_construct(const value_type &init) { return (void*)(new value_type(init)); }
    void destroy_dealloc(pointer p) { delete (value_type *)p; }
    reference dereference(pointer p) const { return *(value_type*)p; }
    const_reference dereference(const_pointer p) const { return *(value_type*)p; }
    pointer nil() const { return 0; }
  };

  /* container based allocator with no destroy */
  template<typename T, typename C>
  class container_allocator {
  public:
    typedef node<T, size_t>    value_type;
    typedef size_t             pointer;
    typedef value_type &       reference;
    typedef const size_t       const_pointer;
    typedef const value_type & const_reference;
    typedef size_t             size_type;
    typedef ssize_t            difference_type;

    container_allocator() { }
    container_allocator(size_type size) : ary(size) { }
    ~container_allocator() { }

    pointer alloc_construct(const value_type &init) { ary.push_back(init); return ary.size() - 1; }
    void destroy_dealloc(pointer& p) { }
    reference dereference(pointer& p) const { return const_cast<reference>(ary[p]); }
    const_reference dereference(const_pointer& p) const { return ary[p]; }
    void clear() { ary.clear(); }
    pointer nil() const { return (pointer)-1; }

  private:
    C ary;
  };

  template<typename T, typename _Pr, typename A>
  class basic_tree {
  public:
    // type detinitions:
    typedef typename A::value_type                  node_type;
    typedef typename A::pointer                     node_pointer;
    typedef typename A::const_pointer               const_node_pointer;
    typedef typename A::reference                   node_reference;
    typedef typename A::const_reference             const_node_reference;
    // stl required type definitions:
    typedef T                                       value_type;
    typedef value_type *                            pointer;
    typedef value_type &                            reference;
    typedef value_type const *                      const_pointer;
    typedef value_type const &                      const_reference;
    typedef typename node_type::size_type           size_type;
    typedef binary_search_tree_iterator<basic_tree> iterator;
    typedef iterator const                          const_iterator;

#define DEREF(x) allocator.dereference(x)
#define NIL      (allocator.nil())
    node_pointer nil() const { return NIL; }
        
    /*    create a new empty instance of basic_tree. */
    basic_tree(): allocator(), m_root(NIL), m_size(0) { }
    virtual ~basic_tree() { destroy_tree(m_root); }

    /*    return the element which holds the biggest value in the tree.   */
    const_iterator max() const { return empty() ? end() : const_iterator(this, max(m_root)); }
    /*    return element which holds the smallest value in the tree.      */
    const_iterator min() const { return empty() ? end() : const_iterator(this, min(m_root)); }
    /**   find a node containing the given value in the tree.     */
    /*    if not found, return end().                                                     */
    const_iterator find(reference value) const { return const_iterator(this, find(m_root, value)); }

    /*    return the first element in the tree.   */
    const_iterator begin() const { return min(); }
    /*    return the element after the last element in the tree.  */
    const_iterator end() const { return const_iterator(this, NIL); }
    /*    reference to first and last elements. Result undefined if called on an empty tree */
    const_reference front() const { return DEREF(min(m_root)).value; }
    const_reference back() const { return DEREF(min(m_root)).value; }
    
    /*    return the root of this tree.           */
    const_iterator root() const { return const_iterator(m_root); }
    /*    return true if the tree is empty.       */
    bool empty() const { return (m_root == NIL); }
    /*    return the size of this tree.   */
    size_type size() const { return m_size; }

    /* insert a new element into the tree. */
    iterator insert(const_reference value) { return iterator(this, insert_node(value)); }
    /* insert a new element into the tree. */
    // iterator insert(iterator & where, const_reference value) {
    //   return insert(value);
    // }
    /* insert a new range of elements into the tree. */
    template<typename InputIter>
    void insert(InputIter _beg, InputIter _end) {
      for(InputIter idx = _beg; idx != _end; ++idx)
        insert(*idx);
    }
    /*    remove an element from the tree.        */
    void remove(value_type value) {
      node_pointer where = find(m_root, value);
      if(where != NIL) {
        node_pointer rem = remove(where);
        allocator.destroy_dealloc(rem);
        --m_size;
      }
    }

    // Check tree consistency
    bool is_consistent() const {
      if(m_root == NIL)
        return true;
      if(!is_black(m_root))
        return false;

      int max_depth = 0;
      int min_depth = std::numeric_limits<int>::max();
      if(!is_consistent(m_root, 0, max_depth, min_depth))
        return false;
      return max_depth == min_depth;
    }

    // methods on nodes
    bool is_leaf(const node_type &n) const { return n.left == NIL && n.right == NIL; }
    bool is_consistent(const_node_pointer n, int depth, int &max_depth, int &min_depth) const {
      //      std::cout << "is_consistent " << n << std::endl;
      assert(n != NIL);
      bool res = DEREF(n).parent == NIL ? is_black(n) : (!is_red(n) || !is_red(DEREF(n).parent));
      if(!res)
        return false;
      if(is_black(n))
        ++depth;
      if(DEREF(n).left == NIL || DEREF(n).right == NIL) {
        max_depth = std::max(max_depth, depth);
        min_depth = std::min(min_depth, depth);
      }
      if(DEREF(n).left != NIL)
        res = res && is_consistent(DEREF(n).left, depth, max_depth, min_depth);
      if(DEREF(n).right != NIL)
        res = res && is_consistent(DEREF(n).right, depth, max_depth, min_depth);
      return res;
    }
    bool is_a_left_son(const_node_pointer np) const { 
      const_node_reference n = DEREF(np);
      return n.parent != NIL && DEREF(n.parent).left == np;
    }
    bool is_a_right_son(const_node_pointer np) const {
      const_node_reference n = DEREF(np);
      return n.parent != NIL && DEREF(n.parent).right == np;
    }
    bool is_black(const_node_pointer np) const {
      return np == NIL || DEREF(np).color == BLACK;
    }
    bool is_red(const_node_pointer np) const {
      return !is_black(np);
    }
    node_pointer max(node_pointer n) const { 
      node_pointer right = DEREF(n).right;
      return right == NIL ? n : max(right);
    }
    node_pointer min(node_pointer n) const {
      node_pointer left = DEREF(n).left;
      return left == NIL ? n : min(left);
    }
    node_pointer successor(node_pointer np) const {
      if(np == NIL)
        return NIL;
      const_node_reference n = DEREF(np);
      if(n.right != NIL)
        return min(n.right);
      if(n.parent == NIL)
        return NIL;
      if(DEREF(n.parent).left == np) // np is a left child
        return n.parent;
      do { np = DEREF(np).parent; }
      while(DEREF(np).parent != NIL && np == DEREF(DEREF(np).parent).right);
      return DEREF(np).parent;
    }
    void destroy_tree(node_pointer np) {
      if(np == NIL)
        return;
      node_reference n = DEREF(np);
      if(n.left != NIL)
        destroy_tree(n.left);
      if(n.right != NIL)
        destroy_tree(n.right);
      allocator.destroy_dealloc(np);
    }
    node_reference dereference(node_pointer& np) const { return DEREF(np); }
    const_node_reference dereference(const_node_pointer& np) const { return DEREF(np); }

  protected:
    A            allocator;
    node_pointer m_root;          // a pointer to the root of this tree.
    size_type    m_size;          // the number of nodes in this tree.

    node_pointer find(node_pointer where, reference value) const {
      if(where == NIL)
        return NIL;
      if(DEREF(where).value == value)
        return where;      // value found!
      if(_Pr()(value, DEREF(where).value))
        return find(DEREF(where).left, value);         // continue searching in the left subtree:
      else
        return find(DEREF(where).right, value);        // continue searching in the right subtree:
    }

    /*    left rotate around <node>.      */
    void left_rotate(node_pointer node) {
      // check if node has a right son (if not do nothing):
      if(DEREF(node).right != NIL) {
        node_pointer right_son = DEREF(node).right;
        // node's new right-son will be it's old right-son's left-son.
        if((DEREF(node).right = DEREF(DEREF(node).right).left) != NIL)
          // update the new right-son of this node (if exists):
          DEREF(DEREF(node).right).parent = node;
        // update this node's parent:
        if(DEREF(node).parent == NIL)
          // node was the root of this tree.
          // the new root is node's right son:
          m_root = right_son;
        else if(is_a_left_son(node))
          DEREF(DEREF(node).parent).left = right_son;
        else // node is a right son.
          DEREF(DEREF(node).parent).right = right_son;       

        // update right_son's parent:
        DEREF(right_son).parent = DEREF(node).parent;
        // the left son of the old right son is this node:
        DEREF(node).parent = right_son;
        DEREF(right_son).left = node;
      }
    }
    /*    right rotate around <node>.     */
    void right_rotate(node_pointer node) {
      // check if node has a left son (if not do nothing):
      if(DEREF(node).left != NIL) {
        node_pointer left_son = DEREF(node).left;
        // node's new left-son will be it's old left-son's right-son.
        if((DEREF(node).left = DEREF(DEREF(node).left).right) != NIL)
          // update the new left-son of this node (if exists):
          DEREF(DEREF(node).left).parent = node;
        // update this node's parent:
        if(DEREF(node).parent == NIL)
          // node was the root of this tree.
          // the new root is node's right son:
          m_root = left_son;
        else if(is_a_left_son(node))
          DEREF(DEREF(node).parent).left = left_son;
        else // node is a right son.
          DEREF(DEREF(node).parent).right = left_son;        

        // update right_son's parent:
        DEREF(left_son).parent = DEREF(node).parent;
        // the right son of the old left son is this node:
        DEREF(node).parent = left_son;
        DEREF(left_son).right = node;
      }
    }
    /**   insert an existing node to the tree. Return true if
     *   rough_insert is enough (no fixing of the tree is
     *   necessary). This happens if the tree was empty of the value
     *   already exists in the tree.
     */
    bool rough_insert(node_pointer &nnode) {
      // check if tree is empty:
      if(empty()) {         // nnode is the new root of this tree:
        DEREF(nnode).parent = NIL;
        DEREF(nnode).color  = BLACK;
        m_root = nnode;
        m_size = 1;
        return true; // don't need to fixup
      }
      return rough_insert(m_root, nnode);
    }
    
    /** Insert the node where it belongs in the sub-tree rooted at
        where. If a node with the same value already exists, nnode is
        deallocated and nnode is changed to the existing node.
     */
    bool rough_insert(node_pointer where, node_pointer &nnode) {
      if(DEREF(nnode).value == DEREF(where).value) { // no insertion needed
        allocator.destroy_dealloc(nnode);
        nnode = where;
        return true; // don't need to fixup
      }
      // check to which side of where should be inserted:
      if(_Pr()(DEREF(nnode).value, DEREF(where).value)) {
        if(DEREF(where).left != NIL)
          return rough_insert(DEREF(where).left, nnode);
        // node is the new left subtree of where:
        DEREF(where).left   = nnode;
        DEREF(nnode).parent = where;
      } else {
        if(DEREF(where).right != NIL)
          return rough_insert(DEREF(where).right, nnode);
        // node is the new right subtree of whe
        DEREF(where).right  = nnode;
        DEREF(nnode).parent = where;
      }
      ++m_size; // Successfully inserted nnod
      return false;
    }

    /**   insert a given node into this tree.     */
    /*    return the given node.                          */
    node_pointer insert_node(const_reference value) {
      node_pointer node = allocator.alloc_construct(node_type(value, NIL, RED, NIL, NIL));
      // first, insert the node to the tree discarding red-black rules:
      if(rough_insert(node))
        return node;

      // now, turn the tree into a valid red-black tree:
      node_pointer x = node;
      DEREF(x).color = RED;
      while(x != m_root && DEREF(DEREF(x).parent).color == RED)  {
        if(is_a_left_son(DEREF(x).parent)) {
          node_pointer y = DEREF(DEREF(DEREF(x).parent).parent).right;
          if(y != NIL && DEREF(y).color == RED) {
            // Case 1 in the book:
            DEREF(DEREF(x).parent).color = BLACK;
            DEREF(y).color = BLACK;
            DEREF(DEREF(DEREF(x).parent).parent).color = RED;
            x = DEREF(DEREF(x).parent).parent;
          } else {  // y's color is black (the color of any leaf is always black).
            if(is_a_right_son(x)) {
              // Case 2 in the book:
              x = DEREF(x).parent;
              left_rotate(x);
            }
            // Case 3 in the book:
            DEREF(DEREF(x).parent).color               = BLACK;
            DEREF(DEREF(DEREF(x).parent).parent).color = RED;
            right_rotate(DEREF(DEREF(x).parent).parent);
          }
        } else if(is_a_right_son(DEREF(x).parent)) {
          node_pointer y = DEREF(DEREF(DEREF(x).parent).parent).left;
          if(y != NIL && DEREF(y).color == RED) {
            // Case 1 in the book:
            DEREF(DEREF(x).parent).color = BLACK;
            DEREF(y).color = BLACK;
            DEREF(DEREF(DEREF(x).parent).parent).color = RED;
            x = DEREF(DEREF(x).parent).parent;
          } else {  // y's color is black (the color of any leaf is always black).
            if(is_a_left_son(x)) {
              // Case 2 in the book:
              x = DEREF(x).parent;
              right_rotate(x);
            }
            // Case 3 in the book:
            DEREF(DEREF(x).parent).color               = BLACK;
            DEREF(DEREF(DEREF(x).parent).parent).color = RED;
            left_rotate(DEREF(DEREF(x).parent).parent);
          }
        }
      }
      // change the color of the root to black:
      DEREF(m_root).color = BLACK;
      return node;
    }

    /**   remove a given node from the subtree.                                   */
    /*    return a pointer to the node replacing the removed one. */
    node_pointer remove(node_pointer z) {
      node_pointer y = z;
      // Using notations from CLR RB-Delete
      if(DEREF(z).left != NIL && DEREF(z).right != NIL)
        y = successor(y);
      node_pointer x = DEREF(y).left != NIL ? DEREF(y).left : DEREF(y).right;
      node_pointer parent_x = DEREF(y).parent;
      if(x != NIL)
        DEREF(x).parent = parent_x;
      if(DEREF(y).parent == NIL) {
        m_root = x;
      } else {
        if(y == DEREF(DEREF(y).parent).left) {
          DEREF(DEREF(y).parent).left  = x;
        } else {
          DEREF(DEREF(y).parent).right = x;
        }
      }
      if(y != z)
        DEREF(z).value = DEREF(y).value;

      if(is_black(y))
        delete_fixup(x, parent_x);
    
      return y;
    }
  
    void delete_fixup(node_pointer x, node_pointer parent_x) {
      while(x != m_root && is_black(x)) {
        if(x == DEREF(parent_x).left) {
          node_pointer w = DEREF(parent_x).right;
          if(is_red(w)) { // Case 1
            DEREF(w).color         = BLACK;
            DEREF(DEREF(w).parent).color = RED;
            left_rotate(parent_x);
            if(x != NIL)
              parent_x = DEREF(x).parent;
            w = DEREF(parent_x).right;
          }
          if(is_black(DEREF(w).left) && is_black(DEREF(w).right)) { // Case 2
            DEREF(w).color = RED;
            x = parent_x;
            parent_x = DEREF(x).parent;
          } else {
            if(is_black(DEREF(w).right)) { // Case 3
              DEREF(DEREF(w).left).color = BLACK;
              DEREF(w).color = RED;
              right_rotate(w);
              w = DEREF(parent_x).right;
            }
            // Case 4
            DEREF(w).color         = DEREF(parent_x).color;
            DEREF(parent_x).color = BLACK;
            DEREF(DEREF(w).right).color  = BLACK;
            left_rotate(parent_x);
            x = m_root;
          }
        } else {
          node_pointer w = DEREF(parent_x).left;
          if(is_red(w)) { // Case 1
            DEREF(w).color         = BLACK;
            DEREF(DEREF(w).parent).color = RED;
            right_rotate(parent_x);
            if(x != NIL)
              parent_x = DEREF(x).parent;
            w = DEREF(parent_x).left;
          }
          if(is_black(DEREF(w).right) && is_black(DEREF(w).left)) { // Case 2
            DEREF(w).color = RED;
            x = parent_x;
            parent_x = DEREF(x).parent;
          } else {
            if(is_black(DEREF(w).left)) { // Case 3
              DEREF(DEREF(w).right).color = BLACK;
              DEREF(w).color = RED;
              left_rotate(w);
              w = DEREF(parent_x).left;
            }
            // Case 4
            DEREF(w).color         = DEREF(parent_x).color;
            DEREF(parent_x).color = BLACK;
            DEREF(DEREF(w).left).color  = BLACK;
            right_rotate(parent_x);
            x = m_root;
          }
        }
      }
      if(x != NIL)
        DEREF(x).color = BLACK;
    }
  };
}


#endif
