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
/*      binary_search_tree_iterator.h                                   */
/*      -----------------------------                                   */
/*      an implementation of an STL iterator for binary search trees.   */
/* original author: Nir Idisis                                          */
/************************************************************************/

#pragma once
#ifndef __BINARY_SEARCH_TREE_ITERATOR
#define __BINARY_SEARCH_TREE_ITERATOR

#include <iterator>
using std::iterator;
using std::bidirectional_iterator_tag;

template<typename T>
class binary_search_tree_iterator:
  public std::iterator<bidirectional_iterator_tag, typename T::value_type, typename T::size_type>
{
protected:
  const T *                tree; // tree pointed by iterator
  typename T::node_pointer node; // a pointer to a binary_node.

public:
  typedef typename T::pointer         pointer; /* pointer typedef (usually T*). */  
  typedef typename T::const_pointer   const_pointer; /* pointer typedef (usually const T*). */  
  typedef typename T::reference       reference; /* reference typedef (usually T&). */  
  typedef typename T::const_reference const_reference; /* reference typedef (usually const T&). */

  // Constructors:
  // -------------
  /* create a new empty instance of binary_search_tree_iterator. */
  explicit binary_search_tree_iterator() : tree(0) { };
  //      /* create a new instance of binary_search_tree_iterator. */
  //      binary_search_tree_iterator(const_reference that): node(that) { }
  /* create a new instance of binary_search_tree_iterator. */
  binary_search_tree_iterator(const T * tree_, typename T::node_pointer node_): tree(tree_), node(node_) { }

  // Increment / Decrement Operators:
  // --------------------------------
  /* go to the next node in the tree. */
  binary_search_tree_iterator operator++() {
    node = tree->successor(node);
    return *this;
 }
  /* go to the next node in the tree. */
  binary_search_tree_iterator operator++(int) {
    pointer temp = node;
    operator++();
    return binary_search_tree_iterator(tree, temp);
 }
  /* go to the previous node in the tree. */
  binary_search_tree_iterator operator--() {
    node = tree->predecessor(node);
    return (*this);
 }
  /* go to the previous node in the tree. */
  binary_search_tree_iterator operator--(int) {
    T * temp = node;
    operator--();
    return binary_search_tree_iterator(tree, temp);
 }
  /* assignment operator. */
  binary_search_tree_iterator & operator=(binary_search_tree_iterator rhs) {
    std::swap(tree, rhs.swap);
    std::swap(node, rhs.node);
    return *this;
 }
  /* compare two general_iterators. */
  bool operator==(const binary_search_tree_iterator & rhs) const { return (node == rhs.node); }
  /* compare two general_iterators (not equal). */
  bool operator!=(const binary_search_tree_iterator & rhs) const { return (node != rhs.node); }
  /* indirect access. */
  reference operator*() {
    typename T::node_reference n = tree->dereference(node);
    reference res = n.value;
    return res; 
  }
  pointer operator->() { return &(tree->dereference(node).value); }
};

#endif
