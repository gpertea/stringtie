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
/*      rb_node.h                                                       */
/*      ---------                                                       */
/*      this class defines a reb-black binary node of a red-black tree. */
/* original author: Nir Idisis                                          */
/************************************************************************/

#pragma once
#ifndef __RB_NODE_H
#define __RB_NODE_H

namespace red_black {
  enum { RED, BLACK };

  template<typename T, typename P>
  struct node {
    typedef T             value_type;
    typedef unsigned long size_type;
    T    value;
    P    parent;
    P    left;
    P    right;
    char color;
    node(value_type _value, P _parent, bool _color, P _left, P _right) :
      value(_value), parent(_parent), left(_left), right(_right), color(_color) { }

    // ----------
    // operators:
    // ----------
    /*      assignment operator for values. */
    value_type & operator=(const value_type & _value) { value = _value; }
    /*      return true if this node is equal to the given value.   */
    bool operator==(const value_type & _value) const { return value == _value; }
    /*      return true if this node is different from the given value.     */
    bool operator!=(const value_type & _value) const { return value != _value; }
    /*      return true if this node is bigger than the given value.        */
    bool operator>(const value_type & _value) const { return value > _value; }
    /*      return true if this node is smaller than the given value.       */
    bool operator<(const value_type & _value) const { return value < _value; }
    /*      return true if this node is bigger or equal to the given value. */
    bool operator>=(const value_type & _value) const { return value >= _value; }
    /*      return true if this node is smaller or equal to the given value.        */
    bool operator<=(const value_type & _value) const { return value <= _value; }
    /*      return true if this node's value is bigger than the value of the given node.    */
    bool operator>(const node & that)  const { return value > that.value; }
    /*      return true if this node's value is smaller than the value of the given node.   */
    bool operator<(const node & that)  const { return value < that.value; }
    /*      return true if this node's value is bigger or equal to the value of the given node.     */
    bool operator>=(const node & that) const { return value >= that.value; }
    /*      return true if this node's value is smaller or equal to the value of the given node.    */
    bool operator<=(const node & that) const { return value <= that.value; }
    /*      return the node which holds the biggest value in the subtree of which this node is the root.    */
    /* get the next node in the search tree. */
  };
}
#endif
