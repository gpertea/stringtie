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


#include <stdlib.h>

#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <memory>

#include <charb.hpp>
#include <jellyfish/err.hpp>
#include <src/sorted_merge_cmdline.hpp>

namespace err = jellyfish::err;

/** Like 'sort -m', but more limited and much faster. Supports sorting
 * only according to 1 key. Support for numerical and string sorting.
 */

template<typename T>
struct heap_elt {
  typedef T key_type;
  key_type       key;
  std::ifstream  is;
  const char *   file_path;
  charb          line;
  charb          tokens;

  heap_elt(const char *path) : is(path), file_path(path) { }
};

template<typename T>
bool greater(heap_elt<T> *e1, heap_elt<T> *e2) { return e1->key > e2->key; }

template<>
bool greater<const char *>(heap_elt<const char *> *e1, heap_elt<const char *> *e2) { 
  return strcmp(e1->key, e2->key) > 0;
}


template<typename T>
T conversion(const char *str);
template<>
long conversion<long>(const char *str) { return atol(str); }
template<>
const char *conversion<const char *>(const char *str) { return str; }


template<typename T>
void check_key(const T &key) { }

template<>
void check_key<const char *>(const char * const &key) {
  assert(key != 0);
  assert(strlen(key) > 0);
}

// Parse a line and set key of the element. Return false if no more lines (eof).
template<typename T>
bool parse_line(heap_elt<T> *elt, uint32_t col) {
  if(!getline(elt->is, elt->line))
    return false;
  strcpy(elt->tokens, elt->line);
  char *saveptr = 0;
  const char *key = strtok_r(elt->tokens, " \t", &saveptr);
  for(uint32_t i = 1; key && i < col; ++i)
    key = strtok_r(0, " \t", &saveptr);
  if(!key)
    err::die(err::msg() << "Invalid file '" << elt->file_path << "', invalid input line: " << elt->line);
  assert(key != 0);
  assert(strlen(key) > 0);
  elt->key = conversion<T>(key);
  check_key(elt->key);
  return true;
}

template<typename T, typename It>
void merge_streams(std::ostream &os, uint32_t col, It begin, It end) {
  std::make_heap(begin, end, greater<T>);

  It last = end - 1;
  while(begin != end) {
    std::pop_heap(begin, end, greater<T>);
    os << (*last)->line << "\n";
    if(parse_line<T>(*last, col)) {
      std::push_heap(begin, end, greater<T>);
    } else {
      end = last;
      --last;
    }
  }
}

template<typename T, typename It>
void merge_files(std::ofstream &os, uint32_t col, It begin, It end) {
  typedef std::vector<heap_elt<T>*> elt_vec;
  typedef typename elt_vec::iterator elt_iterator;
  elt_vec elts;

  for(It it = begin; it != end; ++it) {
    heap_elt<T>* new_elt = new heap_elt<T>(*it);
    if(!new_elt->is.good())
      err::die(err::msg() << "Failed to open input file '" << *it << "': " << err::no);
    if(parse_line(new_elt, col))
      elts.push_back(new_elt);
    else
      delete new_elt;
  }
  merge_streams<T, elt_iterator>(os, col, elts.begin(), elts.end());

   for(elt_iterator it = elts.begin(); it != elts.end(); ++it)
     delete *it;
}

int main(int argc, char *argv[]) {
  cmdline_parse args(argc, argv);
  
  std::ofstream out(args.output_arg);

  if(!out.good())
    err::die(err::msg() << "Failed to open output file '" << args.output_arg << "': " << err::no);
  if(args.numerical_flag)
    merge_files<long, cmdline_parse::input_arg_it>
      (out, args.key_arg, args.input_arg.begin(), args.input_arg.end());
  else
    merge_files<const char *, cmdline_parse::input_arg_const_it>
      (out, args.key_arg, args.input_arg.begin(), args.input_arg.end());

  return 0;
}
