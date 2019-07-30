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


#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <exp_buffer.hpp>
#include <charb.hpp>
#include <jellyfish/err.hpp>
#include <misc.hpp> // for getFldsFromLine
#include <src/sr_names.hpp>
#include <src2/getSuperReadInsertCountsFromReadPlacementFileTwoPasses_cmdline.hpp>
#include <src/bloom_counter2.hpp>

using namespace std;

#define EXCLUDE_MATE 1

struct str_comp {
  bool operator()(const char *s1, const char *s2) {
    return strcmp(s1, s2) < 0;
  }
};

typedef const char*(*coding_fn)(const char* str);

// Function pointer to encode (if so desired) the names.
const char* str_identity(const char* str) { return str; }
const char* str_dup(const char* str) { return strdup(str); }
const char* fib_decode(const char* str) {
  static charb buffer;
  sr_name::to_str(str, buffer);
  return buffer;
}

// Parse the input input stream is. Delegate to s the decision of
// where to store the read name
template<typename store>
size_t parse_store_input(std::istream &is, store &s) {
  const char*       superReadHold = "";
  const char*       prefixHold    = "";
  long long         readNumHold   = 0;
  
  ExpBuffer<char *> flds;
  charb*            lines = new charb[2];
  int               line  = 0;
  charb*            cptr  = &lines[line];

  size_t            inserted = 0;

  while(getline(is, *cptr)) {
    getFldsFromLine(*cptr, flds);
    long long readNum = atoll(*cptr+2);
    if(strcmp(flds[1], superReadHold) == 0 &&
       readNumHold+1 == readNum &&
       readNum % 2 == 1 &&
       strncmp(prefixHold, *cptr, 2) == 0)
      continue; // Exclude second read from a mate-pair
    s(flds[1]);
    ++inserted;
    superReadHold = flds[1];
    prefixHold    = flds[0];
    readNumHold   = readNum;
    cptr = &lines[(line = !line)];
  }

  return inserted;
}

// Store in a bloom filter
struct bloom_store {
  bloom_counter2<const char*> bc;
  size_t existing;
  bloom_store(size_t n) : bc(0.01, n), existing(0) { }
  void operator()(const char *s) { 
    existing += bc.insert(s) >= 1;
  }
};

// Store in a map provided that the entry has been seen more than once
struct map_store {
  bloom_counter2<const char*>&                 bc;
  typedef std::map<const char*, int, str_comp> map_type;
  typedef map_type::const_iterator             iterator;
  map_type                                     map;
  coding_fn                                    encode;
  size_t                                       inserted;
  size_t                                       distinct;
  map_store(bloom_counter2<const char*>& bc_, coding_fn encode_) :
    bc(bc_), map(), encode(encode_), inserted(0), distinct(0) { }
  void operator()(const char *s) {
    if(bc.check(s) > (unsigned int)1) {
      ++inserted;
      const char* to_insert = encode(s);
      auto insert_res = map.insert(std::make_pair(to_insert, 1));
      if(insert_res.second)
        ++distinct;
      else {
        free((void*)to_insert);
        ++insert_res.first->second;
      }
    }
  }
};

int main (int argc, char **argv)
{
  cmdline_parse args(argc, argv);

  std::ofstream output(args.output_arg);
  if(!output.good())
    die(err::msg() << "Can't open output file '" << args.output_arg << "':" << err::no);

  if(args.debug_flag)
    std::cerr << "First pass" << std::endl;
  // Parse input into bloom counter
  bloom_store bs(args.number_reads_arg);
  for(auto file = args.input_arg.begin(); file != args.input_arg.end(); ++file) {
    if(args.debug_flag)
      std::cerr << "Parsing " << *file << std::endl;
    std::ifstream input(*file);
    size_t inserted = parse_store_input(input, bs);
    if(args.debug_flag)
      std::cerr << "Inserted " << inserted << " existing " << bs.existing << std::endl;
  }
  
  if(args.debug_flag)
    std::cerr << "Second pass" << std::endl;
  // Parse input into map, if count > 1
  coding_fn encode = args.fib_flag ? (coding_fn)sr_name::from_str : str_dup;
  map_store ms(bs.bc, encode);
  for(auto file = args.input_arg.begin(); file != args.input_arg.end(); ++file) {
    if(args.debug_flag)
      std::cerr << "Parsing " << *file << std::endl;
    std::ifstream input(*file);
    size_t inserted = parse_store_input(input, ms);
    if(args.debug_flag)
      std::cerr << "Elements " << inserted << " inserted " << ms.inserted << " distinct " << ms.distinct << std::endl;
  }

  // Output result
  coding_fn decode = args.fib_flag ? fib_decode : str_identity;
  size_t output_nb = 0;
  for (map_store::iterator it = ms.map.begin(); it != ms.map.end(); ++it)
    if(it->second > 1) {
      ++output_nb;
      output << it->second << " " << decode(it->first) << "\n";
    }
  output.close();
  if(args.debug_flag)
    std::cerr << "Output " << output_nb << std::endl;

  return (0);
}
