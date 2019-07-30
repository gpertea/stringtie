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
#include <src2/getSuperReadInsertCountsFromReadPlacementFile_cmdline.hpp>
#include <src/bloom_filter.hpp>

using namespace std;

#define EXCLUDE_MATE 1

struct str_comp {
  bool operator()(const char *s1, const char *s2) {
    return strcmp(s1, s2) < 0;
  }
};

typedef const char*(*encode_fn)(const char* str);
typedef const char*(*decode_fn)(const char* str);

// Function pointer to encode (if so desired) the names.
const char* str_identity(const char* str) { return str; }
const char* dup(const char* str) { return strdup(str); }
const char* fib_decode(const char* str) {
  static charb buffer;
  sr_name::to_str(str, buffer);
  return buffer;
}

// Classes to filter names that need to be inserted in the map
class name_filter {
public:
  virtual bool insert(const char* name) = 0;
};
class name_bloom_filter : public name_filter {
  bloom_filter<const char*> filter;
public:
  name_bloom_filter(double fn, size_t n) : filter(fn, n) { }
  virtual bool insert(const char* const name) { return filter.add(name); }
};
class name_true : public name_filter {
public:
  virtual bool insert(const char* const name) { return true; }
};

int main (int argc, char **argv)
{
  cmdline_parse args(argc, argv);

  typedef map<const char *, int, str_comp> read_count_map;
  read_count_map    superReadToCounts;
  charb*            lines = new charb[2];
  charb*            cptr;
  int               line  = 1;
#if EXCLUDE_MATE
  char*             prefixHold;
  long long         readNumHold, readNum;
  char*             superReadHold;
#endif
  ExpBuffer<char *> flds;
  std::ifstream input(args.input_arg);
  if(!input.good())
    err::die(err::msg() << "Can't open input file '" << args.input_arg << "'" << err::no);
  std::ofstream output(args.output_arg);
  if(!output.good())
    err::die(err::msg() << "Can't open output file '" << args.output_arg << "'" << err::no);

  encode_fn encode;
  decode_fn decode;

  if(args.fib_flag) {
    encode = sr_name::from_str;
    decode = fib_decode;
  } else {
    encode = dup;
    decode = str_identity;
  }

  name_filter *filter;
  if(args.bloom_flag)
    filter = new name_bloom_filter(0.01, args.number_reads_arg);
  else
    filter = new name_true();

  cptr = &lines[(line = !line)];
  if(!getline (input, *cptr))
    err::die("Invalid input: empty file");
  getFldsFromLine (*cptr, flds);
  if(filter->insert(flds[1])) {
    const char *nstr = encode(flds[1]);
    ++superReadToCounts[nstr];
  }

#if EXCLUDE_MATE
  superReadHold = flds[1];
  prefixHold    = flds[0];      // Just use the first 2 chars
  readNumHold   = atoll (*cptr+2);
#endif
  cptr          = &lines[(line = !line)];
  while (getline (input, *cptr)) {
    getFldsFromLine (*cptr, flds);
#if EXCLUDE_MATE
    // The next excludes counting the second read of a mate pair from
    // the counts for the super-read. We want super-reads which have
    // reads coming from at least 2 different inserts. Our read counts
    // for the inserts will (probably) be off.
    readNum = atoll(*cptr+2);
    if ((strcmp(flds[1], superReadHold) == 0) &&
        (readNumHold+1 == readNum) &&
        (readNum % 2 == 1) &&
        (strncmp(prefixHold, *cptr, 2) == 0))
      continue;
#endif
    if(filter->insert(flds[1])) {
      const char *nstr = encode(flds[1]);
      ++superReadToCounts[nstr];
    }
#if EXCLUDE_MATE
    superReadHold = flds[1];
    prefixHold    = flds[0];
    readNumHold   = readNum;
#endif
    cptr          = &lines[(line = !line)];
  }
  for (read_count_map::const_iterator it = superReadToCounts.begin(); 
       it != superReadToCounts.end(); ++it)
    output << it->second << " " << decode(it->first) << "\n";

  return (0);
}
