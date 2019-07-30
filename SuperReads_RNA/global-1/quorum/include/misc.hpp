#ifndef __MISC_HPP__
#define __MISC_HPP__

//#define _GNU_SOURCE
#include <string.h>
#include <assert.h>
#include <algorithm>

// IO manipulator for substrings
class substr {

public:
  const char * const _str;
  size_t const       _len;
  substr(const char *str, size_t len) :
    _str(str), _len(len) {}
  substr(const char *str, const char *end) :
    _str(str), _len(end > str ? end - str : 0) {}
};
inline bool is_base(const char c) {
  switch(c) {
  case 'A': case 'C': case 'G': case 'T':
  case 'a': case 'c': case 'g': case 't':
    return true;

  default:
    return false;
  }
}
inline bool not_base(const char c) { return !is_base(c); }
inline std::ostream & operator<<(std::ostream &os, const substr &ss) {
  assert(std::count_if(ss._str, ss._str+ss._len, not_base) == 0);
  return os.write(ss._str, ss._len);
}

template<typename T>
int getFldsFromLine(char *line, T &res) {
  char *saveptr;
  res.clear();

  char *tok = strtok_r(line, " \t\n", &saveptr);
  while(tok) {
    res.push_back(tok);
    tok = strtok_r(0, " \t\n", &saveptr);
  }
  return res.size();
}

template<typename T>
int appendFldsFromLine(char *line, T &res) {
  char *saveptr;
  int numFlds = 0;

  char *tok = strtok_r(line, " \t\n", &saveptr);
  while(tok) {
    res.push_back(tok);
    ++numFlds;
    tok = strtok_r(0, " \t\n", &saveptr);
  }
  return numFlds;
}

#include <signal.h>

#define BREAKPOINT raise(SIGINT);

#endif
