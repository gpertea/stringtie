#ifndef __KMER_HPP__
#define __KMER_HPP__

#include <stdint.h>
#include <assert.h>
#include <iostream>
#include <string>

#include <jellyfish/mer_dna.hpp>

class kmer_t {
  jellyfish::mer_dna _fmer, _rmer;

public:
  bool shift_left(char c) {
    int x = jellyfish::mer_dna::code(c);
    if(x != -1) {
      shift_left(x);
      return true;
    }
    return false;
  }

  void shift_left(int x) {
    _fmer.shift_left(x);
    _rmer.shift_right(mer_dna::complement(x));
  }

  bool shift_right(char c) {
    int x = jellyfish::mer_dna::code(c);
    if(x != -1) {
      shift_right(x);
      return true;
    }
    return false;
  }

  void shift_right(int x) {
    _fmer.shift_right(x);
    _rmer.shift_left(mer_dna::complement(x));
  }

  const jellyfish::mer_dna& canonical() const { return _fmer < _rmer ? _fmer : _rmer; }
  const jellyfish::mer_dna& fmer() const { return _fmer; }
  const jellyfish::mer_dna& rmer() const { return _rmer; }

  void replace(int i, int x) {
    _fmer.base(i)                               = x;
    _rmer.base(jellyfish::mer_dna::k() - i - 1) = jellyfish::mer_dna::complement(x);
  }

  char base(int i) const { return jellyfish::mer_dna::rev_code(code(i)); }
  int code(int i) const { return _fmer.base(i).code(); }
  //  uint64_t rcode(int i) const { assert(i >= 0 && i < _k); return (_rmer >> (2*i)) & c3; }
  std::string str() const { return _fmer.to_str(); }
  std::string rstr() const { return _rmer.to_str(); }

  friend std::ostream &operator<<(std::ostream &os, const kmer_t &mer);
  friend class forward_mer;
  friend class backward_mer;
};

std::ostream &operator<<(std::ostream &os, const kmer_t &mer) {
  return os << mer.str();
}

class forward_mer;
class backward_mer;

class forward_mer {
  kmer_t& _m;
public:
  forward_mer(kmer_t& m) : _m(m) {}
  backward_mer rev_mer() const;
  bool shift(char c) { return _m.shift_left(c); }
  void shift(int x) { _m.shift_left(x); }
  bool rev_shift(char c) { return _m.shift_right(c); }
  void rev_shift(int x) { _m.shift_right(x); }
  const jellyfish::mer_dna& canonical() const { return _m.canonical(); }
  char base(int i) const { return _m.base(i); }
  int code(int i) const { return _m.code(i); }
  void replace(int i, int x) { _m.replace(i, x); }
  const jellyfish::mer_dna& rmer() const { return _m.rmer(); }
  const kmer_t& kmer() const { return _m; }
  friend std::ostream &operator<<(std::ostream &os, const forward_mer &mer);
};
inline std::ostream &operator<<(std::ostream &os, const forward_mer &mer) {
  return os << mer._m.str();
}

class backward_mer {
  kmer_t& _m;
public:
  backward_mer(kmer_t& m) : _m(m) {}
  forward_mer rev_mer() const;
  bool shift(char c) { return _m.shift_right(c); }
  void shift(int x) { _m.shift_right(x); }
  bool rev_shift(char c) { return _m.shift_left(c); }
  void rev_shift(int x) { _m.shift_left(x); }
  const jellyfish::mer_dna& canonical() const { return _m.canonical(); }
  char base(int i) const { return _m.base(jellyfish::mer_dna::k() - i - 1); }
  int code(int i) const { return _m.code(jellyfish::mer_dna::k() - i - 1); }
  void replace(int i, uint64_t c) { _m.replace(jellyfish::mer_dna::k() - i - 1, c); }
  const kmer_t& kmer() const { return _m; }
  friend std::ostream &operator<<(std::ostream &os, const backward_mer &mer);
};
inline std::ostream &operator<<(std::ostream &os, const backward_mer &mer) {
  return os << mer._m.str();
}

inline backward_mer forward_mer::rev_mer() const {
  return backward_mer(_m);
}
inline forward_mer backward_mer::rev_mer() const {
  return forward_mer(_m);
}

#endif
