#ifndef __ERROR_CORRECT_READS_HPP__
#define __ERROR_CORRECT_READS_HPP__

#include <config.h>
#include <jellyfish/err.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/thread_exec.hpp>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <misc.hpp>
#include <src/kmer.hpp>
#include <src/err_log.hpp>

// A forward moving pointer.
template<typename T>
class forward_ptr {
  T *_p;
public:
  forward_ptr(T *p) : _p(p) { }
  forward_ptr operator++(int) {
    return forward_ptr(_p++);
  }
  forward_ptr &operator++() {
    ++_p;
    return *this;
  }
  forward_ptr operator--(int) {
    return forward_ptr(_p--);
  }
  forward_ptr &operator--() {
    --_p;
    return *this;
  }
  forward_ptr operator+(int x) const {
    return forward_ptr(_p + x);
  }
  forward_ptr operator-(int x) const {
    return forward_ptr(_p - x);
  }
  T &operator*() { return *_p; }
  bool operator<(T *e) const { return _p < e; }
  bool operator<(const forward_ptr<T> &e) const { return _p < e._p; }
  bool operator>=(T *e) const { return _p >= e; }
  bool operator>=(const forward_ptr<T> &e) const { return _p >= e._p; }
  ptrdiff_t operator-(T *s) const { return _p - s; }
  T *ptr() const { return _p; }
};

// A backward moving pointer. suffix++ behaves like suffix--. < behave
// likes >. Etc.
template<typename T>
class backward_ptr {
  T *_p;
public:
  backward_ptr(T *p) : _p(p) {}

  backward_ptr operator++(int) {
    return backward_ptr(_p--);
  }
  backward_ptr &operator++() {
    --_p;
    return *this;
  }
  backward_ptr operator--(int) {
    return backward_ptr(_p++);
  }
  backward_ptr &operator--() {
    ++_p;
    return *this;
  }
  backward_ptr operator+(int x) const {
    return backward_ptr(_p - x);
  }
  backward_ptr operator-(int x) const {
    return backward_ptr(_p + x);
  }
  T &operator*() const { return *_p; }
  bool operator<(T *e) const { return _p > e; }
  bool operator<(const backward_ptr<T> &e) const { return _p > e._p; }
  bool operator>=(T *e) const { return _p <= e; }
  bool operator>=(const backward_ptr<T> &e) const { return _p <= e._p; }
  ptrdiff_t operator-(T *s) const { return _p - s; }
  T *ptr() const { return _p; }
};

class forward_counter {
  int _c;
public:
  forward_counter() {}
  forward_counter(int c) : _c(c) {}
  forward_counter &operator++() {
    ++_c;
    return *this;
  }
  forward_counter &operator--() {
    --_c;
    return *this;
  }
  forward_counter operator+(const int x) const {
    return forward_counter(_c + x);
  }
  forward_counter operator-(const int x) const {
    return forward_counter(_c - x);
  }
  int operator-(const forward_counter &c) const {
    return _c - c._c;
  }
  bool operator>(const forward_counter &c) const {
    return _c > c._c;
  }
  bool operator>=(const forward_counter &c) const {
    return _c > c._c || _c == c._c;
  }
  int operator*() const { return _c; }
  friend std::ostream &operator<<(std::ostream &os, const forward_counter &c);
};

class backward_counter {
  int _c;
public:
  backward_counter() {}
  backward_counter(int c) : _c(c) {}
  backward_counter &operator++() {
    --_c;
    return *this;
  }
  backward_counter &operator--() {
    ++_c;
    return *this;
  }
  bool operator>(const backward_counter &c) const {
    return _c < c._c;
  }
  bool operator>=(const backward_counter &c) const {
    return _c < c._c || _c == c._c;
  }
  backward_counter operator+(const int x) const {
    return backward_counter(_c - x);
  }
  backward_counter operator-(const int x) const {
    return backward_counter(_c + x);
  }
  int operator-(const backward_counter &c) const {
    return c._c - _c;
  }
  int operator*() const { return _c; }
  friend std::ostream &operator<<(std::ostream &os, const backward_counter &c);
};

inline std::ostream &operator<<(std::ostream &os, const forward_counter &c) {
  return os << c._c;
}
inline std::ostream &operator<<(std::ostream &os, const backward_counter &c) {
  return os << c._c;
}


class forward_log : public err_log<forward_counter> {
public:
  forward_log(unsigned int window, unsigned int error) :
    err_log<forward_counter>(window, error, "3_trunc") { }
};

class backward_log : public err_log<backward_counter> {
public:
  backward_log(unsigned int window, unsigned int error) :
    err_log<backward_counter>(window, error, "5_trunc") { }
  
  bool truncation(backward_counter pos) {
    return err_log<backward_counter>::truncation(pos - 1);
  }
};

#endif
