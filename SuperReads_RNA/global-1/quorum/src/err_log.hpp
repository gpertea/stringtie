#ifndef __ERR_LOG_HPP__
#define __ERR_LOG_HPP__

#include <config.h>
#include <stdint.h>
#include <stddef.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <assert.h>
#include <misc.hpp>

template<typename T>
class err_log;

template<typename T>
std::ostream &
operator<<(std::ostream &os, const err_log<T> &e);

template<typename T>
class err_log {
public:
  struct sub_t {
    char from;
    char to;
  };

private:
  enum type_t { SUBSTITUTION, TRUNCATION };
  struct entry {
    type_t type;
    T      pos;
    union {
      sub_t sub;
    };

    entry() {};
    entry(type_t t, T p) : type(t), pos(p) {}
  };

  struct greater_than_pos : std::unary_function<entry, bool> {
    T pos;
    bool operator()(entry& e) { return e.pos >= *pos; } // Compare the raw value of pos, not using the operator>=
    greater_than_pos(T& pos_) : pos(pos_) { }
  };

  typedef std::vector<entry>  log_t;
  const unsigned int          _window;
  const unsigned int          _error;
  typename log_t::size_type   _lwin;
  const char                 *_trunc_string;
  log_t                       _log;
public:
  err_log(unsigned int window, unsigned int error, const char *trunc_string) : 
    _window(window), _error(error), _lwin(0), _trunc_string(trunc_string) { 
  }

  bool substitution(T pos, char from, char to) {
    assert(_log.size() == 0 ? true : pos > _log.back().pos);
    struct entry e(SUBSTITUTION, pos);
    e.sub.from = from;
    e.sub.to   = to;
    _log.push_back(e);
    return check_nb_error();
  }

  bool truncation(T pos) {
    assert(_log.size() == 0 ? true : pos >= _log.back().pos);
    _log.push_back(entry(TRUNCATION, pos));
    return check_nb_error();
  }

  // Remove all event log with position >= pos
  bool force_truncate(T pos) {
    greater_than_pos pred(pos);
    auto nend = std::remove_if(_log.begin(), _log.end(), pred);
    // The second argument is useless: we only truncate. But some
    // version of gcc complain if not present (e.g. 4.4).
    _log.resize(nend - _log.begin(), entry(SUBSTITUTION, 0));
    _lwin = 0; // Needs to be recomputed
    return check_nb_error();
  }

  // Check that the number of errors in the window is less than the
  // maximum allowed
  bool check_nb_error() {
    // Update the _lwin member
    if(_log.size() > 0 && _log.back().pos >_window)
      while(_log.back().pos > _log[_lwin].pos + _window) {
        _lwin++;
      }

    return _log.size() - _lwin - 1 >= _error;
  }

  int remove_last_window() {
    if(_log.size() == 0)
      return 0;
    int diff = _log.back().pos - _log[_lwin].pos;
    //    DBG << V(*_log.back().pos) << V(*_log[_lwin].pos) << V(diff);
    _log.erase(_log.begin() + _lwin, _log.end());
    _lwin = 0;
    check_nb_error();
    return diff;
  }

  friend std::ostream &operator<< <> (std::ostream &os, const err_log &l);
};

template<typename T>
std::ostream &operator<<(std::ostream &os, const err_log<T> &l) {
  typename err_log<T>::log_t::const_iterator it;
  bool                                       not_first = false;
  for(it = l._log.begin(); it != l._log.end(); it++) {
    if(not_first)
      os << " ";
    else
      not_first = true;
    switch(it->type) {
    case err_log<T>::SUBSTITUTION:
      os << it->pos << ":sub:" << it->sub.from << "-" << it->sub.to;
      break;

    case err_log<T>::TRUNCATION:
      os << it->pos << ":" << l._trunc_string;
      break;

    default:
      break;
    }
  }

  return os;
}

#endif
