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


#ifndef __MULTIPLEXED_PARSER_HPP__
#define __MULTIPLEXED_PARSER_HPP__

#include <pthread.h>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <jflib/pool.hpp>
// #include <err.hpp>

/** A multi-thread safe mostly lock-free data parser.
 * 
 * This is the base class to implement multi-threaded parsers. One
 * need to implement the virtual function `void parser_loop()` that
 * does the actual parsing of the input into logical blocks
 * (a.k.a. elements) of type `T`, the template parameter. The actual
 * parsing is done in an independent thread which is managed by the
 * multiplexed_parser class.
 *
 * A (consumer) thread then creates a [stream](@ref
 * multiplexed_parser::stream) of type `T` to read logical blocks one
 * by one. Internally, each stream communicates with the parsing
 * thread through a mostly lock-free FIFO.
 *
 * For efficiency, the logical blocks are grouped together. The size
 * of the groups is specified at construction.
 */
template <typename T>
class multiplexed_parser {
  struct group {
    typedef std::vector<T>                 elt_vector;
    typedef typename elt_vector::size_type size_type;
    
    elt_vector elements;
    size_type  nb_filled;
  };

  typedef jflib::pool<group>        pool;
public:
  typedef typename group::size_type size_type;

private:
  const size_type group_size_;
  pool            pool_;
  bool            parser_started_;
  pthread_t       parser_id_;
  const char*     error_;

public:
  /** Constructor for a parser with a given number of threads
   *
   * @param nb_threads Internally sets the size of the lock-free
   * FIFO. Must be larger than the number of consumer threads which
   * will create a [stream](@ref multiplexed_parser::stream).
   *  @param group_size Size of groups of logical units.
   */  explicit multiplexed_parser(int nb_threads = 16,  size_type group_size = 100) :
    group_size_(group_size), pool_(3 * nb_threads), 
    parser_started_(false), error_(0) 
  { }
  virtual ~multiplexed_parser();

  /** Is the multiplexed_parser ready for IO. I.e. it has not encountered an error and it is not at EOF. */
  bool good() const { return !pool_.is_closed_A_to_B() && error() == 0; }
  /** Has the multiplexed_parser reached EOF */
  bool eof() const { return pool_.is_closed_A_to_B(); }
  /** Has an error occurred. EOF is not an error */
  bool fail() const { return error() != 0; }
  /** Error message */
  const char* error() const { return jflib::a_load_ptr((const char*&)error_); }

  /** Start the parsing thread. */
  void start_parsing();

  /** Size of the groups. */
  size_type group_size() const { return group_size_; }

  typedef typename pool::elt elt;
  /** Stream of elements. The next element is obtained with the `++`
   *  operator. The stream itself can be used in the test of a `while`
   *  loop and will test as false when no more elements are available.
   *
   *  Like an iostream, this stream is not copyable. Like an iterator,
   *  the `*` and `->` operator give access to an object of type `T`.
   *
   * The typical usage is as follows:
   * ~~~{.cc}
   * // Create parser. Class file_parser inherits from multiplexed_parser<std::string>
   * file_parser parser("/path/to/file");
   * 
   * // Use parser with a stream: print every element.
   * file_parser::stream stream(parser);
   * for( ; stream; ++stream) {
   *   std::cout << *stream << "\n";
   * }
   * ~~~
   */
  class stream { 
    pool&     pool_;
    elt elt_;
    size_type i_;
  public:
    typedef T value_type;

    /** Construct a stream from a multiplexed_parser. 
     * @param rp The multiplexed parser
     */
    explicit stream(multiplexed_parser& rp) :
      pool_(rp.pool_), elt_(pool_.get_B()), i_(0) 
    { 
      while(!elt_.is_empty() && elt_->nb_filled == (size_type)0)
        elt_ = pool_.get_B();
    }
    // Probably useless
    stream() = default;
    /** Copy constructor is deleted */
    stream(const stream& rhs) = delete;
    /** Copy constructor is deleted */
    stream& operator=(const stream& rhs) = delete;

    /** Get access to the element. */
    T& operator*() { return elt_->elements[i_]; }
    /** Get access to the element. */
    T* operator->() { return &elt_->elements[i_]; }

    /** Get next element. */
    stream& operator++() {
      if(++i_ < elt_->nb_filled)
        return *this;
      i_ = 0;
      do {
        elt_ = pool_.get_B();
      } while(!elt_.is_empty() && elt_->nb_filled == 0);
      return *this;
    }
    /** Evaluate to `(void*)0` if no more element available, to a non-zero pointer otherwise. */
    operator void*() const { return elt_.is_empty() ? (void*)0 : (void*)&elt_; }
  };
  
protected:
  /** Perform the actual parsing.
   *
   * Ideally, this methods does the minimum amount of work to split
   * the input into logical units. The more complex parsing should be
   * done in each individual consumer thread.
   *
   * The parsing thread loops until EOF. It creates an object of type
   * `elt` with `elt_init()`, fills up to `group_size()` element in
   * the `elements` array and update the `nb_filled` member. For
   * example:
   *
   * ~~~{.cc}
   * void parser_loop() {
   *   while(!<EOF>) { // Insert correct test for EOF
   *     elt e(elt_init());
   *     size_type& i = e->nb_filled;
   *     for(i = 0; i < group_size(); ++i) {
   *       // fill element e->elements[i]
   *     }
   *   }
   * }
   * ~~~
   */
  virtual void parser_loop() = 0;

  /** Report an error. The error message will be returned by subsequent call to `error()`. 
   * @param msg The error message.
   */
  void report_error(const char* msg) { jflib::a_store_ptr(error_, msg); }
  typedef typename pool::side side;
  /** Initialization of a group of elements. */
  side& elt_init() { return pool_.get_A(); }

private:
  struct self_pointer {
    multiplexed_parser* self;
  };
  static void* static_start_parser_loop(void* self_);
  void start_parser_loop();
};

template<typename T>
multiplexed_parser<T>::~multiplexed_parser() {
  if(parser_started_) {
    pool_.close_B_to_A();
    // XXX: Do we need to cancel the thread? This seems to leak memory and
    // may not be necessary or desirable.
    //    pthread_cancel(parser_id_);
    pthread_join(parser_id_, 0);
  }
}

template<typename T>
void* multiplexed_parser<T>::static_start_parser_loop(void* self_) {
  std::unique_ptr<self_pointer> self((self_pointer*)self_);
  self->self->start_parser_loop();
  return 0;
}

template<typename T>
void multiplexed_parser<T>::start_parser_loop() {
  parser_loop(); // Call virtual subclass method
  pool_.close_A_to_B();
}

template<typename T>
void multiplexed_parser<T>::start_parsing() {
  // Finish initialization of the read_groups
  for(auto it = pool_.begin(); it != pool_.end(); ++it) {
    it->elements.resize(group_size_);
    it->nb_filled = 0;
  }

  auto self  = new self_pointer;
  self->self = this;
  int  res   = pthread_create(&parser_id_, 0, static_start_parser_loop , (void*)self);
  if(res) {
    delete self;
    throw std::runtime_error("Failed to create the parser thread");
  }
  parser_started_ = true;
}

#endif /* __MULTIPLEXED_PARSER_HPP__ */
