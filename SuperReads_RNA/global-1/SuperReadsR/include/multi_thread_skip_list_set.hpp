#ifndef __MULTI_THREAD_SKIP_LIST_SET_HPP__
#define __MULTI_THREAD_SKIP_LIST_SET_HPP__

#include <algorithm>
#include <utility>
#include <functional>
#include <jellyfish/atomic_field.hpp>
#include <jellyfish/compare_and_swap.hpp>
#include <jellyfish/locks_pthread.hpp>
#include <skip_list_common.hpp>

/** Multi-threaded lock free set based on skip list. In multi-threaded
 * mode, it supports only adding an element to the set and testing for
 * membership. Only the member functions of the class
 * multi_thread_skip_list_set which marked "const" are thread
 * safe. Use the method of the class
 * multi_thread_skip_list_set::thread to insert elements. It is
 * constructed from a multi_thread_skip_list_set class.
 */

/* TODO: There is a lot of code duplication with
   skip_list_set.hpp. Can it be reduced while still being readable?
 */

template <typename Key, typename Compare = std::less<Key>, int p_ = 4, typename Random = xor_random>
class multi_thread_skip_list_set {
protected:
  struct node {
    Key   k;
    int   height;
    node* tower[];
  };
  struct path_node {
    node** ptr;
    node*  val;
  };

  node**                           heads_;
  int                              max_height_;
  Compare                          comp_;
  Random                           rand_;
  jellyfish::locks::pthread::mutex seed_mutex;

public:
  typedef Key                                   key_type;
  typedef Key                                   value_type;
  typedef Compare                               key_compare;
  typedef Compare                               value_compare;
  typedef Key&                                  reference;
  typedef const Key&                            const_reference;
  typedef size_t                                size_type;
  typedef ssize_t                               difference_type;
  typedef Key*                                  pointer;
  typedef const Key*                            const_pointer;
  static const int p = p_;

  class node_iterator {
  protected:
    node* item;
    node_iterator(node* item_) : item(item_) { }
    node_iterator(const node_iterator& rhs) : item(rhs.item) { }
    void next() { item = jflib::a_load(item->tower); }
  public:
    bool operator==(const node_iterator& rhs) const { return item == rhs.item; }
    bool operator!=(const node_iterator& rhs) const { return item != rhs.item; }
  };

  class iterator : 
    public std::iterator<std::forward_iterator_tag, key_type>,
    public node_iterator {
    friend class multi_thread_skip_list_set;
    iterator(node* item_) : node_iterator(item_) { }
  public:
    iterator() : node_iterator(0) { }
    iterator(const node_iterator& rhs) : node_iterator(rhs) { }

    iterator& operator=(iterator rhs) {
      std::swap(node_iterator::item, rhs.item);
      return *this;
    }    
    reference operator*() { return node_iterator::item->k; }
    pointer operator->() { return &node_iterator::item->k; }
    iterator& operator++() {
      node_iterator::next();
      return *this;
    }
    iterator operator++(int) {
      iterator c(*this);
      node_iterator::next();
      return c;
    }
  };
  class const_iterator :
    public std::iterator<std::forward_iterator_tag, key_type>,
    public node_iterator {
    friend class multi_thread_skip_list_set;
    const_iterator(node* item_) : node_iterator(item_) { }
  public:
    const_iterator() : node_iterator(0) { }
    const_iterator(const const_iterator& rhs) : node_iterator(rhs) { }
    const_iterator(const iterator& rhs) : node_iterator(rhs) { }

    const_iterator& operator=(node_iterator rhs) {
      swap(node_iterator::item, rhs.item);
      return *this;
    }    
    const_reference operator*() { return node_iterator::item->k; }
    const_pointer operator->() { return &node_iterator::item->k; }
    const_iterator& operator++() {
      node_iterator::next();
      return *this;
    }
    const_iterator operator++(int) {
      const_iterator c(*this);
      node_iterator::next();
      return c;
    }
  };


  explicit multi_thread_skip_list_set(int max_height = 10,
                                      const Compare& comp = Compare(), 
                                      const Random& rand = Random()) :
    heads_(new node*[max_height]), max_height_(max_height),
    comp_(comp), rand_(rand)
  {
    memset(heads_, '\0', sizeof(node*) * max_height_);
  }

  // template<class InputIterator>
  // multi_thread_skip_list_set(int max_height, InputIterator first, InputIterator last,
  //                            const Compare& comp = Compare(), 
  //                            const Random& rand = Random()) :
  //   heads_(new node*[max_height]),
  //   max_height_(max_height), cur_height_(1), size_(0),
  //   comp_(comp), rh_(rand)
  // {
  //   memset(heads_, '\0', sizeof(node*) * max_height_);    
  //   insert(first, last);
  // }    

  virtual ~multi_thread_skip_list_set() {
    clear();
    delete [] heads_;
  }

  multi_thread_skip_list_set& operator=(multi_thread_skip_list_set rhs) {
    swap(rhs);
    return *this;
  }

  void swap(multi_thread_skip_list_set& rhs) {
    std::swap(heads_, rhs.heads_);
    std::swap(max_height_, rhs.max_height_);
    std::swap(comp_, rhs.comp_);
    std::swap(rand_, rhs.rand_);
  }

  void clear() {
    node* cnode = heads_[0];
    while(cnode) {
      node* nnode = cnode->tower[0];
      delete cnode;
      cnode = nnode;
    }
    memset(heads_, '\0', sizeof(node*) * max_height_);
  }

  /* The following methods are thread safe.
   */
  size_t size() const {
    size_t res = 0;
    for(node* ptr = jflib::a_load(heads_); ptr; ptr = jflib::a_load(ptr->tower))
      ++res;
    return res;
  }
  bool empty() const { return jflib::a_load(heads_) == (node*)0; }
  size_type max_size() const {
    size_type res = 1;
    size_type pp  = p;
    int       n   = max_height_;

    while(n) {
      if(n & 0x1)
        res *= pp;
      pp *= pp;
      n >>= 1;
    }
    return res;
  }
  iterator begin() { return iterator(jflib::a_load(heads_)); }
  const_iterator begin() const { return const_iterator(jflib::a_load(heads_)); }
  iterator end() { return iterator(); }
  const_iterator end() const { return const_iterator(); }

  template<typename T>
  size_type count(const T& x) const {
    return find_node(x) ? 1 : 0;
  }
  template<typename T>
  std::pair<iterator, iterator> equal_range(const T& x) const {
    static_assert(std::is_pod<path_node>::value, "Path_node is a POD");
    path_node path[max_height_];
    node* n = find_node_path(x, path);
    if(n)
      return std::make_pair(iterator(n), ++iterator(n));
    return std::make_pair(iterator(path[0].val), iterator(path[0].val));
  }

  template<typename T>
  iterator find(const T& x) const {
    return iterator(find_node(x));
  }
  template<typename T>
  iterator lower_bound(const T& x) const {
    static_assert(std::is_pod<path_node>::value, "Path_node is a POD");
    path_node path[max_height_];
    find_node_path(x, path);
    return iterator(path[0].val);
  }
  template<typename T>
  iterator upper_bound(const T& x) const {
    static_assert(std::is_pod<path_node>::value, "Path_node is a POD");
    path_node path[max_height_];
    node*     n = find_node_path(x, path);
    if(n)
      return ++iterator(n);
    return iterator(path[0].val);
  }
  

  friend class thread;
  class thread {
    multi_thread_skip_list_set& set_;
    random_height<Random, p_>   rh_;

  public:
    thread(multi_thread_skip_list_set& set) : 
      set_(set), rh_(Random(set.new_seed())) { }
    virtual ~thread() { }

    bool empty() const { return set_.empty(); }
    size_type max_size() const { return set_.max_size(); }
    template<typename T>
    size_type count(const T& x) const { return set_.count(x); }
    template<typename T>
    std::pair<iterator, iterator> equal_range(const T& x) const { 
      return set_.equal_range(x);
    }
    template<typename T>
    iterator find(const T& x) const { return set_.find(x); }
    template<typename T>
    iterator lower_bound(const T& x) const { return set_.lower_bound(x); }
    template<typename T>
    iterator upper_bound(const T& x) const { return set_.upper_bound(x); }

    std::pair<iterator, bool> insert(const value_type& x) {
      static_assert(std::is_pod<path_node>::value, "Path_node is a POD");
      path_node path[set_.max_height_];
      node* n = set_.find_node_path(x, path);
      if(n)
        return std::make_pair(iterator(n), false);
      n = new_node(x);
      for(int i = 0; i < n->height; ++i) {
        node*  cval = path[i].val;
        node** cptr = path[i].ptr;
        n->tower[i] = cval;
        node*  oval = cval;
        while(!jflib::cas(cptr, oval, n, &cval)) {
          if(set_.comp_(cval->k, x)) {
            cptr = &(cval->tower[i]);
          } else if(set_.comp_(x, cval->k)) {
            n->tower[i] = cval;
            oval        = cval;
          } else {
            delete n;
            return std::make_pair(iterator(cval), false);
          }
        }
      }
      return std::make_pair(iterator(n), true);
    }

  private:
    // Allocate a new node. Does raw memory allocation of a node with
    // enough space for the tower. Then in place copy construction of
    // the key from x.
    node* new_node(const value_type& x) {
      //    int height  = std::min(max_height_, std::min(cur_height_ + 1, rh_()));
      int height = std::min(set_.max_height_, rh_());
      node* res   = (node*)operator new(sizeof(node) + height * sizeof(node*));
      res->height = height;
      new ((void*)&res->k) value_type(x);
      return res;
    }
  };

private:
  // Find the path to a node equal to x
  template<typename T>
  node* find_node_path(const T& x, path_node path[]) const {
    int   i     = max_height_ - 1;
    node* cnode = 0;
    node* nnode = 0;

    while(i >= 0) {
      cnode = jflib::a_load(&heads_[i]);
      if(cnode && comp_(cnode->k, x))
        break;
      path[i].ptr = &heads_[i];
      path[i].val = cnode;
      --i;
    }
    if(i < 0)
      return cnode && !comp_(x, cnode->k) ? cnode : 0;

    for( ; i >= 0; --i) {
      nnode = jflib::a_load(&cnode->tower[i]);
      while(nnode && comp_(nnode->k, x)) {
        cnode = nnode;
        nnode = jflib::a_load(&nnode->tower[i]);
      }
      path[i].ptr = &(cnode->tower[i]);
      path[i].val = nnode;
    }

    return nnode && !comp_(x, nnode->k) ? nnode : 0;
  }
  
  // Find a node equal to x. The path is not recorded.
  template<typename T>
  node* find_node(const T& x) const {
    int   i     = max_height_ - 1;
    node* cnode = 0;
    node* nnode = 0;
    
    for( ; i >=0; --i) {
      cnode = jflib::a_load(&heads_[i]);
      if(cnode && comp_(cnode->k, x))
        break;
    }
    if(i < 0)
      return cnode && !comp_(x, cnode->k) ? cnode : 0;

    for( ; i >= 0; --i) {
      nnode = jflib::a_load(&cnode->tower[i]);
      while(nnode && comp_(nnode->k, x)) {
        cnode = nnode;
        nnode = jflib::a_load(&cnode->tower[i]);
      }
    }
    return nnode && !comp_(x, nnode->k) ? nnode : 0;
  }

  uint64_t new_seed() {
    jellyfish::locks::pthread::mutex_lock ml(seed_mutex);
    return rand_();
  }
};

#endif /* __MULTI_THREAD_SKIP_LIST_SET_HPP__ */
