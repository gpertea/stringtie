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


#include <gtest/gtest.h>
#include <unittests/misc.hpp>
#include <multi_thread_skip_list_set.hpp>

#include <string>
#include <vector>

namespace {
  // Compute the power at compile time in a stupid way
  template<long x, long n>
  struct lpow {
    enum { value = x * lpow<x, n-1>::value };
  };
  template<long x>
  struct lpow<x, 0> {
    enum { value = 1 };
  };

  typedef multi_thread_skip_list_set<int> set_type;
  typedef set_type::iterator set_iterator;
  typedef set_type::const_iterator const_set_iterator;
  typedef std::pair<set_iterator, bool> set_ins;

  typedef std::set<int>::iterator std_iterator;
  typedef std::set<int>::const_iterator const_std_iterator;
  typedef std::pair<std_iterator, bool> std_ins;

  TEST(MTSkipListSet, Init) {
    static const int max_height = 5;
    set_type sls(max_height);
    EXPECT_TRUE(sls.empty());
    EXPECT_EQ((size_t)0, sls.size());
    EXPECT_TRUE(sls.end() == sls.begin());
    EXPECT_EQ(((size_t)lpow<set_type::p, max_height>::value), sls.max_size());
  }

  TEST(MTSkipListSet, InsertOneThread) {
    set_type         sls(5, std::less<int>(), xor_random(random()));
    set_type::thread th(sls);
    std::set<int>    set;
  

    for(int i = 0; i < 100; ++i) {
      int     x       = random() % 50;
      std_ins set_res = set.insert(x);
      set_ins sls_res = th.insert(x);
      EXPECT_EQ(set_res.second, sls_res.second);
      EXPECT_EQ(*set_res.first, *sls_res.first);
      EXPECT_EQ(x, *sls_res.first);
      EXPECT_EQ(set.empty(), sls.empty());
      EXPECT_EQ(set.size(), sls.size());
    }

    const_set_iterator set_it = sls.begin();
    for(const_std_iterator std_it = set.begin(); std_it != set.end(); 
        ++std_it, ++set_it) {
      ASSERT_FALSE(set_it == sls.end());
      EXPECT_EQ(*std_it, *set_it);
    }
    ASSERT_TRUE(set_it == sls.end());
  }

  struct thread_insert_data {
    set_type                 set; // The set to add to
    jflib::atomic_field<int> ids; // Counter to get a thread id
    std::vector<int>         v;   // The data to insert
    jflib::atomic_field<int> new_elt; // Nb new element inserted
    jflib::atomic_field<int> exist_elt; // Nb existing element inserted

    static const int         per_th = 10000; // Nb elements per thread
    thread_insert_data(int nb_threads) :
      set(10, std::less<int>(), xor_random(random())), ids(-1),
      new_elt(0), exist_elt(0)
    { 
      const int n = per_th * nb_threads;
      for(int i = 0; i < n; ++i)
        v.push_back(random() % n);
    }
  };
  void* insert_from_thread(void* d) {
    auto data = (thread_insert_data*)d;
    set_type::thread th(data->set);
  
    int tid = (data->ids += 1);
    int new_elt = 0, exist_elt = 0;
    for(int i = tid * data->per_th; i < (tid+1) * data->per_th; ++i) {
      auto res = th.insert(data->v[i]);
      if(res.second)
        ++new_elt;
      else
        ++exist_elt;
    }
    new_elt = (data->new_elt   += new_elt);
    exist_elt = (data->exist_elt += exist_elt);

    return 0;
  }
  TEST(MTSkipListSet, InsertManyThreads) {
    const int nb_threads = 5;
    thread_insert_data data(nb_threads);

    pdo(nb_threads, insert_from_thread, (void*)&data);
    EXPECT_EQ(data.v.size(), (size_t)(data.new_elt + data.exist_elt));
  
    // Do the same single threads into a set and check the statistics
    std::set<int> std_set;
    int new_elt = 0, exist_elt = 0;
    for(auto it = data.v.begin(); it != data.v.end(); ++it) {
      auto res = std_set.insert(*it);
      if(res.second)
        ++new_elt;
      else
        ++exist_elt;
    }
    EXPECT_EQ(new_elt, data.new_elt);
    EXPECT_EQ(exist_elt, data.exist_elt);
    EXPECT_EQ(std_set.size(), data.set.size());

    set_type::thread th(data.set);
    for(int i = -1; i <= (int)data.v.size(); ++i) {
      EXPECT_EQ(std_set.count(i), data.set.count(i));
      { // Test find
        auto std_res = std_set.find(i);
        auto set_res = data.set.find(i);
        if(std_res == std_set.end())
          EXPECT_EQ(data.set.end(), set_res);
        else
          EXPECT_EQ(*std_res, *set_res);
        EXPECT_EQ(set_res, th.find(i));
      }
      { // Test equal_range
        auto std_res = std_set.equal_range(i);
        auto set_res = data.set.equal_range(i);
        auto th_res  = th.equal_range(i);
        EXPECT_EQ(std_res.first == std_res.second, set_res.first == set_res.second);
        EXPECT_EQ(std_res.second == std_set.end(), set_res.second == data.set.end());
        ASSERT_EQ(set_res, th_res);
        if(std_res.first != std_res.second) {
          EXPECT_EQ(i, *set_res.first);
          EXPECT_EQ(set_res.second, ++set_res.first);
        }
      }
      { // Test lower_bound
        auto std_res = std_set.lower_bound(i);
        auto set_res = data.set.lower_bound(i);
        if(std_res == std_set.end())
          EXPECT_EQ(data.set.end(), set_res);
        else
          EXPECT_EQ(*std_res, *set_res);
        EXPECT_EQ(set_res, th.lower_bound(i));
      }
      { // Test upper_bound
        auto std_res = std_set.upper_bound(i);
        auto set_res = data.set.upper_bound(i);
        if(std_res == std_set.end())
          EXPECT_EQ(data.set.end(), set_res);
        else
          EXPECT_EQ(*std_res, *set_res);
        EXPECT_EQ(set_res, th.upper_bound(i));
      }
    }
  }
  // TEST(SkipListSet, InsertIterator) {
  //   skip_list_set<int> sls(5);
  //   std::set<int>      set;
  //   typedef std::set<int>::iterator set_it_type;
  //   typedef skip_list_set<int>::iterator sls_it_type;

  //   std::vector<int> xs;
  //   for(int i = 0; i < 100; ++i)
  //     xs.push_back(random() % 50);
      
  //   set.insert(xs.begin(), xs.end());
  //   sls.insert(xs.begin(), xs.end());
  //   EXPECT_EQ(set.size(), sls.size());
  //   EXPECT_EQ(set.empty(), sls.empty());

  //   sls_it_type sls_it = sls.begin();
  //   for(set_it_type set_it = set.begin(); set_it != set.end(); ++set_it, ++sls_it) {
  //     ASSERT_FALSE(sls.end() == sls_it);
  //     EXPECT_EQ(*set_it, *sls_it);
  //   }
  //   EXPECT_TRUE(sls.end() == sls_it);

  //   for(int i = 0; i < 51; ++i) {
  //     EXPECT_EQ(set.count(i), sls.count(i));
  //     std::pair<set_it_type, set_it_type> set_its = set.equal_range(i);
  //     std::pair<sls_it_type, sls_it_type> sls_its = sls.equal_range(i);
  //     EXPECT_EQ(set_its.first == set_its.second, sls_its.first == sls_its.second);
  //     EXPECT_EQ(set_its.second == set.end(), sls_its.second == sls.end());
  //     if(sls_its.first != sls_its.second) {
  //       EXPECT_EQ(i, *sls_its.first);
  //       EXPECT_EQ(sls_its.second, ++sls_its.first);
  //     }
    
  //     set_it_type set_it = set.find(i);
  //     sls_it_type sls_it = sls.find(i);
  //     if(set_it == set.end()) {
  //       EXPECT_EQ(sls.end(), sls_it);
  //     } else {
  //       ASSERT_NE(sls.end(), sls_it);
  //       EXPECT_EQ(*set_it, *sls_it);
  //     }

  //     set_it = set.lower_bound(i);
  //     sls_it = sls.lower_bound(i);
  //     if(set_it == set.end()) {
  //       EXPECT_EQ(sls.end(), sls_it);
  //     } else {
  //       ASSERT_NE(sls.end(), sls_it);
  //       EXPECT_EQ(*set_it, *sls_it);
  //     }

  //     set_it = set.upper_bound(i);
  //     sls_it = sls.upper_bound(i);
  //     if(set_it == set.end()) {
  //       EXPECT_EQ(sls.end(), sls_it);
  //     } else {
  //       ASSERT_NE(sls.end(), sls_it);
  //       EXPECT_EQ(*set_it, *sls_it);
  //     }
  //   }

  //   set.clear();
  //   sls.clear();
  //   EXPECT_EQ(set.size(), sls.size());
  //   EXPECT_EQ(set.empty(), sls.empty());
  //   EXPECT_EQ(sls.end(), sls.begin());
  // }

  // TEST(SkipListSet, InsertErase) {
  //   skip_list_set<int> sls(5);
  //   std::set<int>      set;
  //   std::vector<int>   xs;

  //   for(int i = 0; i < 100; ++i)
  //     xs.push_back(i * 3);
  //   random_shuffle(xs.begin(), xs.end());

  //   set.insert(xs.begin(), xs.end());
  //   sls.insert(xs.begin(), xs.end());
  //   EXPECT_TRUE(std::equal(set.begin(), set.end(), sls.begin()));
  
  //   random_shuffle(xs.begin(), xs.end());
  //   for(std::vector<int>::iterator it = xs.begin(); it != xs.end(); ++it) {
  //     size_t set_n = set.erase(*it);
  //     size_t sls_n = sls.erase(*it);
  //     EXPECT_EQ(set_n, sls_n);
  //     EXPECT_EQ(set.size(), sls.size());
  //     EXPECT_EQ(set.empty(), sls.empty());
  //     EXPECT_TRUE(std::equal(set.begin(), set.end(), sls.begin()));
  //   }
  //   EXPECT_TRUE(sls.empty());
  // }

  // TEST(SkipListSet, Copy) {
  //   skip_list_set<std::string> sls(7);
  //   std::string str(10, 'A');
  
  //   static const size_t nb_elts = 100;
  //   for(size_t i = 0; i < nb_elts; ++i) {
  //     // Create random string
  //     for(size_t j = 0; j < str.size(); ++j)
  //       str[j] = 'A' + random() % 26;
  //     sls.insert(str);
  //   }
  //   // This could fail with very small probability, if two identical
  //   // random string are generated.
  //   EXPECT_EQ(nb_elts, sls.size());
  
  //   skip_list_set<std::string> nsls(sls);
  //   EXPECT_EQ(sls.size(), nsls.size());
  //   EXPECT_TRUE(std::equal(sls.begin(), sls.end(), nsls.begin()));

  //   skip_list_set<std::string> msls;
  //   EXPECT_TRUE(msls.empty());
  //   msls = sls;
  //   EXPECT_EQ(sls.size(), nsls.size());
  //   EXPECT_TRUE(std::equal(sls.begin(), sls.end(), nsls.begin()));
  // }

} // namespace
