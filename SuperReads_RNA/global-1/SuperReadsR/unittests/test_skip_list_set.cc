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
#include <skip_list_set.hpp>

#include <string>
#include <vector>

// Compute the power at compile time in a stupid way
template<long x, long n>
struct lpow {
  enum { value = x * lpow<x, n-1>::value };
};
template<long x>
struct lpow<x, 0> {
  enum { value = 1 };
};

TEST(SkipListSet, Init) {
  static const int max_height = 5;
  skip_list_set<int> sls(max_height);
  EXPECT_EQ((size_t)0, sls.size());
  EXPECT_TRUE(sls.empty());
  EXPECT_TRUE(sls.end() == sls.begin());
  EXPECT_EQ(((size_t)lpow<skip_list_set<int>::p, max_height>::value), sls.max_size());
}

TEST(SkipListSet, Insert) {
  skip_list_set<int> sls(5);
  std::set<int>      set;
  typedef std::pair<std::set<int>::iterator, bool> set_ins;
  typedef std::pair<skip_list_set<int>::iterator, bool> sls_ins;

  for(int i = 0; i < 100; ++i) {
    int     x       = random() % 50;
    set_ins set_res = set.insert(x);
    sls_ins sls_res = sls.insert(x);
    EXPECT_EQ(set_res.second, sls_res.second);
    EXPECT_EQ(*set_res.first, *sls_res.first);
    EXPECT_EQ(x, *sls_res.first);
    EXPECT_EQ(set.size(), sls.size());
    EXPECT_EQ(set.empty(), sls.empty());
  }
}

TEST(SkipListSet, InsertIterator) {
  skip_list_set<int> sls(5);
  std::set<int>      set;
  typedef std::set<int>::iterator set_it_type;
  typedef skip_list_set<int>::iterator sls_it_type;

  std::vector<int> xs;
  for(int i = 0; i < 100; ++i)
    xs.push_back(random() % 50);
      
  set.insert(xs.begin(), xs.end());
  sls.insert(xs.begin(), xs.end());
  EXPECT_EQ(set.size(), sls.size());
  EXPECT_EQ(set.empty(), sls.empty());

  sls_it_type sls_it = sls.begin();
  for(set_it_type set_it = set.begin(); set_it != set.end(); ++set_it, ++sls_it) {
    ASSERT_FALSE(sls.end() == sls_it);
    EXPECT_EQ(*set_it, *sls_it);
  }
  EXPECT_TRUE(sls.end() == sls_it);

  for(int i = 0; i < 51; ++i) {
    EXPECT_EQ(set.count(i), sls.count(i));
    std::pair<set_it_type, set_it_type> set_its = set.equal_range(i);
    std::pair<sls_it_type, sls_it_type> sls_its = sls.equal_range(i);
    EXPECT_EQ(set_its.first == set_its.second, sls_its.first == sls_its.second);
    EXPECT_EQ(set_its.second == set.end(), sls_its.second == sls.end());
    if(sls_its.first != sls_its.second) {
      EXPECT_EQ(i, *sls_its.first);
      EXPECT_EQ(sls_its.second, ++sls_its.first);
    }
    
    set_it_type set_it = set.find(i);
    sls_it_type sls_it = sls.find(i);
    if(set_it == set.end()) {
      EXPECT_EQ(sls.end(), sls_it);
    } else {
      ASSERT_NE(sls.end(), sls_it);
      EXPECT_EQ(*set_it, *sls_it);
    }

    set_it = set.lower_bound(i);
    sls_it = sls.lower_bound(i);
    if(set_it == set.end()) {
      EXPECT_EQ(sls.end(), sls_it);
    } else {
      ASSERT_NE(sls.end(), sls_it);
      EXPECT_EQ(*set_it, *sls_it);
    }

    set_it = set.upper_bound(i);
    sls_it = sls.upper_bound(i);
    if(set_it == set.end()) {
      EXPECT_EQ(sls.end(), sls_it);
    } else {
      ASSERT_NE(sls.end(), sls_it);
      EXPECT_EQ(*set_it, *sls_it);
    }
  }

  set.clear();
  sls.clear();
  EXPECT_EQ(set.size(), sls.size());
  EXPECT_EQ(set.empty(), sls.empty());
  EXPECT_EQ(sls.end(), sls.begin());
}

TEST(SkipListSet, InsertErase) {
  skip_list_set<int> sls(5);
  std::set<int>      set;
  std::vector<int>   xs;

  for(int i = 0; i < 100; ++i)
    xs.push_back(i * 3);
  random_shuffle(xs.begin(), xs.end());

  set.insert(xs.begin(), xs.end());
  sls.insert(xs.begin(), xs.end());
  EXPECT_TRUE(std::equal(set.begin(), set.end(), sls.begin()));
  
  random_shuffle(xs.begin(), xs.end());
  for(std::vector<int>::iterator it = xs.begin(); it != xs.end(); ++it) {
    size_t set_n = set.erase(*it);
    size_t sls_n = sls.erase(*it);
    EXPECT_EQ(set_n, sls_n);
    EXPECT_EQ(set.size(), sls.size());
    EXPECT_EQ(set.empty(), sls.empty());
    EXPECT_TRUE(std::equal(set.begin(), set.end(), sls.begin()));
  }
  EXPECT_TRUE(sls.empty());
}

TEST(SkipListSet, Copy) {
  skip_list_set<std::string> sls(7);
  std::string str(10, 'A');
  
  static const size_t nb_elts = 100;
  for(size_t i = 0; i < nb_elts; ++i) {
    // Create random string
    for(size_t j = 0; j < str.size(); ++j)
      str[j] = 'A' + random() % 26;
    sls.insert(str);
  }
  // This could fail with very small probability, if two identical
  // random string are generated.
  EXPECT_EQ(nb_elts, sls.size());
  
  skip_list_set<std::string> nsls(sls);
  EXPECT_EQ(sls.size(), nsls.size());
  EXPECT_TRUE(std::equal(sls.begin(), sls.end(), nsls.begin()));

  skip_list_set<std::string> msls;
  EXPECT_TRUE(msls.empty());
  msls = sls;
  EXPECT_EQ(sls.size(), nsls.size());
  EXPECT_TRUE(std::equal(sls.begin(), sls.end(), nsls.begin()));
}
