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


#include<heap.hpp>
#include<gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <exp_buffer.hpp>

static const size_t nb_elts = 1000;

TEST(HeapMax, PushPop) {
  heap<int>::max h;
  std::vector<int> elts;

  for(size_t i = 0; i < nb_elts; ++i)
    elts.push_back((int)i);
  std::random_shuffle(elts.begin(), elts.end());

  int max = std::numeric_limits<int>::min();
  for(int i = 0; i < (int)elts.size(); ++i) {
    max = std::max(max, elts[i]);
    h.push(elts[i]);
    EXPECT_EQ(i+1, (int)h.size());
    EXPECT_EQ(max, h.peek());
  }
  EXPECT_EQ(nb_elts, elts.size());
  for(int i = 0; i < (int)elts.size(); ++i)
    EXPECT_EQ((int)elts.size() - i - 1, h.pop());

  h.clear();
  EXPECT_EQ((size_t)0, h.size());
}

template<typename heap_T>
class HeapMin : public ::testing::Test { };

typedef ::testing::Types<heap<int>::min, heap<int,ExpBuffer<int> >::min> MyTypes;
TYPED_TEST_CASE(HeapMin, MyTypes);

TYPED_TEST(HeapMin, PushPop) {
  TypeParam h;
  std::vector<int> elts;

  for(size_t i = 0; i < nb_elts; ++i)
    elts.push_back((int)i);
  std::random_shuffle(elts.begin(), elts.end());

  int min = std::numeric_limits<int>::max();
  for(int i = 0; i < (int)elts.size(); ++i) {
    min = std::min(min, elts[i]);
    h.push(elts[i]);
    EXPECT_EQ(i+1, (int)h.size());
    EXPECT_EQ(min, h.peek());
  }
  EXPECT_EQ(nb_elts, elts.size());

  for(int i = 0; i < (int)elts.size(); ++i)
    EXPECT_EQ(i, h.pop());

  h.clear();
  EXPECT_EQ((size_t)0, h.size());
}

TYPED_TEST(HeapMin, Heapify) {
  std::vector<int> elts;

  for(int i = 0; i < 10; ++i)
    elts.push_back(i);
  std::random_shuffle(elts.begin(), elts.end());

  heap<int>::min h(elts.begin(), elts.end());
  EXPECT_EQ(elts.size(), h.size());
  for(int i = 0; i < (int)elts.size(); ++i)
    EXPECT_EQ(i, h.pop());
}

TYPED_TEST(HeapMin, Copy) {
  std::vector<int> elts;

  for(int i = 0; i < 10; ++i)
    elts.push_back(i);
  std::random_shuffle(elts.begin(), elts.end());

  heap<int>::min h1(elts.begin(), elts.end());
  heap<int>::min h2;

  EXPECT_EQ(elts.size(), h1.size());
  EXPECT_EQ((size_t)0, h2.size());

  std::swap(h1, h2);
  EXPECT_EQ(elts.size(), h2.size());
  EXPECT_EQ((size_t)0, h1.size());

  heap<int>::min h3(h2);
  EXPECT_EQ(elts.size(), h2.size());
  EXPECT_EQ(elts.size(), h3.size());
  
  heap<int>::min h4;
  h4 = h2;
  EXPECT_EQ(elts.size(), h2.size());
  EXPECT_EQ(elts.size(), h4.size());
}
