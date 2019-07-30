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
#include <exp_buffer.hpp>
#include <algorithm>
#include <misc.hpp>

template<typename T>
class ExpandingBufferInit : public ::testing::Test { };
template<typename T>
class ExpandingBufferDefault : public ::testing::Test { };

typedef ExpandingBuffer<int, reallocator<int> > int_realloc_buf;
typedef ExpandingBuffer<int, reallocator_init<int, -1> > init_realloc_buf;
typedef ExpandingBuffer<int, remaper<int> > int_remap_buf;
typedef ExpandingBuffer<int, remaper_init<int, -1> > init_remap_buf;

typedef ::testing::Types<int_realloc_buf, int_remap_buf> ExpBufferDefaultTypes;
typedef ::testing::Types<init_realloc_buf, init_remap_buf> ExpBufferInitTypes;

TYPED_TEST_CASE(ExpandingBufferDefault, ExpBufferDefaultTypes);
TYPED_TEST_CASE(ExpandingBufferInit, ExpBufferInitTypes);

TYPED_TEST(ExpandingBufferDefault, Initialization) {
  TypeParam b;

  EXPECT_EQ((size_t)0, b.capacity());
  EXPECT_EQ((size_t)0, b.size());

  b[5] = 5;
  EXPECT_EQ((size_t)6, b.capacity());
  EXPECT_EQ((size_t)6, b.size());
  // This holds only for remaper
  // for(int i = 0; i < 5; ++i)
  //   EXPECT_EQ(0, b[i]);
  EXPECT_EQ(5, b[5]);
  b[3] = 3;
  EXPECT_EQ((size_t)6, b.capacity());
  EXPECT_EQ((size_t)6, b.size());
  EXPECT_EQ(3, b[3]);

  b[6] = 6;
  EXPECT_EQ((size_t)12, b.capacity());
  EXPECT_EQ((size_t)7, b.size());

  b[5000] = 5000;
  EXPECT_EQ((size_t)5001, b.size());
  EXPECT_EQ(5000, b.back());
  // Again, only for remaper
  // for(typename TypeParam::iterator it = b.begin(); it != b.end(); ++it)
  //   EXPECT_TRUE(*it == 0 || *it == (it - b.begin()));
}

TYPED_TEST(ExpandingBufferDefault, push_back) {
  TypeParam b;

  EXPECT_TRUE(b.empty());

  for(int i = 0; i < 100; ++i) {
    b.push_back(i);
    EXPECT_FALSE(b.empty());
    EXPECT_EQ((size_t)(i+1), b.size());
    EXPECT_EQ(i, b.back());
  }
  for(int i = 0; i < 100; ++i) {
    EXPECT_EQ(i, b[i]);
  }
}

TYPED_TEST(ExpandingBufferDefault, Swap) {
  TypeParam b(10);

  for(size_t i = 0; i < b.capacity(); ++i)
    b[i] = 2 * i;

  EXPECT_EQ((size_t)10, b.size());

  TypeParam bs;
  b.swap(bs);
  EXPECT_EQ((size_t)0, b.capacity());
  EXPECT_EQ((size_t)10, bs.capacity());
  for(size_t i = 0; i < bs.capacity(); ++i)
    EXPECT_EQ((int)(2 * i), bs[i]);

  std::swap(b, bs);

  EXPECT_EQ((size_t)10, b.capacity());
  EXPECT_EQ((size_t)0, bs.capacity());
  for(size_t i = 0; i < b.capacity(); ++i)
    EXPECT_EQ((int)(2 * i), b[i]);

  ExpBuffer<int> b1, b2(5);
  std::swap(b1, b2);
  EXPECT_EQ((size_t)5, b1.capacity());
  EXPECT_EQ((size_t)0, b2.capacity());
}

TYPED_TEST(ExpandingBufferDefault, ConsDes) {
  TypeParam b(5);

  EXPECT_EQ((size_t)5, b.capacity());
}

template<typename T>
T fill_buffer(bool& parity) {
  T even(10);
  T odd(10);

  for(int i = 0; i < 10; ++i) {
    even[i] = 2 * i;
    odd[i] = 2 * i + 1;
  }

  parity = random_bits(1) == 1;
  return parity ? even : odd;
}
template<typename T>
bool test_parity(const T& b, const bool parity) {
  int add = parity ? 0 : 1;
  for(size_t i = 0; i < b.size(); ++i)
    if(b[i] != 2 * (int)i + add)
        return false;

  return true;
}
template<typename T>
bool test_deep_equal(const T& a, const T& b) {
  if(a.size() != b.size())
    return false;

  for(size_t i = 0; i < a.size(); ++i)
    if(a[i] != b[i])
      return false;

  return true;
}
TYPED_TEST(ExpandingBufferDefault, CopyMove) {
  bool parity;
  TypeParam b(fill_buffer<TypeParam>(parity)); // Call move constructor, or RVO

  EXPECT_TRUE(test_parity(b, parity));

  b = fill_buffer<TypeParam>(parity); // Call move assignment operator
  EXPECT_TRUE(test_parity(b, parity));

  TypeParam b1(b); // Call copy constructor
  EXPECT_TRUE(test_deep_equal(b, b1));

  TypeParam b2(fill_buffer<TypeParam>(parity));
  b1 = b2; // Call copy assignment operator
  EXPECT_TRUE(test_deep_equal(b2, b1));

  TypeParam b3(std::move(b1)); // Call move constructor b1 not valid anymore (but b2 == former b1)
  EXPECT_TRUE(test_deep_equal(b2, b3));
}


TYPED_TEST(ExpandingBufferInit, Initialization) {
  TypeParam b;
  EXPECT_EQ((size_t)0, b.capacity());
  EXPECT_EQ((size_t)0, b.size());
  b[10] = 5;
  EXPECT_EQ(b.size() - 1, (size_t)std::count(b.begin(), b.end(), -1));
}

TYPED_TEST(ExpandingBufferDefault, Resize) {
  TypeParam b;

  EXPECT_EQ((size_t)0, b.size());
  b.resize(10);
  EXPECT_EQ((size_t)10, b.size());
  b.resize(5);
  EXPECT_EQ((size_t)5, b.size());
  b.resize(15, 2345);
  for(auto it = b.begin() + 5; it != b.end(); ++it)
    EXPECT_EQ(2345, *it);
}

TYPED_TEST(ExpandingBufferInit, Resize) {
  TypeParam b;

  EXPECT_EQ((size_t)0, b.size());
  b.resize(10);
  EXPECT_EQ((size_t)10, b.size());
  for(auto it = b.begin(); it != b.end(); ++it)
    EXPECT_EQ(-1, *it);

  b.resize(5);
  EXPECT_EQ((size_t)5, b.size());

  b.resize(15, 2345);
  for(auto it = b.begin(); it != b.begin() + 5; ++it)
    EXPECT_EQ(-1, *it);
  for(auto it = b.begin() + 5; it != b.end(); ++it)
    EXPECT_EQ(2345, *it);
}
