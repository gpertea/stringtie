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


#include <iostream>
#include <algorithm>
#include <gtest/gtest.h>
#include <rb_tree.hpp>
#include <stdlib.h>

typedef std::vector<int> int_vec;
typedef red_black::basic_tree< int, std::less<int>, red_black::new_delete_allocator<int> > int_tree_new;
typedef red_black::basic_tree< int, std::less<int>, red_black::container_allocator<int, std::vector<red_black::node<int, size_t> > > > int_tree_vec;

long my_rand(long m) { return random() % m; }

template<typename T>
class RBTreeTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srandom(time(0));
    for(int i = 0; i < 1024; ++i)
      sorted_elts.push_back(i);
    elts = sorted_elts;
    std::random_shuffle(elts.begin(), elts.end(), my_rand);
  }

  int_vec sorted_elts;
  int_vec elts;
};

typedef ::testing::Types<int_tree_new, int_tree_vec> TreeTypes;
TYPED_TEST_CASE(RBTreeTest, TreeTypes);

TYPED_TEST(RBTreeTest, RandomOrder) {
  int_vec&  elts        = TestFixture::elts;
  int_vec&  sorted_elts = TestFixture::sorted_elts;
  TypeParam tree;

  EXPECT_TRUE(tree.empty());
  EXPECT_EQ((typename TypeParam::size_type)0, tree.size());

  tree.insert(elts.begin(), elts.end());
  EXPECT_FALSE(tree.empty());
  EXPECT_EQ((typename TypeParam::size_type)elts.size(), tree.size());
  typename TypeParam::iterator it_tree = tree.begin();
  int_vec::iterator            it_vec  = sorted_elts.begin();

  for( ; it_tree != tree.end() && it_vec != sorted_elts.end(); ++it_tree, ++it_vec)
    EXPECT_EQ(*it_vec, *it_tree);
  EXPECT_EQ(tree.end(), it_tree);
  EXPECT_EQ(sorted_elts.end(), it_vec);
}

TYPED_TEST(RBTreeTest, Consistency) {
  int_vec&  elts        = TestFixture::elts;
  TypeParam tree;
  EXPECT_TRUE(tree.is_consistent());
  for(int_vec::iterator it_vec = elts.begin(); it_vec != elts.end(); ++it_vec) {
    tree.insert(*it_vec);
    EXPECT_TRUE(tree.is_consistent());
  }

  std::random_shuffle(elts.begin(), elts.end(), my_rand);
  size_t len = elts.size();
  for(int_vec::iterator it_vec = elts.begin(); it_vec != elts.end(); ++it_vec, --len) {
    tree.remove(*it_vec);
    EXPECT_TRUE(tree.is_consistent());
    EXPECT_EQ(len - 1, tree.size());
  }
}
