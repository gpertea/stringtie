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
#include <exp_vector.hpp>

TEST(ExpVector, Subscript) {
  exp_vector<int> ev;
  exp_vector<int> evd(5);

  EXPECT_EQ(int(), ev.default_value());
  EXPECT_EQ(5, evd.default_value());

  ev[10]  = 5;
  evd[10] = 10;
  for(int i = 0; i < 10; ++i) {
    EXPECT_EQ(int(), ev[i]);
    EXPECT_EQ(5, evd[i]);
  }
  EXPECT_EQ(5, ev[10]);
  EXPECT_EQ(10, evd[10]);
}

TEST(ExpVector, Float) {
  exp_vector<float> fv;
  exp_vector<float> fvd(0, 5.0);

  fv[10] = 10.0;
  fvd[10] = 10.0;

  for(int i = 0; i < 10; ++i) {
    EXPECT_FLOAT_EQ(float(), fv[i]);
    EXPECT_FLOAT_EQ((float)5.0, fvd[i]);
  }

  EXPECT_FLOAT_EQ((float)10.0, fv[10]);
  EXPECT_FLOAT_EQ((float)10.0, fvd[10]);
}

TEST(ExpVector, Initialization) {
  static const size_t size = 100;
  exp_vector<int> iv(size);

  EXPECT_EQ(size, iv.size());
  for(size_t i = 0; i < size; ++i)
    EXPECT_EQ(0, iv[i]);
}
