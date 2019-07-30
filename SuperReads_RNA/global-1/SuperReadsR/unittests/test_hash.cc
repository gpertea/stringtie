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
#include <src/MurmurHash3.h>

class rand_block_string {
  std::string str_;

public:
  rand_block_string(int length) : str_(length, '\0') {
    for(int i = 0; i != length; ++i)
      str_[i] = (char)((random() % 255) + 1);
  }

  // Subscript by 64 bit word
  uint64_t operator[](unsigned int i) const {
    uint64_t* ptr = (uint64_t*)str_.c_str() + i;
    return *ptr;
  }
  const char* raw() const { return str_.c_str(); }
  size_t len() const { return str_.size(); }
};

TEST(MurmurHash3, Generic) {
  rand_block_string str(random() % 200);
  uint32_t seed = random();
  uint64_t original[2];
  uint64_t generic[2];
  
  MurmurHash3_x64_128(str.raw(), str.len(), seed, original);
  MurmurHash3_T_128(str, str.len(), seed, generic);
  EXPECT_EQ(original[0], generic[0]);
  EXPECT_EQ(original[1], generic[1]);
}
