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


#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <string>
#include <gtest/gtest.h>
#include <charb.hpp>
#include <src/bloom_counter2.hpp>

typedef std::pair<std::string,unsigned int> elt;
struct parameter {
  size_t nb_strings;
  size_t nb_inserts;
  double error_rate;
};
class BloomCounter : public ::testing::TestWithParam<parameter> { };

TEST_P(BloomCounter, FalsePositive) {
  const size_t nb_inserts = GetParam().nb_inserts;
  const double error_rate = GetParam().error_rate;
  const size_t nb_strings = GetParam().nb_strings;
  static const size_t str_len    = 100;
  std::vector<elt>    counts;
  charb               str(str_len);
  str[str_len] = '\0';

  for(size_t i = 0; i < nb_strings; ++i) {
    for(size_t j = 0; j < str_len; ++j)
      str[j] = 'A' + (random() % 26);
    counts.push_back(std::make_pair((char*)str, (unsigned int)0));
  }

  bloom_counter2<const char *> bc1(error_rate, nb_inserts);
  bloom_counter2<const char *> bc2(error_rate, nb_inserts);
  bloom_counter2<const char *> bc3(error_rate, nb_inserts);

  size_t nb_errors   = 0;
  size_t collisions2 = 0;
  size_t collisions3 = 0;

  for(size_t i = 0; i < nb_inserts; ++i) {
    std::vector<elt>::reference ref = counts[random() % nb_strings];
    ++ref.second;
    nb_errors += bc1.insert(ref.first.c_str()) >= ref.second;
    bc2[ref.first.c_str()]++;
    ++bc3[ref.first.c_str()];
  }

  // Check known strings
  for(size_t i = 0; i < nb_strings; ++i) {
    std::vector<elt>::reference ref = counts[i];
    unsigned int expected = std::min(ref.second, (unsigned int)2);
    unsigned int actual = bc1.check(ref.first.c_str());
    EXPECT_LE(expected, actual);
    if(expected != actual)
      ++nb_errors;
    if(expected != *bc2[ref.first.c_str()])
      ++collisions2;
    if(expected != *bc3[ref.first.c_str()])
      ++collisions3;
  }
  EXPECT_GT(error_rate * nb_strings, nb_errors);
  EXPECT_GT(error_rate * nb_strings, collisions2);
  EXPECT_GT(error_rate * nb_strings, collisions3);

  nb_errors = collisions2 = collisions3 = 0;
  // Check unknown strings
  for(size_t i = 0; i < nb_inserts; ++i) {
    for(size_t j = 0; j < str_len; ++j)
      str[j] = '0' + (random() % 10);
    if(bc1.check(str) > 0)      
      ++nb_errors;
    if(bc2[str] > 0)
      ++collisions2;
    if(bc3[str] > 0)
      ++collisions3;
  }
  EXPECT_GT(2 * error_rate * nb_inserts, nb_errors);
  EXPECT_GT(2 * error_rate * nb_inserts, collisions2);
  EXPECT_GT(2 * error_rate * nb_inserts, collisions3);
}

INSTANTIATE_TEST_CASE_P(BloomCounterTest, BloomCounter,
                        ::testing::Values(parameter({    1024,    2048, 0.01 }),
                                          //                                          parameter({ 1000000, 3000000, 0.01 })));
                                          parameter({   10000,   30000, 0.01 })));
