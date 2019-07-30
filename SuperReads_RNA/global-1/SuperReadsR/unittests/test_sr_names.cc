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


#include <limits>
#include <vector>
#include <gtest/gtest.h>
#include <src/sr_names.hpp>
#include <charb.hpp>

TEST(SR_name, encode_decode) {
  long sr_name_len = random() % 10;
  charb sr_original;
  charb one_entry;
  sprintf(sr_original, "%ld%c", random(), random() % 2 ? 'R' : 'F');
  for(long i = 1; i < sr_name_len; ++i) {
    sprintf(one_entry, "_%ld_%ld%c", random() % 100, random(), random() %2 ? 'R' : 'F');
    strcat(sr_original, one_entry);
  }

  sr_name sr(sr_original);
  charb res;
  sr.to_str(res);
  EXPECT_STREQ(sr_original, res);
}
