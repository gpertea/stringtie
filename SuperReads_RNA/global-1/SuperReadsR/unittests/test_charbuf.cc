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


#include <sstream>
#include <gtest/gtest.h>
#include <charbuf.hpp>

TEST(Charstream, Output) {
  charstream s(10);
  const char *l1 = "Hello\n";
  const char *l2 = "Longer text\n";
  const char *l3 = "Sho\n";

  s << l1;
  EXPECT_EQ(strlen(l1), s.size());
  std::ostringstream os1;
  os1 << s;
  EXPECT_STREQ(l1, os1.str().c_str());
  s.rewind();
  EXPECT_EQ((size_t)0, s.size());

  s << l2;
  EXPECT_EQ(strlen(l2), s.size());
  std::ostringstream os2;
  os2 << s;
  EXPECT_STREQ(l2, os2.str().c_str());
  s.rewind();
  EXPECT_EQ((size_t)0, s.size());

  s << l3;
  EXPECT_EQ(strlen(l3), s.size());
  std::ostringstream os3;
  os3 << s;
  EXPECT_STREQ(l3, os3.str().c_str());
  s.rewind();
  EXPECT_EQ((size_t)0, s.size());
}
