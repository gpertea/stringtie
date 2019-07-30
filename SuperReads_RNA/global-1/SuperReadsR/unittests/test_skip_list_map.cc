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
#include <map>
#include <string>
#include <utility>
#include <skip_list_map.hpp>

std::ostream& operator<<(std::ostream& os, const std::pair<std::string, int>& rhs) {
  return os << "(\"" << rhs.first << "\", " << rhs.second << ")";
}

TEST(SkipListMap, InserIterator) {
  typedef skip_list_map<std::string, int> sls_type;
  typedef std::map<std::string, int> map_type;
  sls_type sls(5);
  map_type map;

  std::string str(100, 'A');
  for(int i = 0; i < 100; ++i) {
    int sum = 0;
    for(int j = 0; j < 100; ++j) {
      int x = random() % 26;
      str[j] = 'A' + x;
      sum += x;
    }
    sls[str] = sum;
    map[str] = sum;
  }
  
  EXPECT_EQ(map.size(), sls.size());
  EXPECT_EQ(map.empty(), sls.empty());
  EXPECT_TRUE(std::equal(map.begin(), map.end(), sls.begin()));

  const sls_type csls(sls);
  for(map_type::const_iterator it = map.begin(); it != map.end(); ++it) {
    ASSERT_FALSE(sls.find(it->first) == sls.end());
    EXPECT_EQ(*it, *sls.find(it->first));
    ASSERT_FALSE(csls.find(it->first) == sls.end());
    EXPECT_EQ(*it, *csls.find(it->first));
    ASSERT_TRUE(sls.lower_bound(it->first) != sls.end());
    EXPECT_EQ(*it, *sls.lower_bound(it->first));
    ASSERT_TRUE(csls.lower_bound(it->first) != sls.end());
    EXPECT_EQ(*it, *csls.lower_bound(it->first));
    
    map_type::iterator map_it = map.upper_bound(it->first);
    sls_type::iterator sls_it = sls.upper_bound(it->first);
    sls_type::const_iterator csls_it = csls.upper_bound(it->first);
    if(map_it == map.end()) {
      EXPECT_TRUE(sls.end() == sls_it);
      EXPECT_TRUE(csls.end() == csls_it);
    } else {
      EXPECT_EQ(*map_it, *sls_it);
      EXPECT_EQ(*map_it, *csls_it);
    }
  }
}
