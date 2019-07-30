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
#include <iostream>
#include <fstream>
#include <gtest/gtest.h>
#include <unittests/test_main_cmdline.hpp>

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  cmdline_parse args(argc, argv);
  
  unsigned int seed;
  if(args.seed_given) {
    seed = args.seed_arg;
  } else {
    std::ifstream urandom("/dev/urandom");
    urandom.read((char*)&seed, sizeof(seed));
    if(!urandom.good()) {
      std::cerr << "Failed to read random seed" << std::endl;
      return 1;
    }
  }

  std::cout << "Using random seed " << seed << std::endl;
  srandom(seed);
  
  return RUN_ALL_TESTS();
}
