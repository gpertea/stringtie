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


#include "fibonacci_coding.hpp"

static const uint64_t _fibs[] = { 
                    1ul,                   2ul,                   3ul,                   5ul,
                    8ul,                  13ul,                  21ul,                  34ul,
                   55ul,                  89ul,                 144ul,                 233ul,
                  377ul,                 610ul,                 987ul,                1597ul,
                 2584ul,                4181ul,                6765ul,               10946ul,
                17711ul,               28657ul,               46368ul,               75025ul,
               121393ul,              196418ul,              317811ul,              514229ul,
               832040ul,             1346269ul,             2178309ul,             3524578ul,
              5702887ul,             9227465ul,            14930352ul,            24157817ul,
             39088169ul,            63245986ul,           102334155ul,           165580141ul,
            267914296ul,           433494437ul,           701408733ul,          1134903170ul,
           1836311903ul,          2971215073ul,          4807526976ul,          7778742049ul,
          12586269025ul,         20365011074ul,         32951280099ul,         53316291173ul,
          86267571272ul,        139583862445ul,        225851433717ul,        365435296162ul,
         591286729879ul,        956722026041ul,       1548008755920ul,       2504730781961ul,
        4052739537881ul,       6557470319842ul,      10610209857723ul,      17167680177565ul,
       27777890035288ul,      44945570212853ul,      72723460248141ul,     117669030460994ul,
      190392490709135ul,     308061521170129ul,     498454011879264ul,     806515533049393ul,
     1304969544928657ul,    2111485077978050ul,    3416454622906707ul,    5527939700884757ul,
     8944394323791464ul,   14472334024676221ul,   23416728348467685ul,   37889062373143906ul,
    61305790721611591ul,   99194853094755497ul,  160500643816367088ul,  259695496911122585ul,
   420196140727489673ul,  679891637638612258ul, 1100087778366101931ul, 1779979416004714189ul,
  2880067194370816120ul, 4660046610375530309ul, 7540113804746346429ul, 12200160415121876738ul
};
const uint64_t *fibonacci::fibs = _fibs;
