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


#include <errno.h>
#include <gtest/gtest.h>
#include <charb.hpp>
#include <ext/stdio_filebuf.h>

namespace {
TEST(CharbBasic, Init) {
  charb empty;
  EXPECT_EQ((size_t)1, empty.capacity());
  EXPECT_EQ((size_t)0, empty.len());

  charb init_len(20);
  EXPECT_EQ((size_t)21, init_len.capacity());
  EXPECT_EQ((size_t)0, init_len.len());
  EXPECT_EQ('\0', init_len[0]);
  EXPECT_EQ('\0', *init_len);
  init_len.reserve(10);
  EXPECT_EQ((size_t)21, init_len.capacity());
  init_len.reserve(30);
  EXPECT_EQ((size_t)42, init_len.capacity());

  charb copy_const(init_len);
  EXPECT_EQ(init_len.len(), copy_const.len());
  EXPECT_EQ('\0', copy_const[0]);
  charb copy_op;
  copy_op = init_len;
  EXPECT_EQ(init_len.len(), copy_op.len());
  EXPECT_EQ('\0', copy_op[0]);
}

TEST(CharbBasic, Copy) {
  static const char *str = "Hello the world";
  const size_t str_len = strlen(str);
  const std::string string(str);

  charb from_str(str);
  EXPECT_EQ(str_len + 1, from_str.capacity());
  EXPECT_EQ(str_len, from_str.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_str[i]);
  EXPECT_STREQ(str, from_str);
  EXPECT_EQ('H', from_str.front());
  EXPECT_EQ('d', from_str.back());

  charb from_str_len(str, str_len);
  EXPECT_EQ(str_len + 1, from_str_len.capacity());
  EXPECT_EQ(str_len, from_str_len.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_str_len[i]);
  EXPECT_STREQ(str, from_str_len);
  EXPECT_EQ('H', from_str_len.front());
  EXPECT_EQ('d', from_str_len.back());

  charb from_string(string);
  EXPECT_EQ(str_len + 1, from_string.capacity());
  EXPECT_EQ(str_len, from_string.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_string[i]);
  EXPECT_STREQ(str, from_string);
  EXPECT_EQ('H', from_string.front());
  EXPECT_EQ('d', from_string.back());

  charb from_copy;
  from_copy = from_string;
  EXPECT_EQ(str_len + 1, from_copy.capacity());
  EXPECT_EQ(str_len, from_copy.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_copy[i]);
  EXPECT_STREQ(str, from_copy);
  EXPECT_EQ('H', from_copy.front());
  EXPECT_EQ('d', from_copy.back());

  charb from_copy_str;
  from_copy_str = str;
  EXPECT_EQ(str_len + 1, from_copy_str.capacity());
  EXPECT_EQ(str_len, from_copy_str.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_copy_str[i]);
  EXPECT_STREQ(str, from_copy_str);
  EXPECT_EQ('H', from_copy_str.front());
  EXPECT_EQ('d', from_copy_str.back());


  charb x("Hello");
  x[3] = '_';
  EXPECT_STREQ("Hel_o", x);
}

TEST(CharbBasic, Cast) {
  static const char *str = "This is a char buffer";
  const size_t str_len = strlen(str);
  charb b(str);

  char *s = b;
  EXPECT_EQ((void*)s, (void*)&b[0]);
  EXPECT_EQ(str_len, strlen(b));
  EXPECT_STREQ(str, b);
}

TEST(CharbBasic, Strerror_r) {
  charb b;
  errno = EINVAL;
#if (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && ! _GNU_SOURCE
  int res = strerror_r(errno, b);
  EXPECT_EQ(0, res);
  EXPECT_LT(0, strlen(b));
#else
  char *res = strerror_r(errno, b);
  EXPECT_EQ(res, (char*)b);
  EXPECT_LT((size_t)0, strlen(b));
#endif
}

class CharbStd : public ::testing::Test {
public:
  virtual void SetUp() {
    s1 = "Hello you";
    s2 = "How are you";
  }
  const char *s1, *s2;
};

TEST_F(CharbStd, Strcat) {
  charb b;
  char *res;

  res = strcat(b, s1);
  EXPECT_STREQ(s1, b);
  EXPECT_STREQ(s1, res);
  EXPECT_EQ(strlen(s1), b.len());

  std::string sres(s1);
  sres += s2;
  res = strcat(b, s2);
  EXPECT_STREQ(sres.c_str(), b);
  EXPECT_EQ(sres.size(), b.len());

  charb b1(s1);
  sres += s1;
  strcat(b, b1);
  EXPECT_EQ(sres.size(),  b.len());
  EXPECT_STREQ(sres.c_str(), b);
}

TEST_F(CharbStd, Strcpy) {
  charb b;
  char *res;

  res = strcpy(b, s1);
  EXPECT_STREQ(s1, b);
  EXPECT_STREQ(s1, res);
  EXPECT_EQ(strlen(s1), b.len());
}

TEST_F(CharbStd, Clear) {
  charb b(s1);

  EXPECT_STREQ(s1, b);
  b.clear();
  EXPECT_STREQ("", b);
  EXPECT_EQ((size_t)0, strlen(b));
  EXPECT_EQ((size_t)0, b.len());
}

class IOTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    tf = tmpfile();
    if(!tf)
      throw std::runtime_error("Can't create tmp file");
    l1 = "Hello there\n";
    l2 = "A very long and meaningless sentence-line.";
    fputs(l1, tf);
    fputs(l2, tf);
    rewind(tf);
  }

  virtual void TearDown() {
    fclose(tf);
  }

  FILE *tf;
  const char *l1;
  const char *l2;
};

TEST_F(IOTest, fgets) {
  charb b(4);
  char *res = fgets(b, 3, tf); // 3 ignored!
  EXPECT_STREQ(l1, res);
  EXPECT_STREQ(l1, b);
  EXPECT_EQ((size_t)20, b.capacity());
  EXPECT_EQ(strlen(l1), b.len());

  res = fgets(b, tf);
  EXPECT_STREQ(l2, res);
  EXPECT_STREQ(l2, b);
  EXPECT_EQ((size_t)80, b.capacity());
  EXPECT_EQ(strlen(l2), b.len());

  res = fgets(b, tf);
  EXPECT_EQ((char *)0, res);
  EXPECT_STREQ(l2, b);
}

TEST_F(IOTest, fgets_empty) {
  charb b; // start empty
  char *res = fgets(b, tf);
  EXPECT_STREQ(l1, res);
  EXPECT_STREQ(l1, b);
}

TEST_F(IOTest, getline) {
  charb b(5);
  ssize_t res;

  res = getline(b, tf);
  EXPECT_EQ(strlen(l1), (size_t)res);
  EXPECT_STREQ(l1, b);
  EXPECT_LE(strlen(l1), b.len());

  res = getline(b, tf);
  EXPECT_EQ(strlen(l2), (size_t)res);
  EXPECT_STREQ(l2, b);
  EXPECT_LE(strlen(l2), b.len());

  res = getline(b, tf);
  EXPECT_EQ(-1, res);
}

TEST_F(IOTest, getline_append) {
  charb                          b(5);
  __gnu_cxx::stdio_filebuf<char> ib(tf, std::ios::in);
  std::istream                   is(&ib);

  getline_append(is, b);
  EXPECT_TRUE(is.good());
  std::string res(l1, strlen(l1) - 1); // Everything but "\n"
  EXPECT_STREQ(res.c_str(), b);

  getline_append(is, b);
  EXPECT_EQ(strlen(l1) + strlen(l2) - 1, b.len());
  EXPECT_EQ(strlen(l1) + strlen(l2) - 1, strlen(b));

}

TEST_F(IOTest, fgets_append) {
  charb b(5);
  char *res;

  res = fgets_append(b, tf);
  EXPECT_STREQ(l1, res);
  EXPECT_STREQ(l1, b);

  res = fgets_append(b, tf);
  EXPECT_STREQ(l2, res);
  std::string str_res(l1);
  str_res += l2;
  EXPECT_STREQ(str_res.c_str(), b);
  EXPECT_EQ(strlen(l1) + strlen(l2), b.len());

  EXPECT_STREQ((char*)0, fgets_append(b, tf));
}

TEST_F(IOTest, getline_cpp) {
  __gnu_cxx::stdio_filebuf<char> ib(tf, std::ios::in);
  std::istream                   is(&ib);
  charb                          b;
  std::string                    str;
  std::string::iterator          it;

  EXPECT_TRUE((bool)getline(is, b));
  str.assign(l1);  // Copy of l1 without newline
  if(*(it = str.end() - 1) == '\n')
    str.erase(it);
  EXPECT_STREQ(str.c_str(), b);

  EXPECT_TRUE((bool)getline(is, b));
  str.assign(l2);
  if(*(it = str.end() - 1) == '\n')
    str.erase(it);
  EXPECT_STREQ(str.c_str(), b);

  EXPECT_FALSE((bool)getline(is, b));
}

TEST(CharbBasic, sprintf) {
  charb b(10);
  const char *fmt = "Hello %d times";
  const char *str_res = "Hello 1000 times";
  int res = sprintf(b, fmt, 1000);
  EXPECT_EQ(strlen(str_res), (size_t)res);
  EXPECT_STREQ(str_res, b);
  EXPECT_EQ(strlen(str_res), b.len());

  res = sprintf_append(b, fmt, 1000);
  EXPECT_EQ(strlen(str_res), (size_t)res);
  EXPECT_EQ(2 * strlen(str_res), b.len());
  EXPECT_EQ(2 * strlen(str_res), strlen(b));
  EXPECT_STREQ(str_res, b + strlen(str_res));
}

TEST(CharbBasic, chomp) {
  const char * no_space_str = "Hello";
  charb no_space(no_space_str);
  no_space.chomp();
  EXPECT_STREQ(no_space_str, no_space);

  const char * empty_str = "";
  charb empty(empty_str);
  empty.chomp();
  EXPECT_STREQ(empty_str, empty);

  const char * some_space_str = "Hello \n";
  charb some_space(some_space_str);
  some_space.chomp();
  EXPECT_STREQ(no_space_str, some_space);
}

TEST(CharbBasic, truncate) {
  const char* text = "Hello the world";
  charb textb(text);

  EXPECT_EQ(strlen(text), strlen(textb));
  EXPECT_EQ(strlen(text), textb.size());

  const size_t nlen = 5;
  textb.truncate(nlen);
  EXPECT_EQ(nlen, strlen(textb));
  EXPECT_EQ(nlen, textb.size());
  EXPECT_STREQ("Hello", (char*)textb);

  const size_t nlen2 = 3 * strlen(text);
  textb.truncate(nlen2);
  EXPECT_EQ(nlen2, textb.size());
  EXPECT_LE(nlen2, textb.capacity());
  EXPECT_NE(strlen(textb), textb.size());
}
}
