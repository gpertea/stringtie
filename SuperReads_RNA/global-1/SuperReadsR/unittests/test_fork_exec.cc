#include <sys/types.h>
#include <sys/wait.h>
#include <gtest/gtest.h>

#include <fork_exec.hpp>

namespace {
TEST(ForkExec,ExecVP) {
  pid_t res;

  {
    char* cmd[2] = { (char*)"/surely this does not exists in root", (char*)0 };
    res = fork_execvp("bad", cmd);
    EXPECT_EQ(-1, res);
  }

  {
    char* cmd[2] = { (char*)"true", 0 };
    res = fork_execvp("true", cmd);
    EXPECT_LT(0, res);
    int status;
    pid_t wpid = waitpid(res, &status, 0);
    EXPECT_EQ(res, wpid);
    ASSERT_EQ(true, WIFEXITED(status));
    EXPECT_EQ(EXIT_SUCCESS, WEXITSTATUS(status));
  }
}
TEST(ForkExec,ExecLP) {
  pid_t res;

  res = fork_execlp("bad", "/this one does not exists either", (char*)0);
  EXPECT_EQ(-1, res);

  res = fork_execlp("false", "false", (char*)0);
  EXPECT_LT(0, res);
  int status;
  pid_t wpid = waitpid(res, &status, 0);
  EXPECT_EQ(res, wpid);
  ASSERT_EQ(true, WIFEXITED(status));
  EXPECT_EQ(EXIT_FAILURE, WEXITSTATUS(status));
}
TEST(ForkExec,Exec) {
  pid_t res;

  res = fork_exec("/surely_you_are_joking_Mr._Feynam is that right");
  EXPECT_EQ(-1, res);

  res = fork_exec("true ignored argument");
  EXPECT_LT(0, res);
  int status;
  pid_t wpid = waitpid(res, &status, 0);
  EXPECT_EQ(res, wpid);
  ASSERT_EQ(true, WIFEXITED(status));
  EXPECT_EQ(EXIT_SUCCESS, WEXITSTATUS(status));
}
}
