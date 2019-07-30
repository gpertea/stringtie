#include <gtest/gtest.h>
#include <thread_pool.hpp>

namespace {
int square(int x) { return x * x; }

TEST(ThreadPool, Squares) {
  static const int nb_threads = 10;
  static const int nb_jobs = 100 * nb_threads;

  int res[nb_jobs];

  thread_pool<int, int> tp(nb_threads, square);
  for(int i = 0; i < nb_jobs; ++i)
    tp.submit_job(&i, &res[i]);
  tp.release_workers();

  for(int i = 0; i < nb_jobs; ++i)
    EXPECT_EQ(i * i, res[i]);

  EXPECT_EQ(nb_jobs, tp.jobs_done());
}
}
