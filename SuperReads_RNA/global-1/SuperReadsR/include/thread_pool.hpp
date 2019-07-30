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

#ifndef _THREAD_POOL_H_
#define _THREAD_POOL_H_

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>

#include <stdexcept>


//#define DEBUG 1
#undef DEBUG

/// Manage a pool of threads. The pool of threads is templatized on
/// two types: the input type `T` and the output type `U` of the
/// function that the threads in the pool execute. The function takes
/// type `T` as argument and return type `U`. The types `T` and `U`
/// must be copyable.
///
/// Here is an example, using a thread pool to compute the squares of
/// numbers. The square of the number `i` will be written in the array
/// cell `res[i]`.
///
/// ~~~~{.cc}
/// static const int nb_threads = 10;
/// static const int nb_jobs = 100 * nb_threads;
///
/// long square(int x) { return (long)x * (long)x; }
/// long res[nb_jobs];
///
/// thread_pool<int, long> tp(nb_threads, square);
/// for(int i = 0; i < nb_jobs; ++i)
///   tp.submit_job(&i, &res[i]);
/// tp.release_workers(); // Wait for computation to be done
/// ~~~~
///
/// @see thread_pool().
template <typename T, typename U>
class thread_pool {
  typedef U (*function_ptr)(T);
  struct arguments {
    T *f_arg;
    U *retVal;
  };
  struct thread_info {
    pthread_t id;
    bool      started;
  };

  int              num_threads;
  thread_info     *thread_id;
  arguments        arg; 
  volatile bool    data_available;
  volatile bool    work_finished;
  volatile int     done;
  struct timespec  t;
  function_ptr     F;

  pthread_mutex_t lock1, lock2, lock3;
  pthread_cond_t  cond_avl, cond_acq;

  static void* __startThread(void* self) {
    ((thread_pool*)self)->worker();
    return 0;
  }

  void worker() {
    U *retValue;
    T  argumentLocal;

#ifdef DEBUG
    printf("Thread %ld is alive!!! %p\n", pthread_self(), this);
#endif
    while(true){
#ifdef DEBUG
      printf("Worker waiting for data %d\n",data_available);
#endif
      pthread_mutex_lock(&lock1);
      if(work_finished)
        goto finished;
      while(!data_available) {
        pthread_cond_wait(&cond_avl,&lock1);
        if(work_finished)
          goto finished;
      }
#ifdef DEBUG
      printf("Worker copying data\n");
#endif
      argumentLocal  = *arg.f_arg;
      retValue       = arg.retVal;
      pthread_mutex_lock(&lock3);
      data_available = false;
      pthread_cond_signal(&cond_acq);
      pthread_mutex_unlock(&lock3);
      pthread_mutex_unlock(&lock1);

      /* call the function */
      *retValue = (*F)(argumentLocal);

      pthread_mutex_lock(&lock2);
      done++;
#ifdef DEBUG
      printf("Thread %d: done = %d\n",(int)pthread_self(),done);
#endif
      pthread_mutex_unlock(&lock2);
    }

  finished:
    pthread_mutex_unlock(&lock1);
  }

public:
  /// Create a thread pool with the given number of threads. The
  /// threads will execute the function F_.
  ///
  /// @param num_threads_ Number of threads to run
  ///
  /// @param F_ Function to run in each thread. The prototype should
  /// be U f_(T). I.e. the return type is T and the parameter is type
  /// U.
  thread_pool(int num_threads_, function_ptr F_) :
    num_threads(num_threads_),
    thread_id(new thread_info[num_threads]),
    arg(),
    data_available(false),
    work_finished(false),
    done(0),
    t({0, 5000}),
    F(F_)
  {
#ifdef DEBUG
    printf("Pool %p\n", this);
#endif
    pthread_mutex_init(&lock1,NULL);
    pthread_mutex_init(&lock2,NULL);
    pthread_mutex_init(&lock3,NULL);
    pthread_cond_init(&cond_avl,NULL);
    pthread_cond_init(&cond_acq,NULL);

    for(int i = 0; i < num_threads; ++i)
      thread_id[i].started = false;

    for(int i=0;i<num_threads;i++) {
      if(pthread_create(&thread_id[i].id, NULL, __startThread, this))
        throw std::runtime_error("Failed to start a thread in the pool");
      thread_id[i].started = true;
    }
  }

  /// Deleted.
  /// A thread pool is not copyable.
  thread_pool(const thread_pool& rhs) = delete;
  /// Deleted.
  /// A thread pool is not copyable.
  thread_pool& operator=(const thread_pool& rhs) = delete;

  /// Delete the thread pool. It first waits for all the jobs to be done.
  /// @see release_workers
  ~thread_pool(void){
    release_workers();

    delete [] thread_id;
    pthread_mutex_destroy(&lock1);
    pthread_mutex_destroy(&lock2);
    pthread_mutex_destroy(&lock3);
    pthread_cond_destroy(&cond_acq);
    pthread_cond_destroy(&cond_avl);
  }

  /// The number of jobs completed so far.
  ///
  /// @return Number of jobs done
  int jobs_done(void){ return done; }

  /// Notify the workers than no more jobs will be submitted. Every
  /// thread will return after it is done its current processing. The
  /// method returns when all threads are done and have returned.
  void release_workers(void){
    pthread_mutex_lock(&lock1);
    work_finished = true;
    pthread_cond_broadcast(&cond_avl);
    pthread_mutex_unlock(&lock1);

    for(int i=0;i<num_threads;i++) {
      if(thread_id[i].started) {
        pthread_join(thread_id[i].id, NULL);
        thread_id[i].started = false;
      }
    }
  }

  /// Submit a job to a thread in the pool.  If no thread is available
  /// at the time of call, the method will block until a thread
  /// finishes its current job and makes itself available.
  /// 
  /// The `*arg` will be copied internally for the thread to work
  /// on. I.e. it does not need to be kept around after the submit_job
  /// method returns. The result of the computation will eventually be
  /// written into `*ret`. Hence `ret` must point to a memory location that is
  /// valid for as long as the worker thread is running.
  ///
  /// @param arg_ Argument to thread
  /// @param ret  Pointer to location for copying result
  void submit_job(T* arg_, U* ret){
#ifdef DEBUG
    printf("Waking worker %d\n",data_available);
#endif
    pthread_mutex_lock(&lock1);
    arg.f_arg      = arg_;
    arg.retVal     = ret;
    data_available = true;
    pthread_cond_signal(&cond_avl);
    pthread_mutex_unlock(&lock1);
#ifdef DEBUG
    printf("Signaled worker %d\n",data_available);
#endif
    /*wait until thread receives the data*/
    pthread_mutex_lock(&lock3);
    while(data_available)
      pthread_cond_wait(&cond_acq,&lock3);
    pthread_mutex_unlock(&lock3);
#ifdef DEBUG
    printf("Worker received data %d\n",data_available);
#endif

  }
};

#endif /* _THREAD_POOL_H_ */
