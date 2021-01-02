/*
Copyright (c) 2010 Marcus Geelnard
(with minor modifications by Geo Pertea)
This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software
    in a product, an acknowledgment in the product documentation would be
    appreciated but is not required.

    2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

    3. This notice may not be removed or altered from any source
    distribution.
*/

#include "GThreads.h"

#if defined(_GTHREADS_POSIX_)
  #include <unistd.h>
//  #include <map>
#elif defined(_GTHREADS_WIN32_)
  #include <process.h>
#endif
#include <string.h>

//------------------------------------------------------------------------------
// condition_variable
//------------------------------------------------------------------------------
// NOTE 1: The Win32 implementation of the condition_variable class is based on
// the corresponding implementation in GLFW, which in turn is based on a
// description by Douglas C. Schmidt and Irfan Pyarali:
// http://www.cs.wustl.edu/~schmidt/win32-cv-1.html
//
// NOTE 2: Windows Vista actually has native support for condition variables
// (InitializeConditionVariable, WakeConditionVariable, etc), but we want to
// be portable with pre-Vista Windows versions, so TinyThread++ does not use
// Vista condition variables.
//------------------------------------------------------------------------------

#if defined(_GTHREADS_WIN32_)
  #define _CONDITION_EVENT_ONE 0
  #define _CONDITION_EVENT_ALL 1
#endif


int GThread::tcounter=0;
int GThread::num_created=0;

#if defined(_GTHREADS_WIN32_)
GConditionVar::GConditionVar() : mWaitersCount(0)
{
  mEvents[_CONDITION_EVENT_ONE] = CreateEvent(NULL, FALSE, FALSE, NULL);
  mEvents[_CONDITION_EVENT_ALL] = CreateEvent(NULL, TRUE, FALSE, NULL);
  InitializeCriticalSection(&mWaitersCountLock);
}
#endif

#if defined(_GTHREADS_WIN32_)
GConditionVar::~GConditionVar()
{
  CloseHandle(mEvents[_CONDITION_EVENT_ONE]);
  CloseHandle(mEvents[_CONDITION_EVENT_ALL]);
  DeleteCriticalSection(&mWaitersCountLock);
}
#endif

#if defined(_GTHREADS_WIN32_)
void GConditionVar::_wait()
{
  // Wait for either event to become signaled due to notify_one() or
  // notify_all() being called
  int result = WaitForMultipleObjects(2, mEvents, FALSE, INFINITE);

  // Check if we are the last waiter
  EnterCriticalSection(&mWaitersCountLock);
  -- mWaitersCount;
  bool lastWaiter = (result == (WAIT_OBJECT_0 + _CONDITION_EVENT_ALL)) &&
                    (mWaitersCount == 0);
  LeaveCriticalSection(&mWaitersCountLock);

  // If we are the last waiter to be notified to stop waiting, reset the event
  if(lastWaiter)
    ResetEvent(mEvents[_CONDITION_EVENT_ALL]);
}
#endif

#if defined(_GTHREADS_WIN32_)
void GConditionVar::notify_one()
{
  // Are there any waiters?
  EnterCriticalSection(&mWaitersCountLock);
  bool haveWaiters = (mWaitersCount > 0);
  LeaveCriticalSection(&mWaitersCountLock);

  // If we have any waiting threads, send them a signal
  if(haveWaiters)
    SetEvent(mEvents[_CONDITION_EVENT_ONE]);
}
#endif

#if defined(_GTHREADS_WIN32_)
void GConditionVar::notify_all()
{
  // Are there any waiters?
  EnterCriticalSection(&mWaitersCountLock);
  bool haveWaiters = (mWaitersCount > 0);
  LeaveCriticalSection(&mWaitersCountLock);

  // If we have any waiting threads, send them a signal
  if(haveWaiters)
    SetEvent(mEvents[_CONDITION_EVENT_ALL]);
}
#endif


//------------------------------------------------------------------------------
// POSIX pthread_t to unique thread::id mapping logic.
// Note: Here we use a global thread safe std::map to convert instances of
// pthread_t to small thread identifier numbers (unique within one process).
// This method should be portable across different POSIX implementations.
//------------------------------------------------------------------------------
/*

#if defined(_GTHREADS_POSIX_)
static thread::id _pthread_t_to_ID(const pthread_t &aHandle)
{
  static mutex idMapLock;
  static std::map<pthread_t, unsigned long int> idMap;
  static unsigned long int idCount(1);

  lock_guard<mutex> guard(idMapLock);
  if(idMap.find(aHandle) == idMap.end())
    idMap[aHandle] = idCount ++;
  return thread::id(idMap[aHandle]);
}
#endif // _GTHREADS_POSIX_
*/

void gthreads_errExit(int err, const char* msg) {
	  if (msg!=NULL)
	    fprintf(stderr, "GThreads Error: %s (%s)\n", msg, strerror(err));
	  else
	    fprintf(stderr, "GThreads Error: %s\n", strerror(err));
	  exit(EXIT_FAILURE);
}



void GThread::update_counter(int inc, GThread* t_update) {
  static GMutex counterLock;
  GLockGuard<GMutex> guard(counterLock);
  if (inc==1) { //joinable thread creation
	  GThread::num_created++;
	  t_update->mId = GThread::num_created;
     }
  GThread::tcounter+=inc;
  if (t_update!=NULL && inc<0)
	  t_update->mId=0; // thread terminated

 }


//------------------------------------------------------------------------------
// thread
//------------------------------------------------------------------------------

/// Information to pass to the new thread (what to run).
struct _thread_start_info {
  /*
  void * mArg;               ///< Function argument for the thread function.
  GThread * mThread;          ///< Pointer to the thread object.
  */
  GThreadData threadData;
  //void (*mFunction)(void *, GThread*);
  void (*mFunction)(void *); ///< Pointer to the function to be executed.
  void (*gtFunction)(GThreadData&); //custom variant, passing GThreadData
  //handy constructors:
  _thread_start_info():threadData() {
      mFunction=NULL;
      gtFunction=NULL;
      }
  _thread_start_info(GThread* t, void (*aFunc)(void *), void* udata):
    threadData(udata, t) {
       mFunction=aFunc;
       gtFunction=NULL;
       }
  _thread_start_info(GThread* t, void (*gtFunc)(GThreadData &), void* udata):
    threadData(udata, t) {
       mFunction=NULL;
       gtFunction=gtFunc;
       }
};

// Thread wrapper function.
#if defined(_GTHREADS_WIN32_)
unsigned WINAPI GThread::wrapper_function(void * aArg)
#elif defined(_GTHREADS_POSIX_)
void * GThread::wrapper_function(void * aArg)
#endif
{
  // Get thread startup information
  _thread_start_info * ti = (_thread_start_info *) aArg;

/*
  try
  {
    // Call the actual client thread function
    ti->mFunction(ti->mArg, ti->mThread);
  }
  catch(...)
  {
    // Uncaught exceptions will terminate the application (default behavior
    // according to the C++11)
    std::terminate();
  }
*/
  //ti->mFunction(ti->mArg, ti->mThread);

  //cheap trick to pass current GThread pointer
  //when the user doesn't pass anything
  if (ti->gtFunction) {
    ti->gtFunction(ti->threadData);
  }
  else {
    if (ti->threadData.udata) {
      ti->mFunction(ti->threadData.udata);
    }
    else {
      ti->mFunction(ti->threadData.thread);
    }
  }
  // The thread is no longer executing
  GLockGuard<GMutex> guard(ti->threadData.thread->mDataMutex);
  ti->threadData.thread->mNotAThread = true;
  GThread::update_counter(-1, ti->threadData.thread);
  // The thread is responsible for freeing the startup information
  delete ti;

  return 0;
}


void GThread::initStart(void* tidata, size_t stacksize) {
 _thread_start_info * ti = (_thread_start_info *) tidata;
   /*ti->mFunction = aFunction;
  ti->mArg = aArg;
  ti->mThread = this;*/

  // The thread is now alive
  mNotAThread = false;

  // Create the thread
#if defined(_GTHREADS_WIN32_)
  mHandle = (HANDLE) _beginthreadex(0, 0, wrapper_function, (void *) ti, 0, &mWin32ThreadID);
#elif defined(_GTHREADS_POSIX_)
  if (stacksize>0) {
	  pthread_attr_t attr;
	  int r=pthread_attr_init(&attr);
	  if (r!=0) gthreads_errExit(r, "pthread_attr_init()");
	  r = pthread_attr_setstacksize(&attr, stacksize);
	  if (r!=0) gthreads_errExit(r, "pthread_attr_setstacksize()");
	  stack_size=stacksize;
	  r=pthread_create(&mHandle, &attr, wrapper_function, (void *) ti);
	  if (r!=0) {
		  gthreads_errExit(r, "pthread_create()");
		  //mHandle = 0;
	  }
	  r=pthread_attr_destroy(&attr);
	  if (r!=0) gthreads_errExit(r, "pthread_attr_destroy()");
  }
  else {
    int r=pthread_create(&mHandle, NULL, wrapper_function, (void *) ti);
    if (r!= 0)
    	gthreads_errExit(r, "pthread_create()");
      //mHandle = 0;
  }

#endif
  // Did we fail to create the thread?
  if(!mHandle)
  {
    mNotAThread = true;
    delete ti;
  }
   else GThread::update_counter(1, this);
}

GThread::GThread(void (*aFunction)(void *), void * aArg, size_t stacksize): mId(0), mHandle(0), mNotAThread(true)
#if defined(_GTHREADS_WIN32_)
    , mWin32ThreadID(0)
#endif
    {
  kickStart(aFunction, aArg, stacksize);
}

void GThread::kickStart(void (*aFunction)(void *), void * aArg, size_t stacksize) {
  // Serialize access to this thread structure
  GLockGuard<GMutex> guard(mDataMutex);
  // Fill out the thread startup information (passed to the thread wrapper,
  // which will eventually free it)
  _thread_start_info * ti = new _thread_start_info(this, aFunction, aArg);
  initStart(ti, stacksize);
}

//custom alternate constructor (non-C++11 compatible), passing GThreadData back to the
//user function in order to easily retrieve current GThread object
//(better alternative to this_thread)
GThread::GThread(void (*gFunction)(GThreadData& thread_data), void * aArg, size_t stacksize) {
	kickStart(gFunction, aArg, stacksize);
}

void GThread::kickStart(void (*gFunction)(GThreadData& thread_data), void * aArg, size_t stacksize) {
  // Serialize access to this thread structure
  GLockGuard<GMutex> guard(mDataMutex);

  // Fill out the thread startup information (passed to the thread wrapper,
  // which will eventually free it)
  _thread_start_info * ti = new _thread_start_info(this, gFunction, aArg);
  initStart(ti, stacksize);
}

GThread::~GThread()
{
  if(joinable()) {
    //std::terminate(); -- why??
    GThread::update_counter(-1, this);
    mDataMutex.lock();
#if defined(_TTHREAD_WIN32_)
    CloseHandle(mHandle);
#elif defined(_TTHREAD_POSIX_)
    pthread_detach(mHandle);
#endif
    mDataMutex.unlock();
    }
}

void GThread::join()
{
  if(joinable())
  {
#if defined(_GTHREADS_WIN32_)
    WaitForSingleObject(mHandle, INFINITE);
    CloseHandle(mHandle);
#elif defined(_GTHREADS_POSIX_)
    pthread_join(mHandle, NULL);
#endif
  }
}


void GThread::detach()
{
  mDataMutex.lock();
  if(!mNotAThread)
  {
#if defined(_TTHREAD_WIN32_)
    CloseHandle(mHandle);
#elif defined(_TTHREAD_POSIX_)
    pthread_detach(mHandle);
#endif
    mNotAThread = true;
  }
  mDataMutex.unlock();
}

void GThread::wait_all() {
  while (GThread::num_running()>0)
	current_thread::sleep_for(2);
}


bool GThread::joinable() const
{
  mDataMutex.lock();
  bool result = !mNotAThread;
  mDataMutex.unlock();
  return result;
}

int GThread::get_id() const
{
  if(!joinable())
    //return id();
    return 0; //FIXME: don't use this
   else
     return mId;
/*
#if defined(_GTHREADS_WIN32_)
  return id((unsigned long int) mWin32ThreadID);
#elif defined(_GTHREADS_POSIX_)
  return _pthread_t_to_ID(mHandle);
#endif
*/
}

unsigned GThread::hardware_concurrency()
{
#if defined(_GTHREADS_WIN32_)
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  return (int) si.dwNumberOfProcessors;
#elif defined(_SC_NPROCESSORS_ONLN)
  return (int) sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROC_ONLN)
  return (int) sysconf(_SC_NPROC_ONLN);
#else
  // The standard requires this function to return zero if the number of
  // hardware cores could not be determined.
  return 0;
#endif
}




//------------------------------------------------------------------------------
// current_thread
//------------------------------------------------------------------------------
/*
int current_thread::get_id() {
#if defined(_GTHREADS_WIN32_)
  return thread::id((unsigned long int) GetCurrentThreadId());
#elif defined(_GTHREADS_POSIX_)
  return _pthread_t_to_ID(pthread_self());
#endif
}
*/

void current_thread::yield() {
#if defined(_GTHREADS_WIN32_)
  Sleep(0);
#else
  sched_yield();
#endif
}
// Blocks the calling thread for a certain time (given in milliseconds)
// Example usage:
// // Sleep for 100 milliseconds:
// current_thread::sleep_for(100);
void current_thread::sleep_for(const int32_t mstime) {
#if defined(_GTHREADS_WIN32_)
  Sleep(mstime);
#else
  usleep((useconds_t)(mstime*1000));
#endif
}

