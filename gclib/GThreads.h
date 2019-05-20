#ifndef _GTHREADS_
#define _GTHREADS_

/*
GThread - multi-platform threads utility class
This is heavily based on the source code of TinyThread++ 1.0 package
by Marcus Geelnard (with only minor modifications),
so all merits for the code go to Marcus, except the bugs which
were probably introduced by me.

Original Copyright notice below
*/

/*
Copyright (c) 2010 Marcus Geelnard

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


/// @file
/// @mainpage TinyThread++ API Reference
///
/// @section intro_sec Introduction
/// TinyThread++ is a minimal, portable implementation of basic threading
/// classes for C++.
///
/// They closely mimic the functionality and naming of the C++11 standard, and
/// should be easily replaceable with the corresponding std:: variants.
///
/// @section port_sec Portability
/// The Win32 variant uses the native Win32 API for implementing the thread
/// classes, while for other systems, the POSIX threads API (pthread) is used.
///
/// @section class_sec Classes
/// In order to mimic the threading API of the C++11 standard, subsets of
/// several classes are provided. The fundamental classes are:
/// @li GThread
/// @li GMutex
/// @li GRecursiveMutex
/// @li GConditionVariable
/// @li GLockGuard
/// @li GFastMutex
///
/// @section misc_sec Miscellaneous
/// The following special keywords are available: #thread_local.
///
/// For more detailed information (including additional classes), browse the
/// different sections of this documentation. A good place to start is:
/// tinythread.h.

// Which platform are we on?
#if !defined(_GTHREADS_PLATFORM_DEFINED_)
  #if defined(_WIN32) || defined(__WIN32__) || defined(__WINDOWS__)
    #define _GTHREADS_WIN32_
  #else
    #define _GTHREADS_POSIX_
  #endif
  #define _GTHREADS_PLATFORM_DEFINED_
#endif

// Check if we can support the assembly language level implementation (otherwise
// revert to the system API)
#if (defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))) || \
    (defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))) || \
    (defined(__GNUC__) && (defined(__ppc__)))
  #define _GFASTMUTEX_ASM_
#else
  #define _FAST_MUTEX_SYS_
#endif

// Platform specific includes
#if defined(_GTHREADS_WIN32_)
 #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
    #define __UNDEF_LEAN_AND_MEAN
  #endif
  #include <windows.h>
  #ifdef __UNDEF_LEAN_AND_MEAN
    #undef WIN32_LEAN_AND_MEAN
    #undef __UNDEF_LEAN_AND_MEAN
  #endif
#else
  #include <pthread.h>
  #include <signal.h>
  #include <sched.h>
  #include <unistd.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

/// TinyThread++ version (major number).
#define TINYTHREAD_VERSION_MAJOR 1
/// TinyThread++ version (minor number).
#define TINYTHREAD_VERSION_MINOR 1
/// TinyThread++ version (full version).
#define TINYTHREAD_VERSION (TINYTHREAD_VERSION_MAJOR * 100 + TINYTHREAD_VERSION_MINOR)

// Do we have a fully featured C++11 compiler?
#if (__cplusplus > 199711L) || (defined(__STDCXX_VERSION__) && (__STDCXX_VERSION__ >= 201001L))
  #define _GTHREADS_CPP11_
#endif

// ...at least partial C++11?
#if defined(_GTHREADS_CPP11_) || defined(__GXX_EXPERIMENTAL_CXX0X__) || defined(__GXX_EXPERIMENTAL_CPP0X__)
  #define _GTHREADS_CPP11_PARTIAL_
#endif

// Macro for disabling assignments of objects.
#ifdef _GTHREADS_CPP11_PARTIAL_
  #define _GTHREADS_DISABLE_ASSIGNMENT(name) \
      name(const name&) = delete; \
      name& operator=(const name&) = delete;
#else
  #define _GTHREADS_DISABLE_ASSIGNMENT(name) \
      name(const name&); \
      name& operator=(const name&);
#endif

/// @def thread_local
/// Thread local storage keyword.
/// A variable that is declared with the \c thread_local keyword makes the
/// value of the variable local to each thread (known as thread-local storage,
/// or TLS). Example usage:
/// @code
/// // This variable is local to each thread.
/// thread_local int variable;
/// @endcode
/// @note The \c thread_local keyword is a macro that maps to the corresponding
/// compiler directive (e.g. \c __declspec(thread)). While the C++11 standard
/// allows for non-trivial types (e.g. classes with constructors and
/// destructors) to be declared with the \c thread_local keyword, most pre-C++11
/// compilers only allow for trivial types (e.g. \c int). So, to guarantee
/// portable code, only use trivial types for thread local storage.
/// @note This directive is currently not supported on Mac OS X (it will give
/// a compiler error), since compile-time TLS is not supported in the Mac OS X
/// executable format. Also, some older versions of MinGW (before GCC 4.x) do
/// not support this directive.
/// @hideinitializer

#if !defined(_GTHREADS_CPP11_) && !defined(thread_local)
 #if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__SUNPRO_CC) || defined(__IBMCPP__)
  #define thread_local __thread
 #else
  #define thread_local __declspec(thread)
 #endif
#endif

// HACK: Mac OS X and early MinGW do not support thread-local storage
#if defined(__APPLE__) || (defined(__MINGW32__) && (__GNUC__ < 4))
 #define GTHREADS_NO_TLS
#endif

void gthreads_errExit(int err, const char* msg=NULL);

#define pthreads_err(msg) \
               do { perror(msg); exit(EXIT_FAILURE); } while (0)
/// GMutex class
/// This is a mutual exclusion object for synchronizing access to shared
/// memory areas for several threads. The mutex is non-recursive (i.e. a
/// program may deadlock if the thread that owns a mutex object calls lock()
/// on that object).
/// @see recursive_mutex
class GMutex {
  public:
    /// Constructor.
    GMutex()
#if defined(_GTHREADS_WIN32_)
      : mAlreadyLocked(false)
#endif
    {
#if defined(_GTHREADS_WIN32_)
      InitializeCriticalSection(&mHandle);
#else
      pthread_mutex_init(&mHandle, NULL);
#endif
    }

    /// Destructor.
    ~GMutex()
    {
#if defined(_GTHREADS_WIN32_)
      DeleteCriticalSection(&mHandle);
#else
      pthread_mutex_destroy(&mHandle);
#endif
    }

    /// Lock the mutex.
    /// The method will block the calling thread until a lock on the mutex can
    /// be obtained. The mutex remains locked until \c unlock() is called.
    /// @see lock_guard
    inline void lock()
    {
#if defined(_GTHREADS_WIN32_)
      EnterCriticalSection(&mHandle);
      while(mAlreadyLocked) Sleep(1000); // Simulate deadlock...
      mAlreadyLocked = true;
#else
      pthread_mutex_lock(&mHandle);
#endif
    }

    /// Try to lock the mutex.
    /// The method will try to lock the mutex. If it fails, the function will
    /// return immediately (non-blocking).
    /// @return \c true if the lock was acquired, or \c false if the lock could
    /// not be acquired.
    inline bool try_lock()
    {
#if defined(_GTHREADS_WIN32_)
      bool ret = (TryEnterCriticalSection(&mHandle) ? true : false);
      if(ret && mAlreadyLocked)
      {
        LeaveCriticalSection(&mHandle);
        ret = false;
      }
      return ret;
#else
      return (pthread_mutex_trylock(&mHandle) == 0) ? true : false;
#endif
    }

    /// Unlock the mutex.
    /// If any threads are waiting for the lock on this mutex, one of them will
    /// be unblocked.
    inline void unlock()
    {
#if defined(_GTHREADS_WIN32_)
      mAlreadyLocked = false;
      LeaveCriticalSection(&mHandle);
#else
      pthread_mutex_unlock(&mHandle);
#endif
    }

    _GTHREADS_DISABLE_ASSIGNMENT(GMutex)

  private:
#if defined(_GTHREADS_WIN32_)
    CRITICAL_SECTION mHandle;
    bool mAlreadyLocked;
#else
    pthread_mutex_t mHandle;
#endif

    friend class GConditionVar;
};

/// Recursive mutex class.
/// This is a mutual exclusion object for synchronizing access to shared
/// memory areas for several threads. The mutex is recursive (i.e. a thread
/// may lock the mutex several times, as long as it unlocks the mutex the same
/// number of times).
/// @see mutex
class GMutexRecursive {
  public:
    /// Constructor.
    GMutexRecursive()
    {
#if defined(_GTHREADS_WIN32_)
      InitializeCriticalSection(&mHandle);
#else
      pthread_mutexattr_t attr;
      pthread_mutexattr_init(&attr);
      pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
      pthread_mutex_init(&mHandle, &attr);
#endif
    }

    /// Destructor.
    ~GMutexRecursive()
    {
#if defined(_GTHREADS_WIN32_)
      DeleteCriticalSection(&mHandle);
#else
      pthread_mutex_destroy(&mHandle);
#endif
    }

    /// Lock the mutex.
    /// The method will block the calling thread until a lock on the mutex can
    /// be obtained. The mutex remains locked until \c unlock() is called.
    /// @see lock_guard
    inline void lock()
    {
#if defined(_GTHREADS_WIN32_)
      EnterCriticalSection(&mHandle);
#else
      pthread_mutex_lock(&mHandle);
#endif
    }

    /// Try to lock the mutex.
    /// The method will try to lock the mutex. If it fails, the function will
    /// return immediately (non-blocking).
    /// @return \c true if the lock was acquired, or \c false if the lock could
    /// not be acquired.
    inline bool try_lock()
    {
#if defined(_GTHREADS_WIN32_)
      return TryEnterCriticalSection(&mHandle) ? true : false;
#else
      return (pthread_mutex_trylock(&mHandle) == 0) ? true : false;
#endif
    }

    /// Unlock the mutex.
    /// If any threads are waiting for the lock on this mutex, one of them will
    /// be unblocked.
    inline void unlock()
    {
#if defined(_GTHREADS_WIN32_)
      LeaveCriticalSection(&mHandle);
#else
      pthread_mutex_unlock(&mHandle);
#endif
    }

    _GTHREADS_DISABLE_ASSIGNMENT(GMutexRecursive)

  private:
#if defined(_GTHREADS_WIN32_)
    CRITICAL_SECTION mHandle;
#else
    pthread_mutex_t mHandle;
#endif

    friend class GConditionVar;
};

/// Fast mutex class.
/// This is a mutual exclusion object for synchronizing access to shared
/// memory areas for several threads. It is similar to the tthread::mutex class,
/// but instead of using system level functions, it is implemented as an atomic
/// spin lock with very low CPU overhead.
///
/// The \c fast_mutex class is NOT compatible with the \c condition_variable
/// class (however, it IS compatible with the \c lock_guard class). It should
/// also be noted that the \c fast_mutex class typically does not provide
/// as accurate thread scheduling as a the standard \c mutex class does.
///
/// Because of the limitations of the class, it should only be used in
/// situations where the mutex needs to be locked/unlocked very frequently.
///
/// @note The "fast" version of this class relies on inline assembler language,
/// which is currently only supported for 32/64-bit Intel x86/AMD64 and
/// PowerPC architectures on a limited number of compilers (GNU g++ and MS
/// Visual C++).
/// For other architectures/compilers, system functions are used instead.

class GFastMutex {
  public:
    /// Constructor.
#if defined(_GFASTMUTEX_ASM_)
    GFastMutex() : mLock(0) {}
#else
    GFastMutex()
    {
  #if defined(_GTHREADS_WIN32_)
      InitializeCriticalSection(&mHandle);
  #elif defined(_GTHREADS_POSIX_)
      pthread_mutex_init(&mHandle, NULL);
  #endif
    }
#endif

#if !defined(_GFASTMUTEX_ASM_)
    /// Destructor.
    ~GFastMutex()
    {
  #if defined(_GTHREADS_WIN32_)
      DeleteCriticalSection(&mHandle);
  #elif defined(_GTHREADS_POSIX_)
      pthread_mutex_destroy(&mHandle);
  #endif
    }
#endif

    /// Lock the mutex.
    /// The method will block the calling thread until a lock on the mutex can
    /// be obtained. The mutex remains locked until \c unlock() is called.
    /// @see lock_guard
    inline void lock()
    {
#if defined(_GFASTMUTEX_ASM_)
      bool gotLock;
      do {
        gotLock = try_lock();
        if(!gotLock)
        {
  #if defined(_GTHREADS_WIN32_)
          Sleep(0);
  #elif defined(_GTHREADS_POSIX_)
          sched_yield();
  #endif
        }
      } while(!gotLock);
#else
  #if defined(_GTHREADS_WIN32_)
      EnterCriticalSection(&mHandle);
  #elif defined(_GTHREADS_POSIX_)
      pthread_mutex_lock(&mHandle);
  #endif
#endif
    }

    /// Try to lock the mutex.
    /// The method will try to lock the mutex. If it fails, the function will
    /// return immediately (non-blocking).
    /// @return \c true if the lock was acquired, or \c false if the lock could
    /// not be acquired.
    inline bool try_lock()
    {
#if defined(_GFASTMUTEX_ASM_)
      int oldLock;
  #if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
      asm volatile (
        "movl $1,%%eax\n\t"
        "xchg %%eax,%0\n\t"
        "movl %%eax,%1\n\t"
        : "=m" (mLock), "=m" (oldLock)
        :
        : "%eax", "memory"
      );
  #elif defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
      int *ptrLock = &mLock;
      __asm {
        mov eax,1
        mov ecx,ptrLock
        xchg eax,[ecx]
        mov oldLock,eax
      }
  #elif defined(__GNUC__) && (defined(__ppc__))
      int newLock = 1;
      asm volatile (
        "\n1:\n\t"
        "lwarx  %0,0,%1\n\t"
        "cmpwi  0,%0,0\n\t"
        "bne-   2f\n\t"
        "stwcx. %2,0,%1\n\t"
        "bne-   1b\n\t"
        "isync\n"
        "2:\n\t"
        : "=&r" (oldLock)
        : "r" (&mLock), "r" (newLock)
        : "cr0", "memory"
      );
  #endif
      return (oldLock == 0);
#else
  #if defined(_GTHREADS_WIN32_)
      return TryEnterCriticalSection(&mHandle) ? true : false;
  #elif defined(_GTHREADS_POSIX_)
      return (pthread_mutex_trylock(&mHandle) == 0) ? true : false;
  #endif
#endif
    }

    /// Unlock the mutex.
    /// If any threads are waiting for the lock on this mutex, one of them will
    /// be unblocked.
    inline void unlock()
    {
#if defined(_GFASTMUTEX_ASM_)
  #if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
      asm volatile (
        "movl $0,%%eax\n\t"
        "xchg %%eax,%0\n\t"
        : "=m" (mLock)
        :
        : "%eax", "memory"
      );
  #elif defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
      int *ptrLock = &mLock;
      __asm {
        mov eax,0
        mov ecx,ptrLock
        xchg eax,[ecx]
      }
  #elif defined(__GNUC__) && (defined(__ppc__))
      asm volatile (
        "sync\n\t"  // Replace with lwsync where possible?
        : : : "memory"
      );
      mLock = 0;
  #endif
#else
  #if defined(_GTHREADS_WIN32_)
      LeaveCriticalSection(&mHandle);
  #elif defined(_GTHREADS_POSIX_)
      pthread_mutex_unlock(&mHandle);
  #endif
#endif
    }

  private:
#if defined(_GFASTMUTEX_ASM_)
    int mLock;
#else
  #if defined(_GTHREADS_WIN32_)
    CRITICAL_SECTION mHandle;
  #elif defined(_GTHREADS_POSIX_)
    pthread_mutex_t mHandle;
  #endif
#endif
};



/// Lock guard class.
/// The constructor locks the mutex, and the destructor unlocks the mutex, so
/// the mutex will automatically be unlocked when the lock guard goes out of
/// scope. Example usage:
/// @code
/// mutex m;
/// int counter;
///
/// void increment()
/// {
///   lock_guard<mutex> guard(m);
///   ++ counter;
/// }
/// @endcode

template <class T>
class GLockGuard {
  public:
    typedef T mutex_type;

    GLockGuard() : mMutex(0) {}

    /// The constructor locks the mutex.
    explicit GLockGuard(mutex_type &aMutex)
    {
      mMutex = &aMutex;
      mMutex->lock();
    }

    /// The destructor unlocks the mutex.
    ~GLockGuard()
    {
      if(mMutex)
        mMutex->unlock();
    }

  private:
    mutex_type * mMutex;
};

/// Condition variable class.
/// This is a signalling object for synchronizing the execution flow for
/// several threads. Example usage:
/// @code
/// // Shared data and associated mutex and condition variable objects
/// int count;
/// mutex m;
/// condition_variable cond;
///
/// // Wait for the counter to reach a certain number
/// void wait_counter(int targetCount)
/// {
///   lock_guard<mutex> guard(m);
///   while(count < targetCount)
///     cond.wait(m);
/// }
///
/// // Increment the counter, and notify waiting threads
/// void increment()
/// {
///   lock_guard<mutex> guard(m);
///   ++ count;
///   cond.notify_all();
/// }
/// @endcode
class GConditionVar {
  public:
    /// Constructor.
#if defined(_GTHREADS_WIN32_)
    GConditionVar();
#else
    GConditionVar()
    {
      pthread_cond_init(&mHandle, NULL);
    }
#endif

    /// Destructor.
#if defined(_GTHREADS_WIN32_)
    ~GConditionVar();
#else
    ~GConditionVar()
    {
      pthread_cond_destroy(&mHandle);
    }
#endif

    /// Wait for the condition.
    /// The function will block the calling thread until the condition variable
    /// is woken by \c notify_one(), \c notify_all() or a spurious wake up.
    /// @param[in] aMutex A mutex that will be unlocked when the wait operation
    ///   starts, an locked again as soon as the wait operation is finished.
    template <class _mutexT>
    inline void wait(_mutexT &aMutex)
    {
#if defined(_GTHREADS_WIN32_)
      // Increment number of waiters
      EnterCriticalSection(&mWaitersCountLock);
      ++ mWaitersCount;
      LeaveCriticalSection(&mWaitersCountLock);

      // Release the mutex while waiting for the condition (will decrease
      // the number of waiters when done)...
      aMutex.unlock();
      _wait();
      aMutex.lock();
#else
      pthread_cond_wait(&mHandle, &aMutex.mHandle);
#endif
    }

    /// Notify one thread that is waiting for the condition.
    /// If at least one thread is blocked waiting for this condition variable,
    /// one will be woken up.
    /// @note Only threads that started waiting prior to this call will be
    /// woken up.
#if defined(_GTHREADS_WIN32_)
    void notify_one();
#else
    inline void notify_one()
    {
      pthread_cond_signal(&mHandle);
    }
#endif

    /// Notify all threads that are waiting for the condition.
    /// All threads that are blocked waiting for this condition variable will
    /// be woken up.
    /// @note Only threads that started waiting prior to this call will be
    /// woken up.
#if defined(_GTHREADS_WIN32_)
    void notify_all();
#else
    inline void notify_all()
    {
      pthread_cond_broadcast(&mHandle);
    }
#endif
    _GTHREADS_DISABLE_ASSIGNMENT(GConditionVar)
  private:
#if defined(_GTHREADS_WIN32_)
    void _wait();
    HANDLE mEvents[2];                  ///< Signal and broadcast event HANDLEs.
    unsigned int mWaitersCount;         ///< Count of the number of waiters.
    CRITICAL_SECTION mWaitersCountLock; ///< Serialize access to mWaitersCount.
#else
    pthread_cond_t mHandle;
#endif
};


class GThread;

struct GThreadData {
  void* udata; //user data
  GThread* thread; //current GThread object
  GThreadData(void* u=NULL, GThread* t=NULL):udata(u),thread(t) {}
};

/// Thread class.
class GThread {
  public:
#if defined(_GTHREADS_WIN32_)
    typedef HANDLE native_handle_type;
#else
    typedef pthread_t native_handle_type;
#endif
 private:
    int mId;
    size_t stack_size; //available only for pthreads
    static int tcounter; //counts live, joinable GThread instances
    static int num_created; //counts all joinable GThread instances ever created by current process

    native_handle_type mHandle;   ///< Thread handle.
    mutable GMutex mDataMutex;     ///< Serializer for access to the thread private data.
    bool mNotAThread;             ///< True if this object is not a thread of execution.
#if defined(_GTHREADS_WIN32_)
    unsigned int mWin32ThreadID;  ///< Unique thread ID (filled out by _beginthreadex).
#endif
public:
    /// Default constructor.
    /// Construct a thread object without an associated thread of execution
    /// (i.e. non-joinable).
    GThread(size_t stacksize=0) : mId(0), stack_size(stacksize), mHandle(0), mNotAThread(true)
#if defined(_GTHREADS_WIN32_)
    , mWin32ThreadID(0)
#endif
    {}

    /// Thread starting constructor.
    /// Construct a @c thread object with a new thread of execution.
    /// @param[in] aFunction A function pointer to a function of type:
    ///          <tt>void fun(void * arg)</tt>
    /// @param[in] aArg Argument to the thread function.
    /// @note This constructor is not fully compatible with the standard C++
    /// thread class. It is more similar to the pthread_create() (POSIX) and
    /// CreateThread() (Windows) functions.
    //GThread(void (*aFunction)(void *, GThread*), void * aArg);
    GThread(void (*aFunction)(void *), void * aArg=NULL, size_t stacksize=0);

    GThread(void (*aFunction)(GThreadData& thread_data), void * aArg, size_t stacksize=0);

    void kickStart(void (*aFunction)(GThreadData& thread_data), void * aArg, size_t stacksize=0);
    void kickStart(void (*aFunction)(void *), void * aArg=NULL, size_t stacksize=0);

    /// Destructor.
    /// @note If the thread is joinable upon destruction, \c std::terminate()
    /// will be called, which terminates the process. It is always wise to do
    /// \c join() before deleting a thread object.
    ~GThread();

    /// Wait for the thread to finish (join execution flows).
    void join();

    /// Check if the thread is joinable.
    /// A thread object is joinable if it has an associated thread of execution.
    bool joinable() const;

    /// Detach from the thread.
    /// After calling @c detach(), the thread object is no longer assicated with
    /// a thread of execution (i.e. it is not joinable). The thread continues
    /// execution without the calling thread blocking, and when the thread
    /// ends execution, any owned resources are released.
    void detach();
    /// Return the thread ID of a thread object.
    int get_id() const; // { return mID; }
    size_t getStackSize() { return stack_size; } //only for pthreads
    /// Get the native handle for this thread.
    /// @note Under Windows, this is a \c HANDLE, and under POSIX systems, this
    /// is a \c pthread_t.
    inline native_handle_type native_handle()
    {
      return mHandle;
    }

    inline void yield() {
#if defined(_GTHREADS_WIN32_)
        Sleep(0);
#else
        sched_yield();
#endif
    }
    static int num_running() {
    	//return number of running (live) threads
    	static GFastMutex vLock;
    	GLockGuard<GFastMutex> guard(vLock);
    	int r=tcounter;
    	return r;
    }
    static size_t defaultStackSize() {
    	pthread_attr_t attr;
    	size_t stacksize;
    	pthread_attr_init(&attr);
    	pthread_attr_getstacksize(&attr, &stacksize);
    	pthread_attr_destroy(&attr);
    	return stacksize;
    }
    static int liveCount() {
      //return number of running (live) threads
      return num_running();
    }
    static void wait_all();
    /// Determine the number of threads which can possibly execute concurrently.
    /// This function is useful for determining the optimal number of threads to
    /// use for a task.
    /// @return The number of hardware thread contexts in the system.
    /// @note If this value is not defined, the function returns zero (0).
    static unsigned hardware_concurrency();

    _GTHREADS_DISABLE_ASSIGNMENT(GThread)

  private:
    void initStart(void* tidata, size_t stacksize=0);
    static void update_counter(int inc=1, GThread* t_update=NULL); //default: increments
    // This is the internal thread wrapper function.
#if defined(_GTHREADS_WIN32_)
    static unsigned WINAPI wrapper_function(void * aArg);
#else
    static void * wrapper_function(void * aArg);
#endif
};


/// The namespace "current_thread" provides methods for dealing with the
/// calling thread.
namespace current_thread {
  /// Return the thread ID of the calling thread.
  //thread::id get_id(); //this can be slow, better not use it

  /// Yield execution to another thread.
  /// Offers the operating system the opportunity to schedule another thread
  /// that is ready to run on the current processor.
  void yield();

  // Blocks the calling thread for a certain time (given in milliseconds)
  // Example usage:
  // // Sleep for 100 milliseconds:
  // current_thread::sleep_for(100);
  void sleep_for(const int32_t mstime);
}

// Define/macro cleanup
#undef _GTHREADS_DISABLE_ASSIGNMENT

#endif // _GTHREADS_
