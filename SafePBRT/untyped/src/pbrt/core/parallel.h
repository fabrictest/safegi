
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef PBRT_CORE_PARALLEL_H
#define PBRT_CORE_PARALLEL_H

// core/parallel.h*
#include "pbrt.h"
#if defined(__APPLE__) && !(defined(__i386__) || defined(__amd64__))
#include <libkern/OSAtomic.h>
#endif // __APPLE__ and not x86
#ifdef WIN32
#include <windows.h>
#endif
#ifdef PBRT_HAS_PTHREADS
#include <pthread.h>
#include <semaphore.h>
#endif
#include "core/probes.h"

// Parallel Declarations
#ifdef WIN32
#if _MSC_VER >= 1300
extern "C" void _ReadWriteBarrier();
#pragma intrinsic(_ReadWriteBarrier)
#else
#define _ReadWriteBarrier()
#endif
typedef volatile LONG AtomicInt32;
#ifdef PBRT_HAS_64_BIT_ATOMICS
typedef volatile LONGLONG AtomicInt64;
#endif // 64-bit
#endif // WIN32
#ifndef WIN32
typedef volatile int32_t AtomicInt32;
#ifdef PBRT_HAS_64_BIT_ATOMICS
typedef volatile int64_t AtomicInt64;
#endif
#endif // !WIN32
inline int32_t AtomicAdd(AtomicInt32 *v, int32_t delta) {
    PBRT_ATOMIC_MEMORY_OP();
#ifdef WIN32
    // Do atomic add with MSVC inline assembly
    int32_t result;
    _ReadWriteBarrier();
    __asm {
        __asm mov edx, v
        __asm mov eax, delta
        __asm lock xadd [edx], eax
        __asm mov result, eax
    }
    _ReadWriteBarrier();
    return result + delta;
#elif defined(__APPLE__) && !(defined(__i386__) || defined(__amd64__))
    return OSAtomicAdd32Barrier(delta, v);
#else
    // Do atomic add with gcc x86 inline assembly
    int32_t origValue;
    __asm__ __volatile__("lock\n"
                         "xaddl %0,%1"
                         : "=r"(origValue), "=m"(*v) : "0"(delta)
                         : "memory");
    return origValue + delta;
#endif
}


inline int32_t AtomicCompareAndSwap(AtomicInt32 *v, int32_t newValue,
    int32_t oldValue);
inline int32_t AtomicCompareAndSwap(AtomicInt32 *v, int32_t newValue, int32_t oldValue) {
    PBRT_ATOMIC_MEMORY_OP();
#if defined(WIN32)
    return InterlockedCompareExchange(v, newValue, oldValue);
#elif defined(__APPLE__) && !(defined(__i386__) || defined(__amd64__))
    return OSAtomicCompareAndSwap32Barrier(oldValue, newValue, v);
#else
    int32_t result;
    __asm__ __volatile__("lock\ncmpxchgl %2,%1"
                          : "=a"(result), "=m"(*v)
                          : "q"(newValue), "0"(oldValue)
                          : "memory");
    return result;
#endif
}


template <typename T>
inline T *AtomicCompareAndSwapPointer(T **v, T *newValue, T *oldValue) {
    PBRT_ATOMIC_MEMORY_OP();
#if defined(WIN32)
    return InterlockedCompareExchange(v, newValue, oldValue);
#elif defined(__APPLE__) && !(defined(__i386__) || defined(__amd64__))
  #ifdef PBRT_HAS_64_BIT_ATOMICS
    return OSAtomicCompareAndSwap64Barrier(oldValue, newValue, v);
  #else
    return OSAtomicCompareAndSwap32Barrier(oldValue, newValue, v);
  #endif
#else
    T *result;
    __asm__ __volatile__("lock\ncmpxchg"
#ifdef PBRT_HAS_64_BIT_ATOMICS
                                       "q"
#else
                                       "l"
#endif // 64 bit atomics
                                        " %2,%1"
                          : "=a"(result), "=m"(*v)
                          : "q"(newValue), "0"(oldValue)
                          : "memory");
    return result;
#endif
}


#ifdef PBRT_HAS_64_BIT_ATOMICS
inline int64_t AtomicAdd(AtomicInt64 *v, int64_t delta) {
    PBRT_ATOMIC_MEMORY_OP();
#ifdef WIN32
    return InterlockedAdd64(v, delta);
#elif defined(__APPLE__) && !(defined(__i386__) || defined(__amd64__))
    return OSAtomicAdd64Barrier(delta, v);
#else
    int64_t result;
    __asm__ __volatile__("lock\nxaddq %0,%1"
                          : "=r"(result), "=m"(*v)
                          : "0"(delta)
                          : "memory");
   return result + delta;
#endif
}



inline int64_t AtomicCompareAndSwap(AtomicInt64 *v, int64_t newValue, int64_t oldValue) {
    PBRT_ATOMIC_MEMORY_OP();
#ifdef WIN32
    return InterlockedCompareExchange64(v, newValue, oldValue);
#elif defined(__APPLE__) && !(defined(__i386__) || defined(__amd64__))
    return OSAtomicCompareAndSwap64Barrier(oldValue, newValue, v);
#else
    int64_t result;
    __asm__ __volatile__("lock\ncmpxchgq %2,%1"
                          : "=a"(result), "=m"(*v)
                          : "q"(newValue), "0"(oldValue)
                          : "memory");
    return result;
#endif
}


#endif // PBRT_HAS_64_BIT_ATOMICS
inline float AtomicAdd(volatile float *val, float delta) {
    PBRT_ATOMIC_MEMORY_OP();
    union bits { float f; int32_t i; };
    bits oldVal, newVal;
    do {
        // On IA32/x64, adding a PAUSE instruction in compare/exchange loops
        // is recommended to improve performance.  (And it does!)
#if (defined(__i386__) || defined(__amd64__))
        __asm__ __volatile__ ("pause\n");
#endif
        oldVal.f = *val;
        newVal.f = oldVal.f + delta;
    } while (AtomicCompareAndSwap(((AtomicInt32 *)val),
                                  newVal.i, oldVal.i) != oldVal.i);
    return newVal.f;
}


struct MutexLock;
class Mutex {
public:
    static Mutex *Create();
    static void Destroy(Mutex *m);
private:
    // Mutex Private Methods
    Mutex();
    ~Mutex();
    friend struct MutexLock;
    Mutex(Mutex &);
    Mutex &operator=(const Mutex &);

    // System-dependent mutex implementation
#ifdef PBRT_HAS_PTHREADS
    pthread_mutex_t mutex;
#elif defined(WIN32)
    CRITICAL_SECTION criticalSection;
#endif
};


struct MutexLock {
    MutexLock(Mutex &m);
    ~MutexLock();
private:
    Mutex &mutex;
    MutexLock(const MutexLock &);
    MutexLock &operator=(const MutexLock &);
};


class RWMutex {
public:
    static RWMutex *Create();
    static void Destroy(RWMutex *m);
private:
    // RWMutex Private Methods
    RWMutex();
    ~RWMutex();
    friend struct RWMutexLock;
    RWMutex(RWMutex &);
    RWMutex &operator=(const RWMutex &);

    // System-dependent rw mutex implementation
#ifdef PBRT_HAS_PTHREADS
    pthread_rwlock_t mutex;
#elif defined(WIN32)
    void AcquireRead();
    void ReleaseRead();
    void AcquireWrite();
    void ReleaseWrite();
    
    LONG numWritersWaiting;
    LONG numReadersWaiting;
    
    // HIWORD is writer active flag;
    // LOWORD is readers active count;
    DWORD activeWriterReaders;
    
    HANDLE hReadyToRead;
    HANDLE hReadyToWrite;
    CRITICAL_SECTION cs;
#endif
};


enum RWMutexLockType { READ, WRITE };
struct RWMutexLock {
    RWMutexLock(RWMutex &m, RWMutexLockType t);
    ~RWMutexLock();
    void UpgradeToWrite();
    void DowngradeToRead();
private:
    RWMutexLockType type;
    RWMutex &mutex;
    RWMutexLock(const RWMutexLock &);
    RWMutexLock &operator=(const RWMutexLock &);
};


class Semaphore {
public:
    // Semaphore Public Methods
    Semaphore();
    ~Semaphore();
    void Post(int count = 1);
    void Wait();
    bool TryWait();
private:
    // Semaphore Private Data
#ifdef PBRT_HAS_PTHREADS
    sem_t *sem;
    static int count;
#endif // PBRT_HAS_PTHREADS
#ifdef WIN32
    HANDLE handle;
#endif // WIN32
};


class ConditionVariable {
public:
    // ConditionVariable Public Methods
    ConditionVariable();
    ~ConditionVariable();
    void Lock();
    void Unlock();
    void Wait();
    void Signal();
private:
    // ConditionVariable Private Data
#ifdef PBRT_HAS_PTHREADS
    pthread_mutex_t mutex;
    pthread_cond_t cond;
#endif // PBRT_HAS_PTHREADS
#ifdef WIN32
    // Count of the number of waiters.
    uint32_t waitersCount;
    // Serialize access to <waitersCount>.
    CRITICAL_SECTION waitersCountMutex, conditionMutex;
    // Signal and broadcast event HANDLEs.
    enum { SIGNAL = 0, BROADCAST=1, NUM_EVENTS=2 };
    HANDLE events[NUM_EVENTS];
#endif // WIN32
};


void TasksInit();
void TasksCleanup();
class Task {
public:
    virtual ~Task();
    virtual void Run() = 0;
};


void EnqueueTasks(const vector<Task *> &tasks);
void WaitForAllTasks();
int NumSystemCores();

#endif // PBRT_CORE_PARALLEL_H
