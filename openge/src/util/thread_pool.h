/*********************************************************************
 *
 * thread_pool.cpp: Various structures and tools for introducing 
 *                  parallelism in OGE classes.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 8 April 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#ifndef OGE_THREAD_POOL_H
#define OGE_THREAD_POOL_H

#include <vector>
#include <queue>
#include <cassert>
#include <pthread.h>
#include <semaphore.h>

#include <algorithm>

#include <unistd.h>

#ifdef __APPLE__
#include <libkern/OSAtomic.h>
#endif

#ifdef __linux__
#include <sys/prctl.h>
#define ogeNameThread(X) prctl(PR_SET_NAME,X,0,0,0)
#else 
#define ogeNameThread(X)
#endif

class ThreadPool;

// The ThreadJob abstract class provides a way to provide jobs to 
// the ThreadPool thread pool. To use, create a thread pool, subclass ThreadJob,
// implement the runJob method of your subclass, and pass an instance to the pool's
// addJob method.
class ThreadJob
{
	friend class ThreadPool;
public:
    virtual ~ThreadJob();
	virtual void runJob() = 0;
protected:
};

class Spinlock
{
private:    //private copy-ctor and assignment operator ensure the lock never gets copied, which might cause issues.
    Spinlock operator=(const Spinlock & asdf);
    Spinlock(const Spinlock & asdf);
#ifdef __APPLE__
    OSSpinLock m_lock;
public:
    Spinlock()
    : m_lock(0)
    {}
    void lock() {
        OSSpinLockLock(&m_lock);
    }
    void unlock() {
        OSSpinLockUnlock(&m_lock);
    }
#else
    pthread_spinlock_t m_lock;
public:
    Spinlock() {
        pthread_spin_init(&m_lock, 0);
    }
    
    void lock() {
        pthread_spin_lock(&m_lock);
    }
    void unlock() {
        pthread_spin_unlock(&m_lock);
    }
    ~Spinlock() {
        pthread_spin_destroy(&m_lock);
    }
#endif
};

// The mutex and condition_variable classes are designed to be
// compatible with C++11. Eventually we will get rid of these
// and use the C++11 equivalents.
class mutex {
    pthread_mutex_t m;
public:
    mutex();
    ~mutex();
    void lock();
    bool try_lock();
    void unlock();
    pthread_mutex_t & native_handle() { return m; };
};

class condition_variable {
    pthread_cond_t c;
public:
    condition_variable();
    ~condition_variable();
    void notify_one();
    void notify_all();
    void wait(mutex & m);
    pthread_cond_t & native_handle() { return c; };
};


template <class T>
class SynchronizedQueue
{
public: 
    SynchronizedQueue() { }
    
    void push(const T & item) {
        lock.lock();
        q.push(item);
        lock.unlock();
    }
    
    size_t size() {
        lock.lock();
        size_t r = q.size();
        lock.unlock();
        return r;
    }
    
    bool empty() {
        lock.lock();
        bool e = q.empty();
        lock.unlock();
        return e;
    }
    
    T pop() {
        lock.lock();
        T ret = q.front();
        q.pop();
        lock.unlock();
        return ret;
    }
    
    T & back() {
        lock.lock();
        T & ret = q.back();
        lock.unlock();
        return ret;
    }
    
    T & front() {
        lock.lock();
        T & ret = q.front();
        lock.unlock();
        return ret;
    }
    
protected:
    
    Spinlock lock;
    std::queue<T> q;
};

template<class T>
class SynchronizedBlockingQueue : public SynchronizedQueue<T>{
public:
    SynchronizedBlockingQueue()
    : cv_waiting(false)
    {}
    
    void push(const T & item) {
        this->lock.lock();
        this->q.push(item);
        bool waiting = cv_waiting;
        this->lock.unlock();
        
        if(waiting) {
            wait_lock.lock();
            wait_cv.notify_one();
            wait_lock.unlock();
        }
    }

    T pop() {
        bool success = false;
        T ret;
        this->lock.lock();
        if(!this->q.empty()) {
            ret = this->q.front();
            this->q.pop();
            success = true;
            this->lock.unlock();
        } else {
            cv_waiting = true;
            wait_lock.lock();
            this->lock.unlock();
            while(true) {
                this->lock.lock();
                if(!this->q.empty()) {
                    success = true;
                    ret = this->q.front();
                    this->q.pop();
                    cv_waiting = false;
                }
                this->lock.unlock();
                if(success) {
                    break;
                }
                wait_cv.wait(wait_lock);
            }
            wait_lock.unlock();
        }
        return ret;
    }
protected:
    bool cv_waiting;
    mutex wait_lock;
    condition_variable wait_cv;
};

class SynchronizedFlag {
    Spinlock s;
    bool b;
public:
    SynchronizedFlag() {}
    SynchronizedFlag(const bool b) {
        s.lock();
        this->b = b;
        s.unlock();
    }
    bool operator=(const bool b) {
        s.lock();
        this->b = b;
        s.unlock();
        return b;
    }
    
    void set() {
        s.lock();
        this->b = true;
        s.unlock();
    }
    
    void clear() {
        s.lock();
        this->b = false;
        s.unlock();
    }
    
    bool isSet() {
        s.lock();
        bool r = b;
        s.unlock();
        return r;
    }
};

class OGEParallelismSettings
{
public:
    // return the number of cores availble in the current system
    static int availableCores();
    
    // Configure the number of threads to use for each thread pool. Set
    // 0 to use all threads on the system.
    static void setNumberThreads(int threads);
    static int getNumberThreads();
    
    // Enable or disable the use of multithreading
    static void disableMultithreading();
    static void enableMultithreading();
    static bool isMultithreadingEnabled() { return m_multithreading_enabled; }
    
protected:
    static int m_configured_threads;
    static bool m_multithreading_enabled;
};

// The shared thread pool should only be used for non-blocking tasks!

class ThreadPool
{
	friend class ThreadJob;
public:
	ThreadPool( int threads = 0);
	virtual ~ThreadPool();
	
	bool addJob(ThreadJob * job);
	int numJobs();
	static int availableCores();
	void waitForJobCompletion();
    static ThreadPool * sharedPool();
    static void closeSharedPool();
    static bool sharedPoolIsStarted() { return NULL != _sharedPool; }
	
protected:
    static ThreadPool * _sharedPool;
	static void * thread_start(void * thread_pool);
	ThreadJob * startJob();
	void stopJob(ThreadJob * job);
    
    pthread_cond_t job_queue_cond, busy_cond;
    pthread_mutex_t job_queue_mutex;
	
    std::queue<ThreadJob *> jobs;	//protected by jobs_mutex
	std::vector<pthread_t> threads;
	SynchronizedFlag threads_exit;
	int jobs_in_process;	//protected by jobs_mutex
	Spinlock jobs_mutex;
	pthread_mutex_t busy_mutex;
    int jobs_current;
};

template<typename _RandomAccessIterator, typename _Compare>
class OGESortJob : public ThreadJob
{
protected:
    _RandomAccessIterator first, last;
    _Compare comp;
    int & completion_ct;
    Spinlock & completion_spinlock;

    virtual void runJob() {
        std::sort(first, last, comp);
        completion_spinlock.lock();
        completion_ct++;
        completion_spinlock.unlock();
    }
public:
    OGESortJob(_RandomAccessIterator __first, _RandomAccessIterator __last,
               _Compare __comp, int & completion_ct, Spinlock & completion_spinlock)
    : first(__first)
    , last(__last)
    , comp(__comp)
    , completion_ct(completion_ct)
    , completion_spinlock(completion_spinlock)
    {}
};

template<typename _RandomAccessIterator, typename _Compare>
inline void
ogeSortMt(_RandomAccessIterator __first, _RandomAccessIterator __last,
	 _Compare __comp)
{
    if(!OGEParallelismSettings::isMultithreadingEnabled())
        return sort(__first, __last, __comp);
    
    size_t job_size = (__last - __first) / ThreadPool::availableCores();
    
    ThreadPool * shared_pool = ThreadPool::sharedPool();
    int completion_ct = 1;
    Spinlock completion_spinlock;
    
    //perform separate sorts
    for(int i = 0; i < ThreadPool::availableCores() - 1; i++) {
        shared_pool->addJob(new OGESortJob<_RandomAccessIterator, _Compare>(__first + i * job_size, __first + (i+1) * job_size, __comp, completion_ct, completion_spinlock));
    }

    std::sort(__first + (ThreadPool::availableCores() - 1) * job_size, __last, __comp);

    //wait for completion:
    while(true) {
        completion_spinlock.lock();
        bool complete = (ThreadPool::availableCores() == completion_ct);
        completion_spinlock.unlock();
        if(complete)
            break;
        usleep(20000);
    }

    //now merge sorted subarrays
    for(int i = 1; i < ThreadPool::availableCores(); i++)  {
        if(i == ThreadPool::availableCores() - 1)   //last
            inplace_merge(__first, __first + i * job_size, __last, __comp);
        else
            inplace_merge(__first, __first + i * job_size, __first + (i+1) * job_size, __comp);
    }
}


#endif
