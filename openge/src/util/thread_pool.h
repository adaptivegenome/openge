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
    T pop() {
        bool success = false;
        T ret;
        while(!success) {
            this->lock.lock();
            if(!this->q.empty()) {
                ret = this->q.front();
                this->q.pop();
                success = true;
            }
            this->lock.unlock();
            if(!success)
                usleep(20000);
        }

        return ret;
    }
protected:
        
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
    static ThreadPool * sharedPool() { if(!_sharedPool) _sharedPool = new ThreadPool(OGEParallelismSettings::getNumberThreads()); return _sharedPool; }
    static void closeSharedPool() {  _sharedPool->waitForJobCompletion(); delete _sharedPool; }
	
protected:
    static ThreadPool * _sharedPool;
	static void * thread_start(void * thread_pool);
	ThreadJob * startJob();
	void stopJob(ThreadJob * job);
    
    pthread_cond_t job_queue_cond, busy_cond;
    pthread_mutex_t job_queue_mutex;
	
    std::queue<ThreadJob *> jobs;	//protected by jobs_mutex
	std::vector<pthread_t> threads;
	bool threads_exit;
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
