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
public:
    Spinlock()
    : lock_holder(0) {}
    // for an explanation of the locks used here, see
    // http://stackoverflow.com/questions/1383363/is-my-spin-lock-implementation-correct-and-optimal
    void lock() {
        while (__sync_lock_test_and_set(&lock_holder, 1))
            while(lock_holder)
                asm("nop");
        
    }
    void unlock() {
        __sync_lock_release(&lock_holder);
    }
protected:
    volatile char lock_holder;
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
	
protected:
	static void * thread_start(void * thread_pool);
	ThreadJob * startJob();
	void stopJob(ThreadJob * job);
	
    std::queue<ThreadJob *> jobs;	//protected by jobs_mutex
	std::vector<pthread_t> threads;
	bool threads_exit;
	int jobs_in_process;	//protected by jobs_mutex
	sem_t * job_semaphore;
    sem_t * job_submission_semaphore;
	Spinlock jobs_mutex;
	pthread_mutex_t busy_mutex;
	char sem_name[48];
	char sem_submission_name[48];
    int jobs_current;
};

#endif
