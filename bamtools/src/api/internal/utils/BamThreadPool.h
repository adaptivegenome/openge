//
//  bamtools_thread_pool.h
//  BamTools
//
//  Created by Lee Baker on 3/7/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//
// ***************************************************************************
// bamtools_thread_poo.h (c) 2012 Lee C. Baker
// ---------------------------------------------------------------------------
// Last modified: 9 March 2012
// ---------------------------------------------------------------------------
// Provides a thread pool designed to execute provided jobs in parallel. Jobs
// should be implemented in a subclass of BamThreadJob.
// ***************************************************************************

#ifndef BAMTOOLS_THREAD_POOL_H
#define BAMTOOLS_THREAD_POOL_H

#include <vector>
#include <queue>
#include <pthread.h>
#include <semaphore.h>
using namespace std;

class BamThreadPool;

// The BamThreadJob abstract class provides a way to provide jobs to 
// the BamThreadPool thread pool. To use, create a thread pool, subclass BamThreadJob,
// implement the runJob method of your subclass, and pass an instance to the pool's
// addJob method.
class BamThreadJob
{
	friend class BamThreadPool;
public:
  virtual ~BamThreadJob();
	virtual void runJob() = 0;
protected:
};

class BamThreadPool
{
	friend class BamThreadJob;
public:
	BamThreadPool();
	virtual ~BamThreadPool();
	
	bool addJob(BamThreadJob * job);
	int numJobs();
	void waitForJobCompletion();
	
protected:
	static void * thread_start(void * thread_pool);
	BamThreadJob * startJob();
	void stopJob(BamThreadJob * job);
	
	queue<BamThreadJob *> jobs;	//protected by jobs_mutex
	vector<pthread_t> threads;
	bool threads_exit;
	int jobs_in_process;	//protected by jobs_mutex
	sem_t * job_semaphore;
    sem_t * job_submission_semaphore;
	pthread_mutex_t jobs_mutex;
	pthread_mutex_t busy_mutex;
	char sem_name[48];
	char sem_submission_name[48];
  int jobs_current;
};

class Spinlock
{
public:
    Spinlock()
    : lock_holder(0) {}
    // for an explanation of the locks used here, see
    // http://stackoverflow.com/questions/1383363/is-my-spin-lock-implementation-correct-and-optimal
    void lock() {
        while (__sync_lock_test_and_set(&lock_holder, 1)) while(lock_holder);
        
    }
    void unlock() {
        __sync_synchronize();
        lock_holder = 0;
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
    queue<T> q;
};

#endif
