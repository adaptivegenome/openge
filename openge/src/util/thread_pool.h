//
//  thread_pool.h
//  OpenGE
//
//  Created by Lee Baker on 4/6/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//
// ***************************************************************************
// thread_poo.h (c) 2012 Lee C. Baker
// ---------------------------------------------------------------------------
// Last modified: 6 April 2012
// ---------------------------------------------------------------------------
// Provides a thread pool designed to execute provided jobs in parallel. Jobs
// should be implemented in a subclass of ThreadJob.
// ***************************************************************************

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <queue>
#include <pthread.h>
#include <semaphore.h>
using namespace std;

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
	
	queue<ThreadJob *> jobs;	//protected by jobs_mutex
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


template <class T>
class SynchronizedQueue
{
public: 
    SynchronizedQueue() 
    : lock_holder(0)
    { }
    
    void push(const T & item) {
        lock();
        q.push(item);
        unlock();
    }
    
    size_t size() {
        lock();
        size_t r = q.size();
        unlock();
        return r;
    }
    
    bool empty() {
        lock();
        bool e = q.empty();
        unlock();
        return e;
    }
    
    T pop() {
        lock();
        T ret = q.front();
        q.pop();
        unlock();
        return ret;
    }
    
protected:
    // for an explanation of the locks used here, see
    // http://stackoverflow.com/questions/1383363/is-my-spin-lock-implementation-correct-and-optimal
    void lock() {
        while (__sync_lock_test_and_set(&lock_holder, 1)) while(lock_holder);
        
    }
    void unlock() {
        __sync_synchronize();
        lock_holder = 0;
    }
    volatile char lock_holder;
    queue<T> q;
};

#endif
