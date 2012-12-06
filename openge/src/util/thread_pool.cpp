/*********************************************************************
 *
 * thread_pool.cpp: Various structures and tools for introducing 
 *                  parallelism in OGE classes.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 7 March 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/


#include "thread_pool.h"

#include <cassert>
#include <string>
#include <iostream>
#include <stdio.h>

#include <fcntl.h>
#include <sys/stat.h>
using namespace std;

// We limit the number of jobs in the queue because in some cases where
// the thread feeding jobs to the pool doesn't have to wait for individual
// jobs to complete, and can submit jobs faster than they can execute, the
// backlog can consume most available memory. To prevent this, we limit the
// job queue size.
const int THREADPOOL_MAX_JOBS_IN_QUEUE = 128;

ThreadPool * ThreadPool::_sharedPool = NULL;

ThreadPool::ThreadPool(int num_threads)
: num_jobs_running(0)
{
    threads_exit.clear();
	if(num_threads == -1 || num_threads == 0)
		num_threads = OGEParallelismSettings::getNumberThreads();

	for(int thread_ctr = 0; thread_ctr < num_threads; thread_ctr++)
	{
		pthread_t thread;
		int error_code = pthread_create(&thread, NULL, ThreadPool::thread_start, this);
		if(0 != error_code) {
			cerr << "Error creating threadpool worker threads. Aborting. (error " << error_code << ")" << endl;
            exit(-1);
        }
		threads.push_back(thread);
	}
}

ThreadPool::~ThreadPool()
{
	//signal threads to return
	threads_exit = true;
	
	//now make the threads check for return signal
    job_queue_mutex.lock();
    job_queue_cond.notify_all();
    job_queue_mutex.unlock();

	//wait for threads to return
	for(size_t thread_ctr = 0; thread_ctr < threads.size(); thread_ctr++)
		pthread_join(threads[thread_ctr], NULL);
}

int ThreadPool::availableCores()
{
    return OGEParallelismSettings::availableCores();
}

//add a job to the queue to be 
bool ThreadPool::addJob(ThreadJob * job)
{
    job_queue_mutex.lock();

	jobs.push(job);

    job_queue_cond.notify_one();
    job_queue_mutex.unlock();

    return true;
}

// When this function is released, it means a thread in the threadpool has had a semaphore
// released, signalling that there is a job available.
// @return A pointer to the task to start, or NULL if the thread should exit.
ThreadJob * ThreadPool::startJob()
{
    ThreadJob * job = NULL;
    
    struct timespec wait_time = {0};
    wait_time.tv_sec = 1;

    job_queue_mutex.lock();

    while(true) {
        //if we are quitting, then allow tht thread to quit
        if(threads_exit.isSet())
            break;
        //if we hit this, there is a job in the queue and we are the only thread that can acquire it
        if(!jobs.empty())
            break;
        
        job_queue_cond.wait(job_queue_mutex);
    }

    jobs_running_mutex.lock();
    num_jobs_running++;
    jobs_running_mutex.unlock();
    
    if(!threads_exit.isSet()) {
        job = jobs.pop();
    }

    job_queue_mutex.unlock();

	if(threads_exit.isSet())
		return NULL;
	
	return job;
}

void ThreadPool::stopJob(ThreadJob * job)
{
    // Determine if all jobs are complete, and there are no more jobs in the queue.
    // If so, we can notify the waitForJobCompletion threads.
    
    jobs_running_mutex.lock();
    num_jobs_running--;
    bool all_jobs_complete = (0 == num_jobs_running) && jobs.empty();
    jobs_running_mutex.unlock();
    
	if(all_jobs_complete) {
        busy_mutex.lock();
        busy_cond.notify_all();
        busy_mutex.unlock();
    }
}

void ThreadPool::waitForJobCompletion()
{
    busy_mutex.lock();
    
    while(true) {
        bool empty = jobs.empty();

        if(empty) break;
        
        busy_cond.wait(busy_mutex);
    }
    
    busy_mutex.unlock();
}

int ThreadPool::numJobs()
{
	return jobs.size();
}

void * ThreadPool::thread_start(void * thread_pool)
{
	ThreadPool * pool = (ThreadPool *)thread_pool;
	while(true)
	{		
		// now pull and do the job
		ThreadJob * job = pool->startJob();
		if(job == NULL)
			return NULL;
		job->runJob();
		pool->stopJob(job);
	}
	
}

ThreadJob::~ThreadJob() {}


ThreadPool * ThreadPool::sharedPool() {
    assert(OGEParallelismSettings::isMultithreadingEnabled());
    if(!_sharedPool)
        _sharedPool = new ThreadPool(OGEParallelismSettings::getNumberThreads());
    return _sharedPool;
}
void ThreadPool::closeSharedPool() {
    if(!_sharedPool)
        return;
    _sharedPool->waitForJobCompletion();
    delete _sharedPool;
    _sharedPool = NULL;
}

#pragma mark OGEParallelismSettings

#include <unistd.h>

int OGEParallelismSettings::m_configured_threads = 0;
bool OGEParallelismSettings::m_multithreading_enabled = false;

int OGEParallelismSettings::availableCores()
{
    
#ifdef OPENMP
    // the openmp function is included here as it is known to be a more reliable way to
    // detect the number of available cores. However, it is only available when compiling
    // with OMP turned on.
	return omp_get_num_procs();
#endif
	
	return sysconf(_SC_NPROCESSORS_ONLN);
}

void OGEParallelismSettings::setNumberThreads(int threads)
{
    m_configured_threads = threads;
}

int OGEParallelismSettings::getNumberThreads()
{
    if(m_configured_threads == 0)
        return availableCores();
    else
        return m_configured_threads;
}

void OGEParallelismSettings::disableMultithreading()
{
    m_multithreading_enabled = false;
}

void OGEParallelismSettings::enableMultithreading()
{
    m_multithreading_enabled = true;
}

mutex::mutex() {
    int ret = pthread_mutex_init(&m, NULL);
    if(0 != ret) {
        cerr << "Error creating mutex (error " << ret << ")." << endl;
        exit(-1);
    }
}
mutex::~mutex() {
    int ret = pthread_mutex_destroy(&m);
    if(0 != ret) {
        cerr << "Error destroying mutex (error " << ret << ")." << endl;
        exit(-1);
    }
}
void mutex::lock() {
    int ret = pthread_mutex_lock(&m);
    if(0 != ret) {
        cerr << "Error locking mutex (error " << ret << ")." << endl;
        exit(-1);
    }
}
bool mutex::try_lock() {
    int ret = pthread_mutex_trylock(&m);
    
    if(ret == 16)   //16 = EBUSY
        return false;
    if(0 != ret) {
        cerr << "Error locking mutex (error " << ret << ")." << endl;
        exit(-1);
    }
    return true;
}
void mutex::unlock() {
    int ret = pthread_mutex_unlock(&m);
    if(0 != ret) {
        cerr << "Error unlocking mutex (error " << ret << ")." << endl;
        exit(-1);
    }
}

condition_variable::condition_variable() {
    int ret = pthread_cond_init(&c, NULL);
    if(0 != ret) {
        cerr << "Error creating CV (error " << ret << ")." << endl;
        exit(-1);
    }
}
condition_variable::~condition_variable() {
    int ret = pthread_cond_destroy(&c);
    if(0 != ret) {
        cerr << "Error destroying CV (error " << ret << ")." << endl;
        exit(-1);
    }
}
void condition_variable::notify_one() {
    int ret = pthread_cond_signal(&c);
    if(0 != ret) {
        cerr << "Error signalling CV (error " << ret << ")." << endl;
        exit(-1);
    }
}
void condition_variable::notify_all() {
    int ret = pthread_cond_broadcast(&c);
    if(0 != ret) {
        cerr << "Error broadcasting CV (error " << ret << ")." << endl;
        exit(-1);
    }
}
void condition_variable::wait(mutex & m) {
    int ret = pthread_cond_wait(&c, &m.native_handle());
    if(0 != ret) {
        cerr << "Error waiting for CV (error " << ret << ")." << endl;
        exit(-1);
    }
}


