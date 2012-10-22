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
#include <errno.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/errno.h>
using namespace std;

// We limit the number of jobs in the queue because in some cases where
// the thread feeding jobs to the pool doesn't have to wait for individual
// jobs to complete, and can submit jobs faster than they can execute, the
// backlog can consume most available memory. To prevent this, we limit the
// job queue size.
const int THREADPOOL_MAX_JOBS_IN_QUEUE = 128;

ThreadPool * ThreadPool::_sharedPool = NULL;

ThreadPool::ThreadPool(int num_threads) :
threads_exit(false),
jobs_current(0)
{
	if(num_threads == -1 || num_threads == 0)
		num_threads = availableCores();
    
	jobs_in_process = 0;
    
    if(0 != pthread_cond_init(&job_queue_cond, NULL))
        perror("Error creating job queue cond variable");
    
    if(0 != pthread_mutex_init(&job_queue_mutex, NULL))
        perror("Error creating job queue cond variable");

	if(0 != pthread_mutex_init(&busy_mutex,NULL))
		perror("Error creating threadpool busy mutex");
    
	for(int thread_ctr = 0; thread_ctr < num_threads; thread_ctr++)
	{
		pthread_t thread;
		int error_code = pthread_create(&thread, NULL, ThreadPool::thread_start, this);
		if(0 != error_code)
			perror("Error creating threadpool worker threads");
		assert(0 == error_code);
		threads.push_back(thread);
	}
}

ThreadPool::~ThreadPool()
{
	//signal threads to return
	threads_exit = true;
	
	//now make the threads check for return signal
	if(0 != pthread_cond_broadcast(&job_queue_cond))
        perror("Error sending threadpool end signal");
	
	//wait for threads to return
	for(size_t thread_ctr = 0; thread_ctr < threads.size(); thread_ctr++)
		pthread_join(threads[thread_ctr], NULL);
    
    int error = pthread_mutex_destroy(&busy_mutex);
	if(0 != error)
        cerr << "Error destroying TP busy mutex (error " << error << ")." << endl ;

    error = pthread_cond_destroy(&job_queue_cond);
    if(0 != error)
        cerr << "Error destroying TP jobpool condition var (error " << error << ")." << endl ;

    error = pthread_mutex_destroy(&job_queue_mutex);
    if(0 != error)
        cerr << "Error destroying TP job queue mutex (error " << error << ")." << endl ;
}

int ThreadPool::availableCores()
{
    return OGEParallelismSettings::availableCores();
}

//add a job to the queue to be 
bool ThreadPool::addJob(ThreadJob * job)
{
    jobs_mutex.lock();

    if(jobs_current == 0)
        pthread_mutex_lock(&busy_mutex);
    jobs_current++;
    
	jobs.push(job);
    jobs_mutex.unlock();

    if(0 != pthread_cond_signal(&job_queue_cond))
        perror("Error signalling job threads on new job.");
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

    if(0 != pthread_mutex_lock(&job_queue_mutex))
        perror("Error locking job queue mutex in startJob()");
    while(!threads_exit && jobs.empty()) {
        if(0 != pthread_cond_wait(&job_queue_cond, &job_queue_mutex) && errno != 0 && errno != ETIMEDOUT)
        //if(0 != pthread_cond_timedwait(&job_queue_cond, &job_queue_mutex, &wait_time) && errno != 0 && errno != ETIMEDOUT)
            perror("Error waiting for job in startJob()");
    }
    
    if(!threads_exit) {
        jobs_mutex.lock();
        job = jobs.front();
        jobs.pop();
        jobs_in_process++;
        jobs_mutex.unlock();
    }
        
    if(0 != pthread_mutex_unlock(&job_queue_mutex))
        perror("Error unlocking job queue mutex in startJob()");

	if(threads_exit)
		return NULL;
	
	return job;
}

void ThreadPool::stopJob(ThreadJob * job)
{
	jobs_mutex.lock();
    
	jobs_in_process--;
    jobs_current--;
	if(jobs_current == 0)
		pthread_mutex_unlock(&busy_mutex);
	jobs_mutex.unlock();
}

void ThreadPool::waitForJobCompletion()
{
	if(0 != pthread_mutex_lock(&busy_mutex))
    {
        perror("Error locking busy mutex");
        assert(0);
    }
    if(0 != pthread_mutex_unlock(&busy_mutex))
    {
        perror("Error unlocking busy mutex");
        assert(0);
    }
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
