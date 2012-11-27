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

ThreadPool::ThreadPool(int num_threads) :
threads_exit(false),
jobs_current(0)
{
	if(num_threads == -1 || num_threads == 0)
		num_threads = OGEParallelismSettings::getNumberThreads();
    
	jobs_in_process = 0;
    
    int ret = pthread_cond_init(&job_queue_cond, NULL);
    if(0 != ret) {
        cerr << "Error creating job queue cond variable. Aborting. (error " << ret << ")" << endl;
        exit(-1);
    }
    
    ret = pthread_cond_init(&busy_cond, NULL);
    if(0 != ret) {
        cerr << "Error creating busy cond variable. Aborting. (error " << ret << ")" << endl;
        exit(-1);
    }
    
    ret = pthread_mutex_init(&job_queue_mutex, NULL);
    if(0 != ret) {
        cerr << "Error creating job queue cond variable. Aborting. (error " << ret << ")" << endl;
        exit(-1);
    }

    ret = pthread_mutex_init(&busy_mutex,NULL);
	if(0 != ret) {
		cerr << "Error creating threadpool busy mutex. Aborting. (error " << ret << ")" << endl;
        exit(-1);
    }

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
    int error = pthread_cond_broadcast(&job_queue_cond);
	if(0 != error)
        cerr << "Error sending threadpool end signal. (error " << error << ")" << endl;
	
	//wait for threads to return
	for(size_t thread_ctr = 0; thread_ctr < threads.size(); thread_ctr++)
		pthread_join(threads[thread_ctr], NULL);
    
    error = pthread_mutex_destroy(&busy_mutex);
	if(0 != error)
        cerr << "Error destroying TP busy mutex (error " << error << ")." << endl ;

    error = pthread_cond_destroy(&job_queue_cond);
    if(0 != error)
        cerr << "Error destroying TP jobpool condition var (error " << error << ")." << endl ;

    error = pthread_cond_destroy(&busy_cond);
    if(0 != error)
        cerr << "Error destroying TP busy condition var (error " << error << ")." << endl ;

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
    int ret = pthread_mutex_lock(&job_queue_mutex);
    if(0 != ret) {
        cerr << "Error locking job queue mutex in addJob(). Aborting. (error " << ret << ")" << endl;
        exit(-1);
    }

    jobs_mutex.lock();
    
    jobs_current++;
    
	jobs.push(job);
    jobs_mutex.unlock();

    ret = pthread_cond_signal(&job_queue_cond);
    if(0 != ret) {
        cerr << "Error signalling job threads on new job. (error " << ret << ")." << endl ;
        exit(-1);
    }
    ret = pthread_mutex_unlock(&job_queue_mutex);
    if(0 != ret) {
        cerr << "Error unlocking job queue mutex in addJob(). Aborting. (error " << ret << ")." << endl;
        exit(-1);
    }
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

    int ret = pthread_mutex_lock(&job_queue_mutex);
    if(0 != ret) {
        cerr << "Error locking job queue mutex in startJob(). Aborting. (error " << ret << ")" << endl;
        exit(-1);
    }

    while(!threads_exit && jobs.empty()) {
        int retval = pthread_cond_wait(&job_queue_cond, &job_queue_mutex);
        if(0 != retval) {
            cerr << "Error waiting for job in startJob(). Aborting. (error " << retval << ")." << endl ;
            exit(-1);
        }
    }
    
    if(!threads_exit) {
        jobs_mutex.lock();
        job = jobs.front();
        jobs.pop();
        jobs_in_process++;
        jobs_mutex.unlock();
    }
        
    ret = pthread_mutex_unlock(&job_queue_mutex);
    if(0 != ret) {
        cerr << "Error unlocking job queue mutex in startJob(). Aborting. (error " << ret << ")." << endl;
        exit(-1);
    }

	if(threads_exit)
		return NULL;
	
	return job;
}

void ThreadPool::stopJob(ThreadJob * job)
{
	jobs_mutex.lock();
    
	jobs_in_process--;
    jobs_current--;
	if(jobs_current == 0) {
        int ret = pthread_mutex_lock(&busy_mutex);
        if(0 != ret) {
            cerr << "Error locking job queue mutex in addJob(). Aborting. (error " << ret << ")" << endl;
            exit(-1);
        }
        ret = pthread_cond_signal(&busy_cond);
        if(0 != ret) {
            cerr << "Error signalling job threads on new job. (error " << ret << ")." << endl ;
            exit(-1);
        }
        ret = pthread_mutex_unlock(&busy_mutex);
        if(0 != ret) {
            cerr << "Error unlocking job queue mutex in addJob(). Aborting. (error " << ret << ")." << endl;
            exit(-1);
        }
    }
	jobs_mutex.unlock();
}

void ThreadPool::waitForJobCompletion()
{
    int ret = pthread_mutex_lock(&busy_mutex);
    if(0 != ret) {
        cerr << "Error locking job queue busy in waitForJobCompletion(). Aborting. (error " << ret << ")" << endl;
        exit(-1);
    }
    
    while(true) {
        jobs_mutex.lock();
        bool empty = jobs.empty();
        jobs_mutex.unlock();

        if(empty) break;
        
        int retval = pthread_cond_wait(&busy_cond, &busy_mutex);
        if(0 != retval) {
            cerr << "Error waiting for busy condition variable in waitForJobCompletion(). Aborting. (error " << retval << ")." << endl ;
            exit(-1);
        }
    }
    
    ret = pthread_mutex_unlock(&busy_mutex);
    if(0 != ret) {
        cerr << "Error unlocking busy mutex in waitForJobCompletion(). Aborting. (error " << ret << ")." << endl;
        exit(-1);
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


