//
//  bamtools_thread_pool.cpp
//  BamTools
//
//  Created by Lee Baker on 3/7/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#include "api/internal/utils/BamThreadPool.h"
#include "api/BamParallelismSettings.h"
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

BamThreadPool::BamThreadPool() :
threads_exit(false),
jobs_current(0)
{
	int num_threads = BamParallelismSettings::getNumberThreads();

	jobs_in_process = 0;
  int32_t sem_id = 0xffffffff & (int64_t) this;

  sprintf(sem_name, "bam_tp_%x",sem_id);

    if(0 != sem_unlink(sem_name) && errno != ENOENT && errno != EINVAL)
       perror("Error unlinking BT sem_name");
  job_semaphore = sem_open(sem_name, O_CREAT | O_EXCL,0700,0);

	if(job_semaphore == SEM_FAILED &&  0 != errno) {
		perror("Error opening threadpool semaphore");
		assert(0);
	}

  sprintf(sem_submission_name, "bam_tpjs_%x",sem_id);
  if(0 != sem_unlink(sem_submission_name) && errno != ENOENT && errno != EINVAL)
    perror("Error unlinking sem_submission_name");
  job_submission_semaphore = sem_open(sem_submission_name, O_CREAT | O_EXCL,0700,THREADPOOL_MAX_JOBS_IN_QUEUE);

	if(job_submission_semaphore == SEM_FAILED &&  0 != errno) {
		perror("Error opening threadpool job submission semaphore");
		assert(0);
	}

	int retval = pthread_mutex_init(&jobs_mutex,NULL);

	if(0 != retval)
		perror("Error opening threadpool mutex");
	assert(0 == retval);

	retval = pthread_mutex_init(&busy_mutex,NULL);

	if(0 != retval)
		perror("Error opening threadpool busy mutex");
	assert(0 == retval);

	for(int thread_ctr = 0; thread_ctr < num_threads; thread_ctr++)
	{
		pthread_t thread;
		int error_code = pthread_create(&thread, NULL, BamThreadPool::thread_start, this);
		if(0 != error_code)
			perror("Error creating threadpool worker threads");
		assert(0 == error_code);
		threads.push_back(thread);
	}
}

BamThreadPool::~BamThreadPool()
{
	//signal threads to return
	threads_exit = true;
	
	//now make the threads check for return signal
	if(0 != sem_post(job_semaphore))
    perror("Error posting job semaphore");
	
	//wait for threads to return
	for(size_t thread_ctr = 0; thread_ctr < threads.size(); thread_ctr++)
		pthread_join(threads[thread_ctr], NULL);
	
	if(0 != pthread_mutex_destroy(&jobs_mutex))
  {
    perror("Error destroying jobs mutex");
    assert(0);
  }
  
	if(0 != pthread_mutex_destroy(&busy_mutex))
  {
    perror("Error destroying busy mutex");
    assert(0);
  }
  if(0 != sem_close(job_semaphore))
    perror("Error closing job_semaphore");
	if(0 != sem_unlink(sem_name))
    perror("Error unlinking sem_name");
  if(0 != sem_close(job_submission_semaphore))
    perror("Error closing job_submission_semaphore");
  if(0 != sem_unlink(sem_submission_name) && errno != EINVAL)
    perror("Error unlinking sem_submission_name");
}

//add a job to the queue to be 
bool BamThreadPool::addJob(BamThreadJob * job)
{
  if(0 != sem_wait(job_submission_semaphore))
    perror("Error waiting for job submission semaphore (add job)");
  
	if(0 != pthread_mutex_lock(&jobs_mutex))
  {
    perror("Error locking jobs mutex");
    assert(0);
  }
  if(jobs_current == 0)
    pthread_mutex_lock(&busy_mutex);
  jobs_current++;
  
	jobs.push(job);
	if(0 != pthread_mutex_unlock(&jobs_mutex))
  {
   perror("Error unlocking jobs mutex");
   assert(0);
  }
	if(0 != sem_post(job_semaphore))
    perror("Error posting job semaphore");
	return true;
}

// When this function is released, it means a thread in the threadpool has had a semaphore
// released, signalling that there is a job available.
// @return A pointer to the task to start, or NULL if the thread should exit.
BamThreadJob * BamThreadPool::startJob()
{
    BamThreadJob * job  = NULL;
    do {
        if(0 != sem_wait(job_semaphore))
            perror("Error waiting for job submission semaphore (start job)");

        if(threads_exit)
        {
            if(0 != sem_post(job_semaphore))
                perror("Error posting job semaphore");
            return NULL;
        }

        if(0 != pthread_mutex_lock(&jobs_mutex))
        {
            perror("Error locking jobs mutex");
            assert(0);
        }

        if(jobs.size() > 0) {
            //assert(jobs.size() > 0);
            job = jobs.front();
            jobs.pop();
            jobs_in_process++;
        }

        if(0 != pthread_mutex_unlock(&jobs_mutex))
        {
            perror("Error unlocking jobs mutex");
            assert(0);
        }
    } while(!job);
	
	return job;
}

void BamThreadPool::stopJob(BamThreadJob * job)
{
	if(0 != pthread_mutex_lock(&jobs_mutex))
  {
    perror("Error locking jobs mutex");
    assert(0);
  }
  
	jobs_in_process--;
  jobs_current--;
	if(jobs_current == 0)
		pthread_mutex_unlock(&busy_mutex);
	if(0 != pthread_mutex_unlock(&jobs_mutex))
  {
    perror("Error unlocking jobs mutex");
    assert(0);
  }
  if(0 != sem_post(job_submission_semaphore))
    perror("Error posting job_submission_semaphore");
}

void BamThreadPool::waitForJobCompletion()
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

int BamThreadPool::numJobs()
{
	return jobs.size();
}

void * BamThreadPool::thread_start(void * thread_pool)
{
	BamThreadPool * pool = (BamThreadPool *)thread_pool;
	while(true)
	{		
		// now pull and do the job
		BamThreadJob * job = pool->startJob();
		if(job == NULL)
			return NULL;
		job->runJob();
		pool->stopJob(job);
    delete job;
	}
	
}

BamThreadJob::~BamThreadJob() {}
