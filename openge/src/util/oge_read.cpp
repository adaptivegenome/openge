/*********************************************************************
 *
 * oge_read.cpp:  Main class for storing read data.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 17 Oct 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial
 * Purpose License. A copy of this license has been provided in
 * the openge/ directory.
 *
 *********************************************************************/

#include "oge_read.h"
#include "thread_pool.h"

#include <vector>

using std::vector;

//////////////
// cache allocations for performance

// this class should be a subclass of BamAlignment, but the way the API exports are structured,
// we would have to export more classes with the API. Instead, we'll just do this temporarily
// until bamtools is replaced.
class OGEReadClearJob : public ThreadJob {
public:
    SynchronizedQueue<OGERead *> * cached_allocations, *cached_allocations_cleared;
    Spinlock * allocator_spinlock;
    bool * clean_thread_running;
    virtual void runJob();
};

void OGEReadClearJob::runJob() {
    while(!cached_allocations->empty()) {
        OGERead * al = cached_allocations->pop();
        al->clear();
        cached_allocations_cleared->push(al);
    }
    *clean_thread_running = false;
}

OGERead * OGERead::allocate() {
    allocate_lock.lock();
    OGERead * read = NULL;
    if(!cached_allocations_cleared.empty()) {
        read = cached_allocations_cleared.pop();
    } else {
        read = new OGERead();
    }
    allocate_lock.unlock();
    return read;
}

void OGERead::deallocate(OGERead * al) {
    cached_allocations.push(al);
    
    if(cached_allocations.size() > 100 && !clean_thread_running) {
        clean_thread_running = true;
        OGEReadClearJob * job = new OGEReadClearJob;
        job->cached_allocations = &cached_allocations;
        job->clean_thread_running = &clean_thread_running;
        job->cached_allocations_cleared = &cached_allocations_cleared;
        ThreadPool::sharedPool()->addJob(job);
    }
}

void OGERead::clearCachedAllocations() {
    ThreadPool::sharedPool()->waitForJobCompletion();

    while(!cached_allocations.empty()) {
        OGERead * al = cached_allocations.pop();
        delete(al);
    }
    
    while(!cached_allocations_cleared.empty()) {
        OGERead * al = cached_allocations_cleared.pop();
        delete(al);
    }
}

SynchronizedQueue<OGERead *> OGERead::cached_allocations;
SynchronizedQueue<OGERead *> OGERead::cached_allocations_cleared;
bool OGERead::clean_thread_running = false;
Spinlock OGERead::allocate_lock;
