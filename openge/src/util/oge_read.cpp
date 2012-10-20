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
    vector<OGERead *> * cached_allocations, *cached_allocations_cleared;
    Spinlock * allocator_spinlock;
    bool * clean_thread_running;
    virtual void runJob();
};

void OGEReadClearJob::runJob() {
    vector<OGERead *> to_clear;    //thread private copy of to-be-cleared alignments
    to_clear.reserve(100);
    static int ct = 0;
    while (true) {
        allocator_spinlock->lock();
        while(!cached_allocations->empty() && to_clear.size() < 100) {
            OGERead * al = cached_allocations->back();
            cached_allocations->pop_back();
            ct++;
            to_clear.push_back(al);
        }
        
        allocator_spinlock->unlock();
        
        if(to_clear.empty()) break;
        
        for(vector<OGERead *>::iterator i = to_clear.begin(); i!= to_clear.end(); i++)
            (*i)->clear();
        
        allocator_spinlock->lock();
        for(vector<OGERead *>::iterator i = to_clear.begin(); i!= to_clear.end(); i++)
            cached_allocations_cleared->push_back(*i);
        allocator_spinlock->unlock();
        
        to_clear.clear();
    }
    *clean_thread_running = false;
}

OGERead * OGERead::allocate() {
    OGERead * ret = NULL;
    allocator_spinlock.lock();
    
    if(!cached_allocations_cleared.empty()) {
        ret = cached_allocations_cleared.back();
        cached_allocations_cleared.pop_back();
    }
    
    allocator_spinlock.unlock();
    if(!ret) ret = new OGERead();
    //else ret->clear();
    
    return ret;
}

void OGERead::deallocate(OGERead * al) {
    allocator_spinlock.lock();
    cached_allocations.push_back(al);
    allocator_spinlock.unlock();
    
    if(cached_allocations.size() > 100 && !clean_thread_running) {
        clean_thread_running = true;
        OGEReadClearJob * job = new OGEReadClearJob;
        job->allocator_spinlock = &allocator_spinlock;
        job->cached_allocations = & cached_allocations;
        job->clean_thread_running = &clean_thread_running;
        job->cached_allocations_cleared = &cached_allocations_cleared;
        ThreadPool::sharedPool()->addJob(job);
    }
}

void OGERead::clearCachedAllocations() {
    ThreadPool::sharedPool()->waitForJobCompletion();
    allocator_spinlock.lock();
    while(!cached_allocations.empty()) {
        OGERead * al = cached_allocations.back();
        delete(al);
        cached_allocations.pop_back();
    }
    allocator_spinlock.unlock();
}

Spinlock OGERead::allocator_spinlock;
std::vector<OGERead *> OGERead::cached_allocations;
std::vector<OGERead *> OGERead::cached_allocations_cleared;
bool OGERead::clean_thread_running = false;