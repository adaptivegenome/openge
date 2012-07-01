/*********************************************************************
 *
 * sam_reader.h: A multithreaded sequential SAM file reader.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 16 March 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#ifndef OGE_SAMREADER_H
#define OGE_SAMREADER_H
#include <string>
#include <vector>
#include <api/BamAlignment.h>
#include <api/SamHeader.h>
#include <iostream>
#include <queue>
#include "thread_pool.h"

class BamThreadPool;

class SamReader;

class SamLine {
public:
    BamTools::BamAlignment * al;
    char * line;
    char line_static[640];  //for optimization, we statically allocate a block so we can avoid allocation for small lines. If line == NULL, then we are using line_static instead of line.
    bool parsed;
    SamLine() : al(NULL), line(NULL), parsed(false) {}
};

// SamReader is capable of sequentially reading a SAM file. It doesn't support
// most of the features that BamReader does, only enough to support converting SAM
// files to the BAM format.
class SamReader
{
public:
    SamReader();
    bool Open(const std::string & filename);
    bool Close();
    bool IsLoaded();
    // retrieves next available alignment
    bool GetNextAlignment(BamTools::BamAlignment& alignment);
    BamTools::BamAlignment * GetNextAlignment();
    
    // returns the current file's header data
    BamTools::SamHeader GetHeader(void) const;
    // get reference data
    const BamTools::RefVector & GetReferenceData(void);
    //read a single line of a SAM file
    BamTools::BamAlignment * ParseAlignment(const char * line_s);
protected:
    // retrieves BAM alignment under file pointer
    // (does no overlap checking or character data parsing)
    BamTools::BamAlignment * LoadNextAlignment();
    
    // retrieves header text from SAM file
    void LoadHeaderData(void);

    std::ifstream file;
    BamTools::SamHeader header;
    BamTools::RefVector m_refData;
    std::string filename;
    
    // multithreading variables
    bool loaded;
    bool finished;
    SynchronizedQueue<SamLine *> jobs;
    
    static void * LineGenerationThread(void * data);
    static void * LineWorkerThread(void * reader_p);
    pthread_t line_generation_thread;
    
    sem_t * sam_worker_sem;
    char sam_worker_sem_name[32];
    SynchronizedQueue<SamLine * > jobs_for_workers;
    std::vector<pthread_t> worker_threads;
    int num_worker_threads;
    int active_workers; //< Keep track of the number of threads processing lines.
    bool workers_finished;
    int lines_since_last_sem_unlock;
    Spinlock worker_jobs_lock;  //we need to syncronize access to the queue, just to make sure that two threads don't try to pull out the last job (and one get invalid data when it should instead be waiting again)
};

#endif
