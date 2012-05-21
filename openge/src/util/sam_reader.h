//
//  SamReader.h
//  BamTools
//
//  Created by Lee Baker on 3/16/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#ifndef BamTools_SamReader_h
#define BamTools_SamReader_h
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
    bool parsed;
    SamLine() : line((char *)3), al((BamTools::BamAlignment *)4), parsed(false) {}
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

    //std::vector<BamTools::BamAlignment> alignments;
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
    
    SynchronizedQueue<SamLine * > jobs_for_workers;
    std::vector<pthread_t> worker_threads;
    bool workers_finished;
    Spinlock worker_jobs_lock;  //we need to syncronize access to the queue, just to make sure that two threads don't try to pull out the last job (and one get invalid data when it should instead be waiting again)
};

#endif
