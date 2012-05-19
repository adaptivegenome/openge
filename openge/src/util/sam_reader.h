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

typedef struct SamLine {
    std::string line;
    BamTools::BamAlignment * al;
    SamReader * reader;
    bool parsed;
    SamLine() : al(NULL), parsed(false) {}
} SamLine_t;

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
    BamTools::BamAlignment * ParseAlignment(const std::string & line_s);
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
    
    // multithreading variables
    bool loaded;
    bool finished;
    SynchronizedQueue<SamLine_t *> jobs;
    ThreadPool pool;
    
    static void * LineGenerationThread(void * data);
    pthread_t line_generation_thread;
};

#endif
