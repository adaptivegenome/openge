//
//  ConstrainedMateFixingManager.h
//  OpenGE
//
//  Created by Lee Baker on 5/5/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#ifndef OpenGE_ConstrainedMateFixingManager_h
#define OpenGE_ConstrainedMateFixingManager_h

#include "api/BamWriter.h"
#include "api/SamHeader.h"
#include "api/BamAlignment.h"
#include "GenomeLoc.h"
#include "GenomeLocParser.h"

#include <list>
#include <map>
#include <set>
#include <string>

class LocalRealignment;

class ConstrainedMateFixingManager {
    
private:
    static const bool DEBUG = false;
    
    /** How often do we check whether we want to emit reads? */
    const static int EMIT_FREQUENCY = 1000;
    
    /** The place where we ultimately write out our records */
    LocalRealignment * output_module;
    
    /**
     * what is the maximum isize of a pair of reads that can move?  Reads with isize > this value
     * are assumes to not be allowed to move in the incoming read stream.
     */
    int maxInsertSizeForMovingReadPairs;
    
    /**
     * How much could a single read move in position from its original position?
     */
    int MAX_POS_MOVE_ALLOWED;
    
    /**
     * How many reads should we store in memory before flushing the queue?
     */
    int MAX_RECORDS_IN_MEMORY;

    GenomeLoc * lastLocFlushed;
    
    int counter;
    GenomeLocParser * loc_parser;
    
    class SAMRecordHashObject {
    public:
        BamTools::BamAlignment * record;
        bool wasModified;
        
        SAMRecordHashObject(BamTools::BamAlignment * record, bool wasModified) 
        : record(record), wasModified(wasModified)
        { }
        SAMRecordHashObject(const SAMRecordHashObject & o)
        : record(o.record), wasModified(o.wasModified) {}
        SAMRecordHashObject() {}
    };
    
    /** read.name -> records */
    std::map<std::string, SAMRecordHashObject> forMateMatching;
    std::set<BamTools::BamAlignment *> waitingReads;
    
private:
    BamTools::BamAlignment * remove(std::set<BamTools::BamAlignment *> & treeSet);
    
    
    //private SimpleTimer timer = new SimpleTimer("ConstrainedWriter");
    //private long PROGRESS_PRINT_FREQUENCY = 10 * 1000;             // in milliseconds
    //private long lastProgressPrintTime = -1;                       // When was the last time we printed progress log?
    
    
    /**
     *
     * @param writer                                 actual writer
     * @param genomeLocParser                        the GenomeLocParser object
     * @param maxInsertSizeForMovingReadPairs        max insert size allowed for moving pairs
     * @param maxMoveAllowed                         max positional move allowed for any read
     * @param maxRecordsInMemory                     max records to keep in memory
     */
public:
    ConstrainedMateFixingManager( LocalRealignment * writer, const int maxInsertSizeForMovingReadPairs, const int maxMoveAllowed, const int maxRecordsInMemory,GenomeLocParser * loc_parser);
    
    int getNReadsInQueue();
    
    bool canMoveReads(const GenomeLoc & earliestPosition);
    
private:
    bool noReadCanMoveBefore(int pos, BamTools::BamAlignment * addedRead);
    
public:
    void addReads(std::vector<BamTools::BamAlignment *> newReads, std::set<BamTools::BamAlignment *> modifiedReads);
    void addRead(BamTools::BamAlignment * newRead, bool readWasModified, bool canFlush = true);
    
private:
    void writeRead(BamTools::BamAlignment * read);
    
    /**
     * @param read  the read
     * @return true if the read shouldn't be moved given the constraints of this SAMFileWriter
     */
public:
    bool iSizeTooBigToMove(BamTools::BamAlignment * read);
    static bool iSizeTooBigToMove(BamTools::BamAlignment * read, int maxInsertSizeForMovingReadPairs);
    
private:
    void purgeUnmodifiedMates();
    bool pairedReadIsMovable(BamTools::BamAlignment * read);
    
public:
    void close();
};

#endif
