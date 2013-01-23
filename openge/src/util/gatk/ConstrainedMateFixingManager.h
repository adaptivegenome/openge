/*********************************************************************
 *
 * ConstrainedMateFixingManager.h: Port of GATK's ConstrainedMateFixingManager.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 5 May 2012
 *
 *********************************************************************
 *
 * This file has been ported from GATK's implementation in Java, and
 * is released under the Virginia Tech Non-Commercial Purpose License.
 * A copy of this license has been provided in  the openge/ directory.
 * 
 * The original file, ConstrainedMateFixinigManager.java, was released 
 * under the following license:
 *
 * Copyright (c) 2010 The Broad Institute. 
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef OpenGE_ConstrainedMateFixingManager_h
#define OpenGE_ConstrainedMateFixingManager_h

#include "GenomeLoc.h"
#include "GenomeLocParser.h"

#include "../bamtools/Sort.h"

#include "../oge_read.h"
#include "../thread_pool.h"

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
    
    typedef struct {
        OGERead * read;
        bool readWasModified;
        bool canFlush;
    } cmfm_read_t;
    SynchronizedBlockingQueue<cmfm_read_t> addReadQueue;
    pthread_t add_read_thread;
    mutex add_read_lock;

    GenomeLoc * lastLocFlushed;
    
    int counter;
    GenomeLocParser * loc_parser;
    
    class SAMRecordHashObject {
    public:
        OGERead * record;
        bool wasModified;
        
        SAMRecordHashObject(OGERead * record, bool wasModified) 
        : record(record), wasModified(wasModified)
        { }
        SAMRecordHashObject(const SAMRecordHashObject & o)
        : record(o.record), wasModified(o.wasModified) {}
        SAMRecordHashObject() {}
    };
    
    /** read.name -> records */
    std::map<std::string, SAMRecordHashObject> forMateMatching;
    typedef std::set<OGERead *, BamTools::Algorithms::Sort::ByPosition> waitingReads_t;
    waitingReads_t waitingReads;

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
    ~ConstrainedMateFixingManager();

    bool canMoveReads(const GenomeLoc & earliestPosition) const;
    
private:
    bool noReadCanMoveBefore(int pos, const OGERead & addedRead) const;
    
public:
    void addReads(const std::vector<OGERead *> & newReads, const std::set<OGERead *> & modifiedReads);
    void addRead(OGERead * newRead, bool readWasModified, bool canFlush);
    // In some situations, the read queue can fill up faster than we can process reads.
    // The class feeding reads into CMFM should throttle itself to ensure we don't end up
    // eating all system memory by making this queue enormous.
    //
    // Ideally we would check how much free memory we have in the system- something like
    // getrsuage()'s ru_maxrss, but unfortunately Linux has chosen not to implement this,
    // so it isn't a workable solution.
    //
    // 10000 is just a number I picked- there really should be more analysis behind this.
    // This is not a hard limit.
    bool isReadAddQueueFull() { return addReadQueue.size() > 10000; }
private:
    static void * addread_threadproc(void * data);
    void addReadInternal( OGERead * newRead, bool readWasModified, bool canFlush);
    OGERead * remove(waitingReads_t & treeSet);
    void writeRead( OGERead * read) const;
    
    /**
     * @param read  the read
     * @return true if the read shouldn't be moved given the constraints of this SAMFileWriter
     */
public:
    bool iSizeTooBigToMove(const OGERead & read) const;
    static bool iSizeTooBigToMove(const OGERead & read, int maxInsertSizeForMovingReadPairs);
    
private:
    void purgeUnmodifiedMates();
    bool pairedReadIsMovable(const OGERead & read) const;
    
public:
    void close();
};

#endif
