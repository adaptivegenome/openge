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

#include "api/BamWriter.h"
#include "api/SamHeader.h"
#include "api/BamAlignment.h"
#include "GenomeLoc.h"
#include "GenomeLocParser.h"

#include "api/algorithms/Sort.h"

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
    //std::less<BamTools::BamAlignment *> cmp;
    BamTools::Algorithms::Sort::AlignmentPtrSortBase cmp;
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
    void addReads(const std::vector<BamTools::BamAlignment *> & newReads, const std::set<BamTools::BamAlignment *> & modifiedReads);
    void addRead(BamTools::BamAlignment * newRead, bool readWasModified, bool canFlush = true);
    
private:
    void writeRead(BamTools::BamAlignment * read);
    
    /**
     * @param read  the read
     * @return true if the read shouldn't be moved given the constraints of this SAMFileWriter
     */
public:
    bool iSizeTooBigToMove(const BamTools::BamAlignment & read);
    static bool iSizeTooBigToMove(const BamTools::BamAlignment & read, int maxInsertSizeForMovingReadPairs);
    
private:
    void purgeUnmodifiedMates();
    bool pairedReadIsMovable(BamTools::BamAlignment * read);
    
public:
    void close();
};

#endif
