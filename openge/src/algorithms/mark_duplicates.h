#ifndef OGE_ALGO_DEDUP_H
#define OGE_ALGO_DEDUP_H

/*********************************************************************
 *
 * mark_duplicates.cpp:  Mark duplicate reads.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 20 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include "algorithm_module.h"
#include "../util/picard_structures.h"

#include <map>
#include <string>
#include <vector>

class MarkDuplicates : public AlgorithmModule
{
protected:
    std::vector<ReadEnds *> pairSort;
    std::vector<ReadEnds *> fragSort;
    std::set<int> duplicateIndexes;
    int numDuplicateIndices;
    
    std::map<std::string,short> libraryIds;
    short nextLibraryId;
    
    std::string bufferFilename;
    
public:
    bool removeDuplicates;
    MarkDuplicates();

    std::string getBufferFileName();
    void setBufferFileName(std::string filename) { bufferFilename = filename; }

protected:
    /////////////////
    // From Samtools' SAMRecord.java:
    int getReferenceLength(const OGERead &rec);
    int getAlignmentStart(const OGERead & rec); 
    int getAlignmentEnd(const OGERead & rec);
    int getUnclippedStart(const OGERead & rec);
    int getUnclippedEnd(const OGERead & rec);

    ////////////////
    // From Picard MarkDuplicates.java
    short getScore(const OGERead & rec);
    ReadEnds * buildReadEnds(BamTools::SamHeader & header, long index, const OGERead & rec);
    readends_orientation_t getOrientationByte(bool read1NegativeStrand, bool read2NegativeStrand);
    void buildSortedReadEndLists();
    short getLibraryId(BamTools::SamHeader & header, const OGERead & rec);
    std::string getLibraryName(BamTools::SamHeader & header, const OGERead & rec);
    void generateDuplicateIndexes();
    bool areComparableForDuplicates(const ReadEnds & lhs, const ReadEnds & rhs, bool compareRead2);
    void addIndexAsDuplicate(long bamIndex);
    void markDuplicatePairs(const std::vector<ReadEnds *>& list);
    void markDuplicateFragments(const std::vector<ReadEnds *>& list, bool containsPairs);

    int runInternal();
};

#endif