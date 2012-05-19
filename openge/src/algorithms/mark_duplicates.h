#ifndef DEDUP_H
#define DEDUP_H

#include "algorithm_module.h"
#include "api/BamAlignment.h"
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
    bool verbose;
    MarkDuplicates();

    std::string getBufferFileName();
    void setBufferFileName(std::string filename) { bufferFilename = filename; }

protected:
    /////////////////
    // From Samtools' SAMRecord.java:
    int getReferenceLength(const BamTools::BamAlignment &rec);
    int getAlignmentStart(const BamTools::BamAlignment & rec); 
    int getAlignmentEnd(const BamTools::BamAlignment & rec);
    int getUnclippedStart(const BamTools::BamAlignment & rec);
    int getUnclippedEnd(const BamTools::BamAlignment & rec);

    ////////////////
    // From Picard MarkDuplicates.java
    short getScore(const BamTools::BamAlignment & rec);
    ReadEnds * buildReadEnds(BamTools::SamHeader & header, long index, const BamTools::BamAlignment & rec);
    readends_orientation_t getOrientationByte(bool read1NegativeStrand, bool read2NegativeStrand);
    void buildSortedReadEndLists();
    short getLibraryId(BamTools::SamHeader & header, const BamTools::BamAlignment & rec);
    std::string getLibraryName(BamTools::SamHeader & header, const BamTools::BamAlignment & rec);
    void generateDuplicateIndexes();
    bool areComparableForDuplicates(const ReadEnds & lhs, const ReadEnds & rhs, bool compareRead2);
    void addIndexAsDuplicate(long bamIndex);
    void markDuplicatePairs(const std::vector<ReadEnds *>& list);
    void markDuplicateFragments(const std::vector<ReadEnds *>& list, bool containsPairs);

    int runInternal();
};

#endif