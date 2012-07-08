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

#include "mark_duplicates.h"
#include <api/BamMultiReader.h>
#include <api/BamWriter.h>

#include <algorithm>
using namespace BamTools;
using namespace std;

MarkDuplicates::MarkDuplicates()
: numDuplicateIndices(0)
, nextLibraryId(1)
, removeDuplicates(false)
{
    char filename[48];
    sprintf(filename, "/tmp/dedup_%8x.bam",  (uint32_t)(0xffffffff & (uint64_t)this));

    setBufferFileName(string(filename));
}

string MarkDuplicates::getBufferFileName()
{
    return bufferFilename;
}

/////////////////
// From Samtools' SAMRecord.java:

int MarkDuplicates::getReferenceLength(const BamAlignment &rec) {
    int length = 0;
    for(vector<CigarOp>::const_iterator i = rec.getCigarData().begin(); i != rec.getCigarData().end(); i++) {
        switch (i->Type) {
            case 'M':
            case 'D':
            case 'N':
            case '=':
            case 'X':
                length += i->Length;
            default:
                break;
        }
    }
    return length;
}

/**
 * @return 1-based inclusive leftmost position of the clippped sequence, or 0 if there is no position.
 */
int MarkDuplicates::getAlignmentStart(const BamAlignment & rec) {
    return rec.getPosition();
}

/**
 * @return 1-based inclusive rightmost position of the clippped sequence, or 0 read if unmapped.
 */
int MarkDuplicates::getAlignmentEnd(const BamAlignment & rec) {
    if (!rec.IsMapped() ) {
        return -1;
    } else {
        return getAlignmentStart(rec) + getReferenceLength(rec)-1;
    }
}

/**
 * @return the alignment start (1-based, inclusive) adjusted for clipped bases.  For example if the read
 * has an alignment start of 100 but the first 4 bases were clipped (hard or soft clipped)
 * then this method will return 96.
 *
 * Invalid to call on an unmapped read.
 */
int MarkDuplicates::getUnclippedStart(const BamAlignment & rec) {
    int pos = getAlignmentStart(rec);
    
    for (vector<CigarOp>::const_iterator op = rec.getCigarData().begin(); op != rec.getCigarData().end(); op++ ) {
        if (op->Type == 'S' || op->Type == 'H') {
            pos -= op->Length;
        }
        else {
            break;
        }
    }
    
    return pos;
}

/**
 * @return the alignment end (1-based, inclusive) adjusted for clipped bases.  For example if the read
 * has an alignment end of 100 but the last 7 bases were clipped (hard or soft clipped)
 * then this method will return 107.
 *
 * Invalid to call on an unmapped read.
 */
int MarkDuplicates::getUnclippedEnd(const BamAlignment & rec) {
    int pos = getAlignmentEnd(rec);
    
    for (int i=rec.getCigarData().size() - 1; i>=0; --i) {
        const CigarOp & op = rec.getCigarData()[i];
        
        if (op.Type == 'S' || op.Type =='H') {
            pos += op.Length;
        }
        else {
            break;
        }
    }
    
    return pos;               
}

// end Samtools
/////////////////

/** Calculates a score for the read which is the sum of scores over Q20. */
short MarkDuplicates::getScore(const BamAlignment & rec) {
    short score = 0;
    for (int i = 0; i < rec.getQualities().size(); i++) {
        uint8_t b = rec.getQualities()[i]-33;   //33 comes from the conversion in BamAlignment
        if (b >= 15) score += b;
    }
    
    return score;
}

/** Builds a read ends object that represents a single read. */
ReadEnds * MarkDuplicates::buildReadEnds(SamHeader & header, long index, const BamAlignment & rec) {
    ReadEnds * ends = new ReadEnds();
    ends->read1Sequence    = rec.getRefID();
    ends->read1Coordinate  = rec.IsReverseStrand() ? getUnclippedEnd(rec) : getUnclippedStart(rec);
    ends->orientation = rec.IsReverseStrand() ? RE_R : RE_F;
    ends->read1IndexInFile = index;
    ends->score = getScore(rec);
    
    // Doing this lets the ends object know that it's part of a pair
    if (rec.IsPaired() && rec.IsMateMapped()) {
        ends->read2Sequence = rec.getMateRefID();
    }
    
    // Fill in the library ID
    ends->libraryId = getLibraryId(header, rec);
    
    return ends;
}

/**
 * Returns a single uint8_t that encodes the orientation of the two reads in a pair.
 */
readends_orientation_t MarkDuplicates::getOrientationByte(bool read1NegativeStrand, bool read2NegativeStrand) {
    if (read1NegativeStrand) {
        if (read2NegativeStrand)  return RE_RR;
        else return RE_RF;
    }
    else {
        if (read2NegativeStrand)  return RE_FR;
        else return RE_FF;
    }
}

/**
 * Goes through all the records in a file and generates a set of ReadEnds objects that
 * hold the necessary information (reference sequence, 5' read coordinate) to do
 * duplication, caching to disk as necssary to sort them.
 */
void MarkDuplicates::buildSortedReadEndLists() {
    
    ReadEndsMap tmp;
    long index = 0;

    SamHeader header = source->getHeader();
    
    BamWriter writer;

    writer.SetCompressionMode(BamWriter::Uncompressed);
    writer.Open(getBufferFileName(), header, source->getReferences());
    
    while (true) {
        BamAlignment * prec = getInputAlignment();
        if(!prec)
            break;

        BamAlignment & rec = *prec;

        if (!rec.IsMapped() || rec.getRefID() == -1) {
            // When we hit the unmapped reads or reads with no coordinate, just write them.
        }
        else if (rec.IsPrimaryAlignment()){
            ReadEnds * fragmentEnd = buildReadEnds(header, index, rec);
            fragSort.push_back(fragmentEnd);
            
            if (rec.IsPaired() && rec.IsMateMapped()) {
                string read_group;
                rec.GetTag("RG", read_group);
                string key = read_group + ":" + rec.getName();
                ReadEnds * pairedEnds = tmp.remove(rec.getRefID(), key);
                
                // See if we've already seen the first end or not
                if (pairedEnds == NULL) {
                    pairedEnds = buildReadEnds(header, index, rec);
                    tmp.put(pairedEnds->read1Sequence, key, pairedEnds);
                }
                else {
                    int sequence = fragmentEnd->read1Sequence;
                    int coordinate = fragmentEnd->read1Coordinate;
                    
                    // If the second read is actually later, just add the second read data, else flip the reads
                    if (sequence > pairedEnds->read1Sequence ||
                        (sequence == pairedEnds->read1Sequence && coordinate >= pairedEnds->read1Coordinate)) {
                        pairedEnds->read2Sequence    = sequence;
                        pairedEnds->read2Coordinate  = coordinate;
                        pairedEnds->read2IndexInFile = index;
                        pairedEnds->orientation = getOrientationByte(pairedEnds->orientation == RE_R, rec.IsReverseStrand());
                    }
                    else {
                        pairedEnds->read2Sequence    = pairedEnds->read1Sequence;
                        pairedEnds->read2Coordinate  = pairedEnds->read1Coordinate;
                        pairedEnds->read2IndexInFile = pairedEnds->read1IndexInFile;
                        pairedEnds->read1Sequence    = sequence;
                        pairedEnds->read1Coordinate  = coordinate;
                        pairedEnds->read1IndexInFile = index;
                        pairedEnds->orientation = getOrientationByte(rec.IsReverseStrand(), pairedEnds->orientation == RE_R);
                    }
                    
                    pairedEnds->score += getScore(rec);
                    pairSort.push_back(pairedEnds);
                }
            }
        }
        
        // Print out some stats every 1m reads
        if (verbose && ++index % 100000 == 0) {
            cerr << "\rRead " << index << " records. Tracking " << tmp.size() << " as yet unmatched pairs. Last sequence index: " << rec.getPosition() << std::flush;
        }
        
        writer.SaveAlignment(rec);
        delete prec;
    }

    writer.Close();
    
    if(verbose)
        cerr << "Read " << index << " records. " << tmp.size() << " pairs never matched." << endl << "Sorting pairs..." << flush;
    if(nothreads)
        sort(pairSort.begin(), pairSort.end(), compareReadEnds());
    else
        ogeSortMt(pairSort.begin(), pairSort.end(), compareReadEnds());

    if(verbose) cerr << "fragments..." << flush;
    if(nothreads)
        sort(fragSort.begin(), fragSort.end(), compareReadEnds());
    else
        ogeSortMt(fragSort.begin(), fragSort.end(), compareReadEnds());
    cerr << "done." << endl;
    
    vector<ReadEnds *>contents = tmp.allReadEnds();
    
    // delete unmatched read ends
    for(vector<ReadEnds *>::const_iterator i = contents.begin(); i != contents.end(); i++)
        delete *i;
}

/** Get the library ID for the given SAM record. */
short MarkDuplicates::getLibraryId(SamHeader & header, const BamAlignment & rec) {
    string library = getLibraryName(header, rec);
    
    short libraryId;
    
    if(!libraryIds.count(library)) {
        libraryId = nextLibraryId++;
        libraryIds[library] = libraryId;
    } else
        libraryId = libraryIds[library];
    
    return libraryId;
}

/**
 * Gets the library name from the header for the record. If the RG tag is not present on
 * the record, or the library isn't denoted on the read group, a constant string is
 * returned.
 */
string MarkDuplicates::getLibraryName(SamHeader & header, const BamAlignment & rec) {     
    
    string read_group;
    static const string RG("RG");
    static const string unknown_library("Unknown Library");
    rec.GetTag(RG, read_group);
    
    if (read_group.size() > 0 && header.ReadGroups.Contains(read_group)) {
        SamReadGroupDictionary & d = header.ReadGroups;
        const SamReadGroup & rg = d[read_group];
        
        if(rg.HasLibrary()) {
            return rg.Library;
        }
    }
    
    return unknown_library;
}

/**
 * Goes through the accumulated ReadEnds objects and determines which of them are
 * to be marked as duplicates.
 *
 * @return an array with an ordered list of indexes into the source file
 */
void MarkDuplicates::generateDuplicateIndexes() {
    
    ReadEnds * firstOfNextChunk = NULL;
    vector<ReadEnds *> nextChunk;
    nextChunk.reserve(200);
    
    // First just do the pairs
    if(verbose)
        cerr << "Finding duplicate pairs..." << flush;
    
    for (int i = 0;i < pairSort.size(); i++) {
        ReadEnds * next = pairSort[i];
        if (firstOfNextChunk == NULL) {
            firstOfNextChunk = next;
            nextChunk.push_back(firstOfNextChunk);
        }
        else if (areComparableForDuplicates(*firstOfNextChunk, *next, true)) {
            nextChunk.push_back(next);
        }
        else {
            if (nextChunk.size() > 1) {
                markDuplicatePairs(nextChunk);
            }
            
            nextChunk.clear();
            nextChunk.push_back(next);
            firstOfNextChunk = next;
        }
    }
    markDuplicatePairs(nextChunk);
    
    for (int i = 0;i < pairSort.size(); i++) {
        delete pairSort[i];
        pairSort[i] = NULL;
    }
    pairSort.clear();
    
    // Now deal with the fragments
    if(verbose)
        cerr << "duplicate fragments..." << flush;
    
    bool containsPairs = false;
    bool containsFrags = false;
    
    for (int i = 0;i < fragSort.size(); i++) {
        ReadEnds * next = fragSort[i];
        if (firstOfNextChunk != NULL && areComparableForDuplicates(*firstOfNextChunk, *next, false)) {
            nextChunk.push_back(next);
            containsPairs = containsPairs || next->isPaired();
            containsFrags = containsFrags || !next->isPaired();
        }
        else {
            if (nextChunk.size() > 1 && containsFrags) {
                markDuplicateFragments(nextChunk, containsPairs);
            }
            
            nextChunk.clear();
            nextChunk.push_back(next);
            firstOfNextChunk = next;
            containsPairs = next->isPaired();
            containsFrags = !next->isPaired();
        }
    }
    markDuplicateFragments(nextChunk, containsPairs);
    
    for (int i = 0;i < fragSort.size(); i++) {
        delete fragSort[i];
        fragSort[i] = NULL;
    }
    fragSort.clear();
    
    if(verbose)
        cerr << "done." << endl << "Sorting list of duplicate records." << endl;
}

bool MarkDuplicates::areComparableForDuplicates(const ReadEnds & lhs, const ReadEnds & rhs, bool compareRead2) {
    bool retval = (lhs.libraryId  == rhs.libraryId) &&
    (lhs.read1Sequence   == rhs.read1Sequence) &&
    (lhs.read1Coordinate == rhs.read1Coordinate) &&
    (lhs.orientation     == rhs.orientation);
    
    if (retval && compareRead2) {
        retval = (lhs.read2Sequence   == rhs.read2Sequence) &&
        (lhs.read2Coordinate == rhs.read2Coordinate);
    }
    
    return retval;
}

/**
 * Main work method.  Reads the BAM file once and collects sorted information about
 * the 5' ends of both ends of each read (or just one end in the case of pairs).
 * Then makes a pass through those determining duplicates before re-reading the
 * input file and writing it out with duplication flags set correctly.
 */
int MarkDuplicates::runInternal() {
    
    ogeNameThread("am_MarkDuplicates");

    if(verbose)
        cerr << "Reading input file and constructing read end information." << endl;
    buildSortedReadEndLists();
    
    generateDuplicateIndexes();
    
    if(verbose)
        cerr << "Marking " << numDuplicateIndices << " records as duplicates." << endl;
    
    BamReader in;

    in.Open(getBufferFileName());

    // Now copy over the file while marking all the necessary indexes as duplicates
    long recordInFileIndex = 0;
    
    long written = 0;
    while (true) {
        BamAlignment * prec = in.GetNextAlignment();
        if(!prec)
            break;
        
        if (prec->IsPrimaryAlignment()) {
            if (duplicateIndexes.count(recordInFileIndex) == 1)
                prec->SetIsDuplicate(true);
            else
                prec->SetIsDuplicate(false);
        }
        recordInFileIndex++;
        
        if (removeDuplicates && prec->IsDuplicate()) {
            // do nothing
        }
        else {
            putOutputAlignment(prec);
            if (verbose && read_count && ++written % 100000 == 0) {
                cerr << "\rWritten " << written << " records (" << written * 100 / read_count <<"%)." << std::flush;
            }
        }
    }

    if (verbose && read_count)
        cerr << "\rWritten " << written << " records (" << written * 100 / read_count <<"%)." << endl;

    in.Close();

    remove(getBufferFileName().c_str());
    
    return 0;
}

void MarkDuplicates::addIndexAsDuplicate(long bamIndex) {
    duplicateIndexes.insert(bamIndex);
    ++numDuplicateIndices;
}

/**
 * Takes a list of ReadEnds objects and removes from it all objects that should
 * not be marked as duplicates.
 *
 * @param list
 */
void MarkDuplicates::markDuplicatePairs(const vector<ReadEnds *>& list) {
    short maxScore = 0;
    ReadEnds * best = NULL;
    
    for (int i = 0; i < list.size(); i++) {
        ReadEnds * end = list[i];
        if (end->score > maxScore || best == NULL) {
            maxScore = end->score;
            best = end;
        }
    }
    
    for (int i = 0; i < list.size(); i++) {
        ReadEnds * end = list[i];
        if (end != best) {
            addIndexAsDuplicate(end->read1IndexInFile);
            addIndexAsDuplicate(end->read2IndexInFile);
        }
    }
}

/**
 * Takes a list of ReadEnds objects and removes from it all objects that should
 * not be marked as duplicates.
 *
 * @param list
 */
void MarkDuplicates::markDuplicateFragments(const vector<ReadEnds *>& list, bool containsPairs) {
    if (containsPairs) {
        for (int i = 0; i < list.size(); i++) {
            ReadEnds * end = list[i];
            if (!end->isPaired()) 
                addIndexAsDuplicate(end->read1IndexInFile);
        }
    }
    else {
        short maxScore = 0;
        ReadEnds * best = NULL;
        for (int i = 0; i < list.size(); i++) {
            ReadEnds * end = list[i];
            if (end->score > maxScore || best == NULL) {
                maxScore = end->score;
                best = end;
            }
        }
        
        for (int i = 0; i < list.size(); i++) {
            ReadEnds * end = list[i];
            if (end != best)
                addIndexAsDuplicate(end->read1IndexInFile);
        }
    }
}

/** Comparator for ReadEnds that orders by read1 position then pair orientation then read2 position. */
class ReadEndsComparator {
    int compare(ReadEnds lhs, ReadEnds rhs) {
        int retval = lhs.libraryId - rhs.libraryId;
        if (retval == 0) retval = lhs.read1Sequence - rhs.read1Sequence;
        if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
        if (retval == 0) retval = lhs.orientation - rhs.orientation;
        if (retval == 0) retval = lhs.read2Sequence   - rhs.read2Sequence;
        if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
        if (retval == 0) retval = (int) (lhs.read1IndexInFile - rhs.read1IndexInFile);
        if (retval == 0) retval = (int) (lhs.read2IndexInFile - rhs.read2IndexInFile);
        
        return retval;
    }
};
