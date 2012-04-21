//
//  command_dedup.cpp
//  OpenGE
//
//  Created by Lee Baker on 4/14/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#include "commands.h"
#include <iostream>

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
using namespace BamTools;
namespace po = boost::program_options;
using namespace std;

typedef enum
{
    RE_NONE, RE_F, RE_R, RE_FF, RE_RR, RE_FR, RE_RF
} readends_orientation_t;

class ReadEnds{
public:

    short libraryId;
    short score;
    readends_orientation_t orientation;
    int read1Sequence;
    int read1Coordinate;
    long read1IndexInFile;
    int read2Sequence;
    int read2Coordinate;
    long read2IndexInFile;

    ReadEnds()
    : libraryId(-1)
    , score(-1)
    , orientation(RE_NONE)
    , read1Sequence(-1)
    , read1Coordinate(-1)
    , read1IndexInFile(-1)
    , read2Sequence(-1)
    , read2Coordinate(-1)
    , read2IndexInFile(-1)
    {}
    
    bool isPaired() { return read2Sequence != -1; }
    
    static int compare(const ReadEnds & lhs, const ReadEnds & rhs) {
        int retval = 0;
        if (retval == 0) retval = lhs.libraryId - rhs.libraryId;
        if (retval == 0) retval = lhs.read1Sequence - rhs.read1Sequence;
        if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
        if (retval == 0) retval = lhs.orientation - rhs.orientation;
        if (retval == 0) retval = lhs.read2Sequence   - rhs.read2Sequence;
        if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
        if (retval == 0) retval = (int) (lhs.read1IndexInFile - rhs.read1IndexInFile);
        if (retval == 0) retval = (int) (lhs.read2IndexInFile - rhs.read2IndexInFile);
        
        return retval;
    }
    
    bool operator<(const ReadEnds & a) const
    {
        return 0 > compare(*this, a);
    }
};
ostream& operator<< (ostream& out, const ReadEnds & re )
{
    out << "ReadEnds (LID " << re.libraryId << ")" << endl;
    out << " Seq: " << re.read1Sequence << "/" << re.read2Sequence << endl;
    out << " Coord: " << re.read1Coordinate << "/" << re.read2Coordinate << endl;
    out << " Orientation: " << re.orientation << endl;
    out << " Score: " << re.score << endl;
    
    return out;
}

struct compareReadEnds {
    bool operator ()(const ReadEnds *lhs, const ReadEnds *rhs) { return *lhs < *rhs; }
    bool operator ()(const ReadEnds & lhs, const ReadEnds & rhs) { return lhs < rhs; }
};

class ReadEndsMap
{
protected:
    map<string, ReadEnds *> m;
public:
    void put(int index, string key, ReadEnds * val) {
        m[key] = val;
    }

    ReadEnds * remove(int index, string key)
    {
        ReadEnds * ret = m[key];
        m.erase(key);
        return ret;
    }
    /*
    ~ReadEndsMap()
    {
        for (map<pair<int, string>, ReadEnds *>::iterator iter = m.begin(); iter != m.end(); ++iter)
            cerr << "Missing: " << iter->first.first << "::" << iter->first.second << endl;
    }*/
    
    size_t size() { return m.size();};
    
};


class MarkDuplicates {
protected:
    vector<ReadEnds *> pairSort;
    vector<ReadEnds *> fragSort;
    set<int> duplicateIndexes;
    int numDuplicateIndices;
    
    bool removeDuplicates;
    
    map<string,short> libraryIds;
    short nextLibraryId;
    
public:
    string filename_out;
    vector<string> filenames_in;
    MarkDuplicates()
    : numDuplicateIndices(0)
    , removeDuplicates(false)
    , nextLibraryId(1)
    {}
    /////////////////
    // From Samtools' SAMRecord.java:
#if 1
    
    int getReferenceLength(const BamAlignment * rec) {
        int length = 0;
        for(vector<CigarOp>::const_iterator i = rec->CigarData.begin(); i != rec->CigarData.end(); i++) {
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
    int getAlignmentStart(const BamAlignment * rec) {
        return rec->Position;
    }
    
    /**
     * @return 1-based inclusive rightmost position of the clippped sequence, or 0 read if unmapped.
     */
    int getAlignmentEnd(const BamAlignment * rec) {
        if (!rec->IsMapped() ) {
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
    int getUnclippedStart(const BamAlignment * rec) {
        int pos = getAlignmentStart(rec);
        
        for (vector<CigarOp>::const_iterator op = rec->CigarData.begin(); op != rec->CigarData.end(); op++ ) {
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
    int getUnclippedEnd(const BamAlignment * rec) {
        int pos = getAlignmentEnd(rec);
        
        for (int i=rec->CigarData.size() - 1; i>=0; --i) {
            const CigarOp & op = rec->CigarData[i];
            
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
#endif

    /** Calculates a score for the read which is the sum of scores over Q20. */
    short getScore(const BamAlignment * rec) {
        short score = 0;
        for (int i = 0; i < rec->Qualities.size(); i++) {
            uint8_t b = rec->Qualities[i]-33;   //33 comes from the conversion in BamAlignment
            if (b >= 15) score += b;
        }
        
        return score;
    }
    
    /** Builds a read ends object that represents a single read. */
    ReadEnds * buildReadEnds(SamHeader & header, long index, BamAlignment * rec) {
        ReadEnds * ends = new ReadEnds();
        ends->read1Sequence    = rec->RefID;
        ends->read1Coordinate  = rec->IsReverseStrand() ? getUnclippedEnd(rec) : getUnclippedStart(rec);
        ends->orientation = rec->IsReverseStrand() ? RE_R : RE_F;
        ends->read1IndexInFile = index;
        ends->score = getScore(rec);
        
        // Doing this lets the ends object know that it's part of a pair
        if (rec->IsPaired() && rec->IsMateMapped()) {
            ends->read2Sequence = rec->MateRefID;
        }
        
        // Fill in the library ID
        ends->libraryId = getLibraryId(header, rec);
        
        return ends;
    }
    
    /**
     * Returns a single uint8_t that encodes the orientation of the two reads in a pair.
     */
    readends_orientation_t getOrientationByte(bool read1NegativeStrand, bool read2NegativeStrand) {
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
    void buildSortedReadEndLists() {

        ReadEndsMap tmp;
        long index = 0;

        BamReader reader;
        reader.Open(filenames_in.front());
        SamHeader header = reader.GetHeader();
        
        while (true) {
            BamAlignment * rec = reader.GetNextAlignment();
            if(!rec) break;
            
            if (!rec->IsMapped()) {
                if (rec->RefID == -1) {
                    // When we hit the unmapped reads with no coordinate, no reason to continue.
                    break;
                }
                // If this read is unmapped but sorted with the mapped reads, just skip it.
            }
            else if (rec->IsPrimaryAlignment()){
                ReadEnds * fragmentEnd = buildReadEnds(header, index, rec);
                fragSort.push_back(fragmentEnd);
                
                if (rec->IsPaired() && rec->IsMateMapped()) {
                    string read_group;
                    rec->GetTag("RG", read_group);
                    string key = read_group + ":" + rec->Name;
                    ReadEnds * pairedEnds = tmp.remove(rec->RefID, key);
                    
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
                            pairedEnds->orientation = getOrientationByte(pairedEnds->orientation == RE_R, rec->IsReverseStrand());
                        }
                        else {
                            pairedEnds->read2Sequence    = pairedEnds->read1Sequence;
                            pairedEnds->read2Coordinate  = pairedEnds->read1Coordinate;
                            pairedEnds->read2IndexInFile = pairedEnds->read1IndexInFile;
                            pairedEnds->read1Sequence    = sequence;
                            pairedEnds->read1Coordinate  = coordinate;
                            pairedEnds->read1IndexInFile = index;
                            pairedEnds->orientation = getOrientationByte(rec->IsReverseStrand(), pairedEnds->orientation == RE_R);
                        }

                        pairedEnds->score += getScore(rec);
                        pairSort.push_back(pairedEnds);
                    }
                }
            }
            
            // Print out some stats every 1m reads
            if (++index % 100000 == 0) {
                cerr << "Read " << index << " records. Tracking " << tmp.size() << " as yet unmatched pairs. Last sequence index: " << rec->Position << "\r" << endl;
            }
            
            delete rec;
        }
        
        reader.Close();
        
        cerr << "Read " << index << " records. " << tmp.size() << " pairs never matched." << endl;
        sort(pairSort.begin(), pairSort.end(), compareReadEnds());
        sort(fragSort.begin(), fragSort.end(), compareReadEnds());
    }
    
    /** Get the library ID for the given SAM record. */
    short getLibraryId(SamHeader & header, const BamAlignment * rec) {
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
    string getLibraryName(SamHeader & header, const BamAlignment * rec) {     

        string read_group;
        rec->GetTag("RG", read_group);

        if (read_group.size() > 0 && header.ReadGroups.Contains(read_group)) {
            SamReadGroupDictionary & d = header.ReadGroups;
            const SamReadGroup & rg = d[read_group];
            
            if(rg.HasLibrary()) {
                return rg.Library;
            }
        }

        return "Unknown Library";
    }
    
    /**
     * Goes through the accumulated ReadEnds objects and determines which of them are
     * to be marked as duplicates.
     *
     * @return an array with an ordered list of indexes into the source file
     */
    void generateDuplicateIndexes() {

        ReadEnds * firstOfNextChunk = NULL;
        vector<ReadEnds *> nextChunk;
        nextChunk.reserve(200);

        // First just do the pairs
        cerr << "Traversing read pair information and detecting duplicates." << endl;
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
        cerr << "Traversing fragment information and detecting duplicates." << endl;
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

        cerr << "Sorting list of duplicate records." << endl;
    }
    
    bool areComparableForDuplicates(const ReadEnds & lhs, const ReadEnds & rhs, bool compareRead2) {
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
    int doWork() {
        
        cerr << "Reading input file and constructing read end information." << endl;
        buildSortedReadEndLists();
        
        generateDuplicateIndexes();
        
        cerr << "Marking " << numDuplicateIndices << " records as duplicates." << endl;
        
        BamMultiReader in;
        in.Open(filenames_in);
        
        BamWriter out;
        out.Open(filename_out,in.GetHeader(),in.GetReferenceData());
        
        // Now copy over the file while marking all the necessary indexes as duplicates
        long recordInFileIndex = 0;
        
        long written = 0;
        while (true) {
            BamAlignment rec;
            if(!in.GetNextAlignmentCore(rec))
                break;
            
            if (rec.IsPrimaryAlignment()) {
                if (duplicateIndexes.count(recordInFileIndex) == 1)
                    rec.SetIsDuplicate(true);
                else
                    rec.SetIsDuplicate(false);
            }
            recordInFileIndex++;
            
            if (removeDuplicates && rec.IsDuplicate()) {
                // do nothing
            }
            else {
                out.SaveAlignment(rec);
                if (++written % 100000 == 0) {
                    cerr << "Written " << written << " records.\r" << endl;
                }
            }
        }
        
        in.Close();
        out.Close();
        
        return 0;
    }
    
    void addIndexAsDuplicate(long bamIndex) {
        duplicateIndexes.insert(bamIndex);
        ++numDuplicateIndices;
    }
    
    /**
     * Takes a list of ReadEnds objects and removes from it all objects that should
     * not be marked as duplicates.
     *
     * @param list
     */
    void markDuplicatePairs(const vector<ReadEnds *>& list) {
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
    void markDuplicateFragments(const vector<ReadEnds *>& list, bool containsPairs) {
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
};

void DedupCommand::getOptions()
{
    options.add_options()
    ("out,o", po::value<string>()->default_value("stdout"), "Output filename. Omit for stdout.")
    ;
}

int DedupCommand::runCommand()
{
    MarkDuplicates md;
    md.filename_out = vm["out"].as<string>();
    md.filenames_in = input_filenames;
    
    if(input_filenames.size() == 1 && input_filenames[0] == string("stdin"))
    {
        cerr << "The current deduplication algorithm does not support reading from" << endl << "stdin. Please provide a filename to read from." << endl;
        return -1;
    }
    
    if(input_filenames.size() != 1)
    {
        cerr << "One input file is required. You supplied " << input_filenames.size() << endl;
        return -1;
    }

    
    md.doWork();
    
    return 0;
}
