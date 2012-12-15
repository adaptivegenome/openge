/*********************************************************************
 *
 * local_realignment.cpp:  Local realignment algorithm.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 31 May 2012
 *
 *********************************************************************
 *
 * This file has been ported from GATK's implementation in Java, and
 * is distributed under the license listed below.
 *
 *********************************************************************/
/*
 * Copyright (c) 2010 The Broad Institute. 
 * Ported to C++ by Lee C. Baker, Virginia Bioinformatics Institute
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

#include "local_realignment.h"
#include "../util/bamtools/Sort.h"

#include <cassert>
#include <ctime>
#include <sys/time.h>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <cstring>
#include <string>
#include <sstream>
#include <cstdlib>  //for rand()
#include <cassert>
#include <algorithm>
using namespace std;

/**
 * Performs local realignment of reads based on misalignments due to the presence of ;.
 *
 * <p>
 * The local realignment tool is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases
 * is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion
 * or deletion (indels) in the individual's genome with respect to the reference genome.  Such alignment artifacts result in many bases mismatching
 * the reference near the misalignment, which are easily mistaken as SNPs.  Moreover, since read mapping algorithms operate on each read independently,
 * it is impossible to place reads on the reference genome such at mismatches are minimized across all reads.  Consequently, even when some reads are
 * correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel,
 * also requiring realignment.  Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus
 * indel suitable for standard variant discovery approaches.  Unlike most mappers, this walker uses the full alignment context to determine whether an
 * appropriate alternate reference (i.e. indel) exists.  Following local realignment, the GATK tool Unified Genotyper can be used to sensitively and
 * specifically identify indels.
 * <p>
 *     <ol>There are 2 steps to the realignment process:
 *     <li>Determining (small) suspicious intervals which are likely in need of realignment (see the RealignerTargetCreator tool)</li>
 *     <li>Running the realigner over those intervals (LocalRealignment)</li>
 *     </ol>
 *     <p>
 * An important note: the input bam(s), reference, and known indel file(s) should be the same ones used for the RealignerTargetCreator step.
 * <p>
 * Another important note: because reads produced from the 454 technology inherently contain false indels, the realigner will not currently work with them
 * (or with reads from similar technologies).
 *
 * <h2>Input</h2>
 * <p>
 * One or more aligned BAM files and optionally one or more lists of known indels.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A realigned version of your input BAM file(s).
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -I input.bam \
 *   -R ref.fasta \
 *   -T LocalRealignment \
 *   -targetIntervals intervalListFromRTC.intervals \
 *   -o realignedBam.bam \
 *   [--known /path/to/indels.vcf] \
 *   [-compress 0]    (this argument recommended to speed up the process *if* this is only a temporary file; otherwise, use the default value)
 * </pre>
 *
 * @author ebanks
 */

static const string ORIGINAL_CIGAR_TAG("OC");
static const string ORIGINAL_POSITION_TAG("OP");
static const string PROGRAM_RECORD_NAME("OpenGE LocalRealignment");


const int LocalRealignment::MAX_QUAL = 99;
const double LocalRealignment::MISMATCH_COLUMN_CLEANED_FRACTION = 0.75;
const int LocalRealignment::REFERENCE_PADDING = 30;
    
#pragma mark IndelRealinger::AlignedRead


int LocalRealignment::AlignedRead::MAX_POS_MOVE_ALLOWED;
int LocalRealignment::AlignedRead::NO_ORIGINAL_ALIGNMENT_TAGS;

string  LocalRealignment::AlignedRead::getReadBases() {
    if ( readBases.size() == 0 )
        getUnclippedBases();
    return readBases;
}

string  LocalRealignment::AlignedRead::getBaseQualities() {
    if ( baseQuals.size() == 0 )
        getUnclippedBases();
    return baseQuals;
}

// pull out the bases that aren't clipped out
void LocalRealignment::AlignedRead::getUnclippedBases() {
    readBases.clear();
    readBases.reserve(getReadLength());
    baseQuals.clear();
    baseQuals.reserve(getReadLength());
    const string & actualReadBases = read->getQueryBases();
    const string & actualBaseQuals = read->getQualities();

    int fromIndex = 0, toIndex = 0;
    
    vector<CigarOp> cigar = read->getCigarData();
    for ( vector<CigarOp>::const_iterator ce = cigar.begin(); ce != cigar.end(); ce++ ) {
        uint32_t elementLength = ce->length;
        switch ( ce->type ) {
            case 'S':
                fromIndex += elementLength;
                break;
            case 'M':
            case 'I':
                readBases.append(actualReadBases,fromIndex, elementLength);
                baseQuals.append(actualBaseQuals, fromIndex, elementLength);
                fromIndex += elementLength;
                toIndex += elementLength;
            default:
                break;
        }
    }
}

// pull out the bases that aren't clipped out
vector<CigarOp> LocalRealignment::AlignedRead::reclipCigar(const vector<CigarOp> & cigar) const {
    return LocalRealignment::reclipCigar(cigar, read);
}

const vector<CigarOp> LocalRealignment::AlignedRead::getCigar() const {
    return (newCigar.size() > 0 ? newCigar : read->getCigarData());
}

// tentatively sets the new Cigar, but it needs to be confirmed later
void LocalRealignment::AlignedRead::setCigar(const vector<CigarOp> & cigar_in, bool fixClippedCigar) {

    bool reclip = false;
    if ( fixClippedCigar && getReadBases().length() < read->getLength() ) {
        reclip = true;
    }
    
    const vector<CigarOp> cigar = reclip ? reclipCigar(cigar_in) : cigar_in;

    // no change?
    if ( read->getCigarData() == cigar ) {
        clearCigar();
        return;
    }

    // no indel?
    bool has_i_or_d = false;
    for (vector<CigarOp>::const_iterator i = cigar.begin(); i != cigar.end(); i++)
        if(i->type == 'I' || i->type == 'D')
            has_i_or_d = true;

    if ( !has_i_or_d) {
        cerr << "Modifying a read with no associated indel; although this is possible, it is highly unlikely.  Perhaps this region should be double-checked: " << read->getName() << " near " << (*sequences)[read->getRefID()].getName() << ":" << read->getPosition() << endl;
    }
    newCigar = cigar;
}

// tentatively sets the new start, but it needs to be confirmed later
void LocalRealignment::AlignedRead::setAlignmentStart(int start) {
    newStart = start;
}

int LocalRealignment::AlignedRead::getAlignmentStart() const{
    return (newStart != -1 ? newStart : read->getPosition());
}

int LocalRealignment::AlignedRead::getOriginalAlignmentStart() const {
    return read->getPosition();
}

// constizes the changes made.
// returns true if this record actually changes, false otherwise
bool LocalRealignment::AlignedRead::constizeUpdate() {
    // if we haven't made any changes, don't do anything
    if ( newCigar.size() == 0 ) {
        return false;
    }
    if ( newStart == -1 )
        newStart = read->getPosition();
    else if ( abs(newStart - read->getPosition()) > MAX_POS_MOVE_ALLOWED ) {
        cerr << "Attempting to realign read " << read->getName() << " at " << read->getPosition() << " more than " << MAX_POS_MOVE_ALLOWED << " bases to " << newStart << "." << endl;
        return false;
    }
    
    // annotate the record with the original cigar (and optionally the alignment start)
    if ( !NO_ORIGINAL_ALIGNMENT_TAGS ) {
        read->AddTag(ORIGINAL_CIGAR_TAG, "Z", read->cigarString());
        if ( newStart != read->getPosition() )
            read->AddTag(ORIGINAL_POSITION_TAG, "i", read->getPosition()+1);
    }
        
    read->setCigarData(newCigar);
    read->setPosition(newStart);
    
    return true;
}

void LocalRealignment::AlignedRead::setMismatchScoreToReference(int score) {
    mismatchScoreToReference = score;
}

int LocalRealignment::AlignedRead::getMismatchScoreToReference() const {
    return mismatchScoreToReference;
}

void LocalRealignment::AlignedRead::setAlignerMismatchScore(long score) {
    alignerMismatchScore = score;
}

long LocalRealignment::AlignedRead::getAlignerMismatchScore() const {
    return alignerMismatchScore;
}

#pragma mark IndelRealinger::ReadBin

// Return false if we can't process this read bin because the reads are not correctly overlapping.
// This can happen if e.g. there's a large known indel with no overlapping reads.
void LocalRealignment::ReadBin::add(OGERead * read) {
    
    GenomeLoc locForRead = loc_parser->createGenomeLoc(*read);
    if ( loc == NULL )
        loc = new GenomeLoc(locForRead);
    else if ( locForRead.getStop() > loc->getStop() ) {
        GenomeLoc * old_loc = loc;
        loc = new GenomeLoc(loc_parser->createGenomeLoc(loc->getContig(), loc->getStart(), locForRead.getStop()));
        delete old_loc;
    }

    reads.push_back(read);
}

string LocalRealignment::ReadBin::getReference(FastaReader & referenceReader) {
    // set up the reference if we haven't done so yet
    if ( reference.size() == 0 ) {
        // first, pad the reference to handle deletions in narrow windows (e.g. those with only 1 read)
        int padLeft = max(loc->getStart()-REFERENCE_PADDING, 0);    //LCB 1->0
        int padRight = min(loc->getStop()+REFERENCE_PADDING, (int)referenceReader.getSequenceDictionary()[loc->getContig()].getLength()-1);
        GenomeLoc * old_loc = loc;
        loc = new GenomeLoc(loc_parser->createGenomeLoc(loc->getContig(), padLeft, padRight));
        delete old_loc;
        reference = referenceReader.readSequence(loc->getContig(), loc->getStart(), loc->getStop() - loc->getStart() + 1);

        //transform to upper case. STL sometimes makes things harder than they should be.
        std::transform(reference.begin(), reference.end(), reference.begin(), (int (*)(int))std::toupper); 
    }
    
    return reference;
}

void LocalRealignment::ReadBin::clear() {
    reads.clear();
    reference.clear();
    delete loc;
    loc = NULL;
}

void LocalRealignment::initialize() {
    if ( LOD_THRESHOLD < 0.0 ) {
        cerr << "LOD threshold cannot be a negative number" << endl;
        exit(-1);
    }
    if ( MISMATCH_THRESHOLD <= 0.0 || MISMATCH_THRESHOLD > 1.0 ) {
        cerr << "Entropy threshold must be a fraction between 0 and 1" << endl;
        exit(-1);
    }

    referenceReader = new FastaReader();
    assert(reference_filename.size() > 0);
    if(!referenceReader->open(reference_filename)) {
        cerr << "Error opening reference file " << reference_filename << endl;
        exit(-1);
    }
    
    loc_parser = new GenomeLocParser(sequence_dictionary);

    // load intervals file
    {
        int line_counter = 0;
        ifstream intervals_file(intervals_filename.c_str());
        
        while(!intervals_file.eof())
        {
            string line;
            getline(intervals_file, line); 
            if(line.size() == 0)
                break;
            GenomeLoc g = loc_parser->parseGenomeLoc(line);
            
            intervalsFile.push_back(new GenomeLoc(g));
            line_counter++;
        }
        if(verbose)
            cerr << "Parsed " << line_counter << " intervals from " << intervals_filename << "." << endl;
    }
    
    if(intervalsFile.empty()) {
        cerr << "Error parsing intervals file. Aborting." << endl;
        exit(-1);
    }

    
    AlignedRead::NO_ORIGINAL_ALIGNMENT_TAGS = NO_ORIGINAL_ALIGNMENT_TAGS;
    AlignedRead::MAX_POS_MOVE_ALLOWED = MAX_POS_MOVE_ALLOWED;

    // set up the output writer
    manager = new ConstrainedMateFixingManager(this, MAX_ISIZE_FOR_MOVEMENT, MAX_POS_MOVE_ALLOWED, MAX_RECORDS_IN_MEMORY, loc_parser);
    
#ifdef LR_SUPPORT_ADDITIONAL_OUTPUT_FILES
    if ( write_out_indels )
        indelOutput.open(OUT_INDELS.c_str());
    if ( write_out_stats )
        statsOutput.open(OUT_STATS.c_str());
    if ( write_out_snps )
        snpsOutput.open(OUT_SNPS.c_str());
#endif
}

void LocalRealignment::emit(IntervalData & interval_data, OGERead * read) {
    
    // check to see whether the read was modified by looking at the temporary tag
    bool wasModified = interval_data.readsActuallyCleaned.count(read) > 0;
    manager->addRead(read, wasModified, true);
}

void LocalRealignment::emitReadLists(IntervalData & interval_data) {
    // pre-merge lists to sort them in preparation for constrained SAMFileWriter
    vector<OGERead *> to_insert = interval_data.readsToClean.getReads();
    interval_data.readsNotToClean.insert(interval_data.readsNotToClean.end(), to_insert.begin(), to_insert.end());
    std::stable_sort( interval_data.readsNotToClean.begin(), interval_data.readsNotToClean.end(), BamTools::Algorithms::Sort::ByPosition() );
    manager->addReads(interval_data.readsNotToClean, interval_data.readsActuallyCleaned);
}

LocalRealignment::Emittable::~Emittable() {}

class LocalRealignment::EmittableRead : public LocalRealignment::Emittable
{
    LocalRealignment::IntervalData * id;
    OGERead * read;
    LocalRealignment & lr;
public:
    EmittableRead(LocalRealignment & lr, IntervalData * id, OGERead * read)
    : id(id)
    , read(read)
    , lr(lr)
    {}
    virtual bool canEmit() { return true; }
    virtual void emit() {
        lr.emit(*id, read);
    }
    virtual ~EmittableRead() {
    }
};

class LocalRealignment::EmittableReadList : public LocalRealignment::Emittable
{
    LocalRealignment::IntervalData * id;
    LocalRealignment & lr;
public:
    EmittableReadList(LocalRealignment & lr, IntervalData * id)
    : id(id)
    , lr(lr)
    {}
    virtual bool canEmit() { return true; }
    virtual void emit() {
        lr.emitReadLists(*id);
    }
    virtual ~EmittableReadList() {
        delete id;
    }
};

class LocalRealignment::CleanAndEmitReadList : public LocalRealignment::Emittable
{
public:
    bool clean_done;
    LocalRealignment & lr;
    LocalRealignment::IntervalData * id;

    CleanAndEmitReadList(LocalRealignment & lr, IntervalData * id)
    : clean_done(false)
    , lr(lr)
    , id(id)
    { }
    virtual bool canEmit() { return clean_done; }
    virtual void emit() {
        lr.emitReadLists(*id);
    }
    virtual ~CleanAndEmitReadList() {
        delete id;
    }
};

class LocalRealignment::CleanJob : public ThreadJob
{
    CleanAndEmitReadList * c;
    bool self_delete;
public:
    CleanJob(CleanAndEmitReadList * c, bool self_delete = false)
    : c(c)
    , self_delete(self_delete)
    {}
    
    virtual bool deleteOnCompletion() { return self_delete; }

    virtual void runJob()
    {
        
        ogeNameThread("LRCleanJob");
        IntervalData * interval_data = c->id;
        if ( interval_data->readsToClean.size() > 0 ) {
            GenomeLoc earliestPossibleMove = c->lr.loc_parser->createGenomeLoc(*(interval_data->readsToClean.getReads()[0]));
            if ( c->lr.manager->canMoveReads(earliestPossibleMove) )
                c->lr.clean(*interval_data);
        }
        interval_data->knownIndelsToTry.clear();
        interval_data->indelRodsSeen.clear();
        
        c->clean_done = true;
        
        c->lr.flushEmitQueue();
    }
};

int LocalRealignment::map_func(OGERead * read, const ReadMetaDataTracker & metaDataTracker) {
    const int NO_ALIGNMENT_REFERENCE_INDEX = -1;
    if ( !loading_interval_data->current_interval_valid ) {
        pushToEmitQueue(new EmittableRead(*this, loading_interval_data, read));
        loading_interval_data = new IntervalData(NULL);
        loading_interval_data->readsToClean.initialize(loc_parser, &sequence_dictionary);
        return 0;
    }
    
    // edge case: when the last target interval abuts the end of the genome, we'll get one of the
    //   unmapped reads while the currentInterval still isn't null.  We need to trigger the cleaning
    //   at this point without trying to create a GenomeLoc.
    if ( read->getRefID() == NO_ALIGNMENT_REFERENCE_INDEX ) {
        CleanAndEmitReadList * read_list = new CleanAndEmitReadList(*this, loading_interval_data);
        pushToEmitQueue(read_list);
        if(nothreads)
            CleanJob(read_list).runJob();
        else
            ThreadPool::sharedPool()->addJob(new CleanJob(read_list));
        flushEmitQueue();

        do {
            ++interval_it;
        } while ( interval_it != intervalsFile.end() );
        current_interval = NULL;
        loading_interval_data = new IntervalData(NULL);
        loading_interval_data->readsToClean.initialize(loc_parser, &sequence_dictionary);
        
        sawReadInCurrentInterval = false;
        map_func(read, metaDataTracker);

        return 0;
    }
    GenomeLoc readLoc = loc_parser->createGenomeLoc(*read);
    // hack to get around unmapped reads having screwy locations
    if ( readLoc.getStop() == 0 )
        readLoc = loc_parser->createGenomeLoc(readLoc.getContig(), readLoc.getStart(), readLoc.getStart());
    
    if ( readLoc.isBefore(loading_interval_data->current_interval) ) {
        if ( !sawReadInCurrentInterval ) {
            pushToEmitQueue(new EmittableRead(*this, loading_interval_data, read));
        } else
            loading_interval_data->readsNotToClean.push_back(read);
    }
    else if ( readLoc.overlapsP(loading_interval_data->current_interval) ) {
        sawReadInCurrentInterval = true;
        
        if ( doNotTryToClean(*read) ) {
            loading_interval_data->readsNotToClean.push_back(read);
        } else {
            loading_interval_data->readsToClean.add(read);
            
            // add the rods to the list of known variants
            populateKnownIndels(*loading_interval_data, metaDataTracker);
        }
        
        if ( loading_interval_data->readsToClean.size() + loading_interval_data->readsNotToClean.size() >= MAX_READS ) {
            cerr << "Not attempting realignment in interval " << loading_interval_data->current_interval.toString().c_str() << " because there are too many reads(more than " << MAX_READS << ")." << endl;
            pushToEmitQueue(new EmittableReadList(*this, loading_interval_data));
            flushEmitQueue();
            ++interval_it;
            if(interval_it != intervalsFile.end()) {
                current_interval = *interval_it;
            } else {
                current_interval = NULL;
            }
            
            loading_interval_data = new IntervalData(current_interval);
            loading_interval_data->readsToClean.initialize(loc_parser, &sequence_dictionary);

            sawReadInCurrentInterval = false;
        }
    }
    else {  // the read is past the current interval
        CleanAndEmitReadList * read_list = new CleanAndEmitReadList(*this, loading_interval_data);
        pushToEmitQueue(read_list);
        if(nothreads)
            CleanJob(read_list).runJob();
        else
            ThreadPool::sharedPool()->addJob(new CleanJob(read_list, true));
        flushEmitQueue();

        do {
            ++interval_it;
            if(interval_it != intervalsFile.end()) {
                current_interval = *interval_it;
            } else {
                current_interval = NULL;
            }
        } while ( current_interval != NULL && current_interval->isBefore(readLoc) );
        loading_interval_data = new IntervalData(current_interval);
        loading_interval_data->readsToClean.initialize(loc_parser, &sequence_dictionary);

        sawReadInCurrentInterval = false;
        map_func( read, metaDataTracker);
    }
    
    return 0;
}

bool LocalRealignment::doNotTryToClean( const OGERead & read) {
    const int NO_ALIGNMENT_START = -1;
    
    bool is454 = false;
    if(read.HasTag("RG")) {
        string rg_val;
        read.GetTag("RG", rg_val);
        is454 = NULL != strstr( rg_val.c_str(), "454");
    }

    return !read.IsMapped() ||
    !read.IsPrimaryAlignment() ||
    read.IsFailedQC() ||
    read.getMapQuality() == 0 ||
    read.getPosition() == NO_ALIGNMENT_START ||

    ConstrainedMateFixingManager::iSizeTooBigToMove(read, MAX_ISIZE_FOR_MOVEMENT) ||
    is454;

    // TODO -- it would be nice if we could use indels from 454 reads as alternate consenses
}

void LocalRealignment::onTraversalDone(IntervalData & interval_data, int result) {
    //wait for emits to finish
    bool finished = false;
    while(!finished) {
        //if locking this mutex fails, we are busy elsewhere, so we know that we are for sure not finished.
        if(emit_mutex.try_lock()) {
            finished = emit_queue.empty();
            emit_mutex.unlock();
        }
        
        if(!finished) {
            usleep(20000);
            flushEmitQueue();
        }
    }
    
    if ( interval_data.readsToClean.size() > 0 ) {
        OGERead * read = interval_data.readsToClean.getReads()[0];
        GenomeLoc earliestPossibleMove = loc_parser->createGenomeLoc(*read);
        if ( manager->canMoveReads(earliestPossibleMove) )
            clean(interval_data);
        emitReadLists(interval_data);
    } else if ( interval_data.readsNotToClean.size() > 0 ) {
        emitReadLists(interval_data);
    }
    
    interval_data.knownIndelsToTry.clear();
    interval_data.indelRodsSeen.clear();
    if(!nothreads)
    ThreadPool::sharedPool()->waitForJobCompletion();
    manager->close();

#ifdef LR_SUPPORT_ADDITIONAL_OUTPUT_FILES
    if ( write_out_indels )
        indelOutput.close();
    if ( write_out_stats )
        statsOutput.close();
    if ( write_out_snps ) 
        snpsOutput.close();
#endif
    
    delete manager;
    
    for(vector<GenomeLoc *>::const_iterator interval_it = intervalsFile.begin(); interval_it != intervalsFile.end(); interval_it++)
        delete *interval_it;
}

void LocalRealignment::populateKnownIndels(IntervalData & interval_data, const ReadMetaDataTracker & metaDataTracker) {
    map<int, set<GATKFeature> > contigOffsetMapping = metaDataTracker.getContigOffsetMapping();

    for ( map<int, set<GATKFeature> >::iterator rods_it  = contigOffsetMapping.begin(); rods_it != contigOffsetMapping.end(); rods_it++) {
        for(set<GATKFeature>::iterator rodIter = rods_it->second.begin(); rodIter != rods_it->second.end();rodIter++) {
            const GATKFeature & feature = *rodIter;
            const GATKFeature & rod = feature;    // LCB TODO verify the correctness of this function vs original source code
            if ( interval_data.indelRodsSeen.count(rod) > 0 )
                continue;
            interval_data.indelRodsSeen.insert(rod);
            //LCB TODO do we need to include the next two lines? What is the significance of GATK's child being a VariantContext??
            //if ( dynamic_cast<VariantContext *>(&rod) ) //we need to statically decide this in OGE
            //    knownIndelsToTry.push_back(*(dynamic_cast<const VariantContext*>(&rod)));
        }
    }
}

int LocalRealignment::mismatchQualitySumIgnoreCigar(AlignedRead & aRead, const string & refSeq, int refIndex, int quitAboveThisValue) {

    const string & readSeq = aRead.getReadBases();
    const string & quals = aRead.getBaseQualities();
    int sum = 0;
    size_t readSeq_size = readSeq.size();
    size_t refSeq_size = refSeq.size();
    const char * p_refSeq = refSeq.c_str();
    const char * p_readSeq = readSeq.c_str();
    const char * p_quals = quals.c_str();
    
    int common_length = min(readSeq_size, refSeq_size - refIndex - 1);
    int readIndex = 0;

    for ( ; readIndex < common_length && sum <= quitAboveThisValue ; refIndex++, readIndex++ ) {
        char refChr = p_refSeq[refIndex];
        char readChr = p_readSeq[readIndex];
        if ( !BaseUtils::isRegularBase(readChr) || !BaseUtils::isRegularBase(refChr) ) continue; // do not count Ns/Xs/etc ?
        
        if ( readChr != refChr ) {
            sum += (int)p_quals[readIndex] -33;   // Our qualities are still stored in ASCII
        }
    }

    for (; readIndex < readSeq_size && sum <= quitAboveThisValue; refIndex++, readIndex++ ) {
        if ( refIndex >= refSeq_size ) {
            sum += MAX_QUAL;
        } else {
            char refChr = p_refSeq[refIndex];
            char readChr = p_readSeq[readIndex];
            if ( !BaseUtils::isRegularBase(readChr) || !BaseUtils::isRegularBase(refChr) ) continue; // do not count Ns/Xs/etc ?

            if ( readChr != refChr ) {
                sum += (int)p_quals[readIndex] -33;   // Our qualities are still stored in ASCII
            }
        }
    }
    return sum;
}

void LocalRealignment::clean(IntervalData & interval_data) const {
    
    vector<OGERead *> reads = interval_data.readsToClean.getReads();
    if ( reads.size() == 0 )
        return;
    
    string reference = interval_data.readsToClean.getReference(*referenceReader);
    int leftmostIndex = interval_data.readsToClean.getLocation().getStart();
    
    vector<OGERead *> refReads;                 // reads that perfectly match ref
    vector<AlignedRead *> altReads;               // reads that don't perfectly match
    vector<AlignedRead *> altAlignmentsToTest;  // should we try to make an alt consensus from the read?
    set<Consensus *> altConsenses;               // list of alt consenses
    
    // if there are any known indels for this region, get them and create alternate consenses
    generateAlternateConsensesFromKnownIndels(interval_data, altConsenses, leftmostIndex, reference);

    // decide which reads potentially need to be cleaned;
    // if there are reads with a single indel in them, add that indel to the list of alternate consenses
    long totalRawMismatchSum = determineReadsThatNeedCleaning(reads, refReads, altReads, altAlignmentsToTest, altConsenses, leftmostIndex, reference);
    
    if ( verbose  && !altConsenses.empty()) cerr << "Checking consenses for " << interval_data.current_interval.toString() << "...(" << altConsenses.size() << " consensuses across " << altReads.size() << " reads)" << endl;
    
    timeval start_time;
    gettimeofday(&start_time, NULL);
    
    Consensus * bestConsensus = NULL;

    // Randomize order of consensuses to try to avoid situations where we evaluate all the most expensive consenses first.
    // Instead, we would like to increase our chances of encountering an easy consensus early to provide a lower early-out value more often.
    vector<Consensus *> altConsensusRandomOrder(altConsenses.begin(), altConsenses.end());
    random_shuffle(altConsensusRandomOrder.begin(), altConsensusRandomOrder.end() );
    for (vector<Consensus *>::iterator iter = altConsensusRandomOrder.begin(); iter != altConsensusRandomOrder.end(); iter++) {
        Consensus &consensus = **iter;
        //if(verbose)
            //cerr << "Trying new consensus: " << cigarToString( consensus.cigar) /*<< " " << consensus.str*/ << endl;
        
        consensus.cigar.size();
        if ( false ) {
            cerr << "Checking consensus with alignment at " << consensus.positionOnReference << " cigar " << cigarToString(consensus.cigar) << endl;
            //cerr << consensus.str << endl;
            int z = 0;
            for ( ; z < consensus.positionOnReference; z++ ) 
                cerr << ".";
            for ( z=0 ; z < consensus.cigar[0].length ; z++ )
                cerr << ".";
            if ( consensus.cigar[1].type == 'I' )
                for ( z= 0; z < consensus.cigar[1].length; z++ )  
                    cerr << "I";
            cerr << endl;
        }
        
        //cerr << "Consensus: "<< consensus.str << endl;

        for ( int j = 0; j < altReads.size(); j++ ) {
            AlignedRead & toTest = *altReads[j];
            pair<int, int> altAlignment = findBestOffset(consensus.str, toTest, leftmostIndex);
            
            // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
            int myScore = altAlignment.second;
            
            //cerr << "Orig myscore " << myScore << " aligner " << toTest.getAlignerMismatchScore() << " ref " << toTest.getMismatchScoreToReference() << endl;
            
            if ( myScore > toTest.getAlignerMismatchScore() || myScore >= toTest.getMismatchScoreToReference() ) {
                myScore = toTest.getMismatchScoreToReference();
            }
            // keep track of reads that align better to the alternate consensus.
            // By pushing alignments with equal scores to the alternate, it means we'll over-call (het -> hom non ref) but are less likely to under-call (het -> ref, het non ref -> het)
            else
                consensus.readIndexes.push_back( pair<int, int>(j, altAlignment.first));
            
            //cerr << cigarToString(consensus.cigar) << " vs. " << toTest.getRead()->Name << "-" << toTest.getRead()->QueryBases << " => " << myScore << " vs. " << toTest.getMismatchScoreToReference() << endl;
            if ( !toTest.getRead()->IsDuplicate() )
                consensus.mismatchSum += myScore;
            
            // optimization: once the mismatch sum is higher than the best consensus, quit since this one can't win
            //  THIS MUST BE DISABLED IF WE DECIDE TO ALLOW MORE THAN ONE ALTERNATE CONSENSUS!
            if ( bestConsensus != NULL && consensus.mismatchSum > bestConsensus->mismatchSum )
                break;
        }
        
        //if(verbose)
        //    cerr << "Mismatch sum of new consensus " << cigarToString(consensus.cigar) << ": " << consensus.mismatchSum << endl;
        if ( bestConsensus == NULL || bestConsensus->mismatchSum > consensus.mismatchSum) {
            // we do not need this alt consensus, release memory right away!!
            if ( bestConsensus != NULL ) {
                bestConsensus->readIndexes.clear();
            }
            bestConsensus = &consensus;
            //if(verbose)
            //    cerr << "New consensus " << cigarToString(bestConsensus->cigar) <<  " is now best consensus" << endl;
        } else {
            // we do not need this alt consensus, release memory right away!!
            consensus.readIndexes.clear();
        }
    }
    
    if(verbose && !altConsenses.empty()) {
        cerr << "Scores: ";
        set<Consensus *, ConsensusScoreComparator> consenses_sorted_by_score(altConsenses.begin(), altConsenses.end());
        for(set<Consensus *>::const_iterator i = consenses_sorted_by_score.begin(); i != consenses_sorted_by_score.end(); i ++)
            cerr << cigarToString((*i)->cigar) << " (" << (*i)->mismatchSum << ") ";
        cerr << endl;
        
        timeval stop_time, real_time;
        gettimeofday(&stop_time, NULL);
        real_time.tv_sec = stop_time.tv_sec - start_time.tv_sec;
        real_time.tv_usec = stop_time.tv_usec - start_time.tv_usec;
        
        float time = float(real_time.tv_sec ) + (1.e-6 * real_time.tv_usec);
        if(verbose && time > 0.5)
            cerr << "Elapsed time = " << time << "s" << endl;
    }
    // if:
    // 1) the best alternate consensus has a smaller sum of quality score mismatches than the aligned version of the reads,
    // 2) beats the LOD threshold for the sum of quality score mismatches of the raw version of the reads,
    // 3) didn't just move around the mismatching columns (i.e. it actually reduces entropy), 
    // then clean!
    const double improvement = (bestConsensus == NULL ? -1 : ((double)(totalRawMismatchSum - bestConsensus->mismatchSum))/10.0);
    if ( improvement >= LOD_THRESHOLD ) {
        assert(bestConsensus != NULL);
        bestConsensus->cigar = AlignmentUtils::leftAlignIndel(bestConsensus->cigar, reference, bestConsensus->str, bestConsensus->positionOnReference, bestConsensus->positionOnReference);
        
        // start cleaning the appropriate reads
        for ( vector<pair<int, int> >::iterator indexPair = bestConsensus->readIndexes.begin(); indexPair !=  bestConsensus->readIndexes.end(); indexPair++) {
            AlignedRead & aRead = *altReads[indexPair->first];
            if ( !updateRead(bestConsensus->cigar, bestConsensus->positionOnReference, indexPair->second, aRead, leftmostIndex) )
                return;
        }
        if ( consensusModel != KNOWNS_ONLY && !alternateReducesEntropy(altReads, reference, leftmostIndex) ) {
#ifdef LR_SUPPORT_ADDITIONAL_OUTPUT_FILES
            if ( write_out_stats ) {
                statsOutput << interval_data.current_interval.toString();
                statsOutput << "\tFAIL (bad indel)\t"; // if improvement > LOD_THRESHOLD *BUT* entropy is not reduced (SNPs still exist)
                statsOutput << improvement << endl;
            }
#endif
        } else {
            assert (bestConsensus->cigar.size() < 100);
            //cerr << "CLEAN: " << cigarToString(bestConsensus->cigar) << " " << string(bestConsensus->str) << " " << bestConsensus->cigar.size() << endl;
#ifdef LR_SUPPORT_ADDITIONAL_OUTPUT_FILES
            if ( outputIndels && bestConsensus->cigar.size() > 1 ) {
                // NOTE: indels are printed out in the format specified for the low-coverage pilot1
                //  indel calls (tab-delimited): chr position size type sequence
                stringstream str;

                str << sequence_dictionary[reads[0]->getRefID()].Name;

                int position = bestConsensus->positionOnReference + bestConsensus->cigar[0].length;
                str << "\t" << (leftmostIndex + position - 1);
                CigarOp ce = bestConsensus->cigar[1];
                str << "\t"  << ce.length << "\t" << ce.type << "\t";
                int length = ce.length;
                if ( ce.type == 'D' ) {
                    for ( int i = 0; i < length; i++)
                        str << reference[position+i] ;
                } else {
                    for ( int i = 0; i < length; i++)
                        str << (char)bestConsensus->str[position+i];
                }
                str << "\t" << (((double)(totalRawMismatchSum - bestConsensus->mismatchSum))/10.0) << endl;
                indelOutput << str;
            }
            if ( write_out_stats ) {
                statsOutput << interval_data.current_interval.toString() << "\tCLEAN"; // if improvement > LOD_THRESHOLD *AND* entropy is reduced
                if ( bestConsensus->cigar.size() > 1 )
                    statsOutput << " (found indel)\t" << improvement << endl;
            }
#endif
            
            // finish cleaning the appropriate reads
            for ( vector<pair<int, int> >::iterator indexPair = bestConsensus->readIndexes.begin(); indexPair !=  bestConsensus->readIndexes.end(); indexPair++ ) {
                AlignedRead & aRead = *altReads[indexPair->first];
                if ( aRead.constizeUpdate() ) {
                    // We need to update the mapping quality score of the cleaned reads;
                    // however we don't have enough info to use the proper MAQ scoring system.
                    // For now, we will just arbitrarily add 10 to the mapping quality. [EB, 6/7/2010].
                    // TODO -- we need a better solution here
                    OGERead * read = aRead.getRead();  //TODO LCB verify the below mapping quality, since BamTools and SAM format are a bit different.
                    if ( read->getMapQuality() != 255 ) // 255 == Unknown, so don't modify it
                        read->setMapQuality(min(aRead.getRead()->getMapQuality() + 10, 254));
                    
                    // before we fix the attribute tags we first need to make sure we have enough of the reference sequence
                    int neededBasesToLeft = leftmostIndex - read->getPosition();
                    int neededBasesToRight = read->getPosition() + read->getLength() - leftmostIndex - reference.size() + 1;
                    int neededBases = max(neededBasesToLeft, neededBasesToRight);
                    if ( neededBases > 0 ) {
                        int padLeft = max(leftmostIndex-neededBases, 1);
                        int padRight = min(leftmostIndex+reference.size()+neededBases, (unsigned long)referenceReader->getSequenceDictionary()[interval_data.current_interval.getContig()].getLength());
                        reference = referenceReader->getSubsequenceAt(interval_data.current_interval.getContig(), padLeft, padRight);
                        leftmostIndex = padLeft;
                    }
                    
                    // now, fix the attribute tags
                    if ( read->HasTag("NM") )
                        read->EditTag("NM", "i", SequenceUtil::calculateSamNmTag(read, reference, leftmostIndex-1));
                    if ( read->HasTag("UQ") )
                        read->EditTag("UQ", "i", SequenceUtil::sumQualitiesOfMismatches(read, reference, leftmostIndex-1));

                    // TODO -- this is only temporary until Tim adds code to recalculate this value
                    if ( read->HasTag("MD") )
                        read->RemoveTag("MD");
                    
                    // mark that it was actually cleaned
                    interval_data.readsActuallyCleaned.insert(read);
                }
            }
        }
        
        // END IF ( improvement >= LOD_THRESHOLD )
        
    }
    
    for (set<Consensus *>::iterator iter = altConsenses.begin(); iter != altConsenses.end(); iter++)
        delete *iter;
    for(vector<AlignedRead *>::const_iterator iter = altReads.begin(); iter != altReads.end(); iter++)
        delete *iter;
#ifdef LR_SUPPORT_ADDITIONAL_OUTPUT_FILES
    else if ( write_out_stats ) {
        statsOutput << interval_data.current_interval.toString() << "\tFAIL\t" << improvement << endl;
    }
#endif
}

void LocalRealignment::generateAlternateConsensesFromKnownIndels(IntervalData & interval_data, set<Consensus *> & altConsensesToPopulate, const int leftmostIndex, const string reference) const {
    for ( vector<VariantContext>::iterator knownIndelIt  = interval_data.knownIndelsToTry.begin(); knownIndelIt != interval_data.knownIndelsToTry.end(); knownIndelIt++ ) {
        VariantContext * knownIndel = &*knownIndelIt;
        if ( knownIndel == NULL || !knownIndel->isIndel() || knownIndel->isComplexIndel() )
            continue;
        string indelStr = knownIndel->isSimpleInsertion() ? knownIndel->getAlternateAllele(0).getBases() : string( knownIndel->getReference().length(), '-');  //string constructor used here makes n copies of one character
        int start = knownIndel->getStart() - leftmostIndex + 1;
        Consensus * c = createAlternateConsensus(start, reference, indelStr, *knownIndel);
        if ( c != NULL )
            altConsensesToPopulate.insert(c);
    }
}

long LocalRealignment::determineReadsThatNeedCleaning( vector<OGERead *> & reads,
                                            vector<OGERead *> & refReadsToPopulate,
                                            vector<AlignedRead *> & altReadsToPopulate,
                                            vector<AlignedRead *> & altAlignmentsToTest,
                                            set<Consensus *> & altConsenses,
                                            int leftmostIndex,
                                            string & reference) const{
    
    long totalRawMismatchSum = 0L;
    
    for ( int read_ctr = 0; read_ctr < reads.size(); read_ctr++ ) {
        OGERead * read = reads[read_ctr];
        
        // we can not deal with screwy records
        if ( read->getCigarData().size() == 0 ) {
            refReadsToPopulate.push_back(read);
            continue;
        }
        
        AlignedRead * aRead = new AlignedRead(read, &sequence_dictionary);
        
        // first, move existing indels (for 1 indel reads only) to leftmost position within identical sequence
        int numBlocks = 0;

        vector<CigarOp> cigar = read->getCigarData();
        for(vector<CigarOp>::const_iterator i = cigar.begin(); i != cigar.end(); i++)
            if(i->type == 'M' || i->type == '=' || i->type == 'X')
                numBlocks++;

        if ( numBlocks == 2 ) {
            vector<CigarOp> newCigar = AlignmentUtils::leftAlignIndel(unclipCigar(read->getCigarData()), reference, read->getQueryBases(), read->getPosition()-leftmostIndex, 0);
            aRead->setCigar(newCigar, false);
        }
        
        const int startOnRef = read->getPosition()-leftmostIndex;
        const int rawMismatchScore = mismatchQualitySumIgnoreCigar(*aRead, reference, startOnRef, INT_MAX);
        
        // if this doesn't match perfectly to the reference, let's try to clean it
        if ( rawMismatchScore > 0 ) {
            //if(verbose)
            //    cerr << "Adding " << read->Name << " with raw mismatch score " << rawMismatchScore << " to non-ref reads" << endl;
            
            if ( !read->IsDuplicate() )
                totalRawMismatchSum += rawMismatchScore;
            aRead->setMismatchScoreToReference(rawMismatchScore);
            aRead->setAlignerMismatchScore(AlignmentUtils::mismatchingQualities(aRead->getRead(), reference, startOnRef));
            
            // if it has an indel, let's see if that's the best consensus
            if ( consensusModel != KNOWNS_ONLY && numBlocks == 2 )  {
                Consensus * c = createAlternateConsensus(startOnRef, aRead->getCigar(), reference, aRead->getReadBases());
                if ( c != NULL ) {
                    bool string_exists_in_other_consensus = false;
                    for(set<Consensus *>::const_iterator i = altConsenses.begin(); i != altConsenses.end(); i++)
                        if(c->str == (*i)->str)
                            string_exists_in_other_consensus = true;
                    
                    //if(verbose)
                    //    cerr << "New consensus " << cigarToString(c->cigar) << endl;//" with string " << c->str << endl;
                    if(!string_exists_in_other_consensus)
                        altConsenses.insert(c);
                    else
                        delete c;
                } //else
                    //if(verbose)
                    //    cerr << "No new consensus" << endl;
            } else {
                altAlignmentsToTest.push_back(aRead);
            }
            
            altReadsToPopulate.push_back(aRead);
        }
        // otherwise, we can emit it as is
        else {
            delete aRead;
            //if(verbose)
            //    cerr << "Adding " << read->Name << " with raw mismatch score " << rawMismatchScore << " to ref reads" << endl;
            refReadsToPopulate.push_back(read);
        }
    }
    
    return totalRawMismatchSum;
}

void LocalRealignment::generateAlternateConsensesFromReads( vector<AlignedRead> & altAlignmentsToTest,
                                                 set<Consensus *> & altConsensesToPopulate,
                                                 const string & reference,
                                                 const int leftmostIndex) {
    
    // if we are under the limit, use all reads to generate alternate consenses
    if ( altAlignmentsToTest.size() <= MAX_READS_FOR_CONSENSUSES ) {
    }
    // otherwise, choose reads for alternate consenses randomly
    else {
        int readsSeen = 0;
        while ( readsSeen++ < MAX_READS_FOR_CONSENSUSES && altConsensesToPopulate.size() <= MAX_CONSENSUSES) {
            int index = rand() % altAlignmentsToTest.size();
            vector<AlignedRead>::iterator to_remove = altAlignmentsToTest.begin() + index;
            AlignedRead aRead = *to_remove;
            altAlignmentsToTest.erase(to_remove);
        }
    }
}

// create a Consensus from cigar/read strings which originate somewhere on the reference
LocalRealignment::Consensus * LocalRealignment::createAlternateConsensus(const int indexOnRef, const vector<CigarOp> & c, const string reference, const string readStr) const {
    if ( indexOnRef < 0 )
        return NULL;
    
    // if there are no indels, we do not need this consensus, can abort early:
    if ( c.size() == 1 && c[0].type == 'M' ) return NULL;
    
    // create the new consensus
    vector<CigarOp> elements;
    elements.reserve(c.size()-1);
    stringstream sb;
    for (int i = 0; i < indexOnRef; i++)
        sb << reference[i];
    
    int indelCount = 0;
    int altIdx = 0;
    int refIdx = indexOnRef;
    bool ok_flag = true;
    for ( int i = 0 ; i < c.size() ; i++ ) {
        CigarOp ce = c[i];
        int elementLength = ce.length;
        switch( ce.type ) {
            case 'D':
                refIdx += elementLength;
                indelCount++;
                elements.push_back(ce);
                break;
            case 'M':
                altIdx += elementLength;
            case 'N':
                if ( reference.size() < refIdx + elementLength )
                    ok_flag = false;
                else  {
                    for (int j = 0; j < elementLength; j++)
                        sb << reference[refIdx+j];
                }
                refIdx += elementLength;
                elements.push_back(CigarOp('M', elementLength));
                break;
            case 'I':
                for (int j = 0; j < elementLength; j++) {
                    if ( ! BaseUtils::isRegularBase(readStr[altIdx+j]) ) {
                        // Insertions with N's in them cause real problems sometimes; it's better to drop them altogether
                        ok_flag=false;
                        break;
                    }
                    sb << readStr[altIdx + j];
                }
                altIdx += elementLength;
                indelCount++;
                elements.push_back(ce);
                break;
            case 'S':
            default:
                break;
        }
    }
    // make sure that there is at most only a single indel and it aligns appropriately!
    if ( !ok_flag || indelCount != 1 || reference.size() < refIdx )
        return NULL;
    
    for (int i = refIdx; i < reference.size(); i++)
        sb << reference[i];
    string altConsensus =  sb.str(); // alternative consensus sequence we just built from the current read
    
    return new Consensus(altConsensus, vector<CigarOp> (elements), indexOnRef);
}

// create a Consensus from just the indel string that falls on the reference
LocalRealignment::Consensus * LocalRealignment::createAlternateConsensus(const int indexOnRef, const string & reference, const string & indelStr, VariantContext indel) const {
    if ( indexOnRef < 0 || indexOnRef >= reference.size() )
        return NULL;
    
    // create the new consensus
    stringstream sb;
    vector<CigarOp> cigar;
    int refIdx;
    
    for (refIdx = 0; refIdx < indexOnRef; refIdx++)
        sb << reference[refIdx];
    if ( indexOnRef > 0 )
        cigar.push_back(CigarOp('M', indexOnRef));
    
    if ( indel.isSimpleDeletion() ) {
        refIdx += indelStr.size();
        cigar.push_back(CigarOp('D',indelStr.size()));
    }
    else if ( indel.isSimpleInsertion() ) {
        sb << indelStr;
        cigar.push_back( CigarOp('I', indelStr.size()));
    } else {
        cerr << "Creating an alternate consensus from a complex indel is not allows" << endl;
    }
    
    if ( reference.size() - refIdx > 0 )
        cigar.push_back( CigarOp('M',reference.size() - refIdx));        
    for (; refIdx < reference.size(); refIdx++)
        sb << reference[refIdx];
    string altConsensus = sb.str(); // alternative consensus sequence we just built from the current read
    
    return new LocalRealignment::Consensus(altConsensus, cigar, 0);
}


pair<int, int> LocalRealignment::findBestOffset(const string & ref, AlignedRead read, const int leftmostIndex) const {
    
    // optimization: try the most likely alignment first (to get a low score to beat)
    int originalAlignment = read.getOriginalAlignmentStart() - leftmostIndex;
    int bestScore = mismatchQualitySumIgnoreCigar(read, ref, originalAlignment, INT_MAX);
    int bestIndex = originalAlignment;

    // optimization: we can't get better than 0, so we can quit now
    if ( bestScore == 0 )
        return pair<int, int>(bestIndex, 0);
    
    // optimization: the correct alignment shouldn't be too far from the original one (or else the read wouldn't have aligned in the first place)
    for ( int i = 0; i < originalAlignment; i++ ) {
        int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
        if ( score < bestScore ) {
            bestScore = score;
            bestIndex = i;
        }
        // optimization: we can't get better than 0, so we can quit now
        if ( bestScore == 0 )
            return pair<int, int>(bestIndex, 0);
    }

    const size_t length = read.getCigarLength();
    const int maxPossibleStart = ref.size() - length;

    for ( int i = originalAlignment + 1; i <= maxPossibleStart; i++ ) {
        int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
        if ( score < bestScore ) {
            bestScore = score;
            bestIndex = i;
        }
        // optimization: we can't get better than 0, so we can quit now
        if ( bestScore == 0 )
            return pair<int, int>(bestIndex, 0);
    }

    return pair<int, int>(bestIndex, bestScore);
}

bool LocalRealignment::updateRead(const vector<CigarOp> & altCigar, const int altPosOnRef, const int myPosOnAlt, AlignedRead & aRead, const int leftmostIndex) const {
    vector<CigarOp> readCigar;
    
    // special case: there is no indel
    if ( altCigar.size() == 1 ) {
        aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
        readCigar.push_back( CigarOp('M', aRead.getReadLength()));
        aRead.setCigar(readCigar);
        return true;
    }
    
    CigarOp altCE1 = altCigar[0];
    CigarOp altCE2 = altCigar[1];
    
    int leadingMatchingBlockLength = 0; // length of the leading M element or 0 if the leading element is I
    
    CigarOp indelCE;
    if ( altCE1.type == 'I'  ) {
        indelCE=altCE1;
        if ( altCE2.type != 'M' ) {
            cerr << "When the first element of the alt consensus is I, the second one must be M. Actual: " << cigarToString(altCigar) << ".  Skipping this site..." << endl;
            return false;
        }
    }
    else {
        if ( altCE1.type != 'M'  ) {
            cerr << "First element of the alt consensus cigar must be M or I. Actual: " << cigarToString(altCigar) << ".  Skipping this site..." << endl;
            return false;
        }
        if ( altCE2.type == 'I'  || altCE2.type == 'D' ) {
            indelCE=altCE2;
        } else {
            cerr << "When first element of the alt consensus is M, the second one must be I or D. Actual: " << cigarToString(altCigar) << ".  Skipping this site..." << endl;
            return false;
        }
        leadingMatchingBlockLength = altCE1.length;
    }
    
    // the easiest thing to do is to take each case separately
    int endOfFirstBlock = altPosOnRef + leadingMatchingBlockLength;
    bool sawAlignmentStart = false;
    
    // for reads starting before the indel
    if ( myPosOnAlt < endOfFirstBlock) {
        aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
        sawAlignmentStart = true;
        
        // for reads ending before the indel
        if ( myPosOnAlt + aRead.getReadLength() <= endOfFirstBlock) {
            //readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            //aRead.setCigar(readCigar);
            aRead.clearCigar(); // reset to original alignment
            return true;
        }
        readCigar.push_back( CigarOp('M', endOfFirstBlock - myPosOnAlt));
    }
    
    // forward along the indel
    //int indelOffsetOnRef = 0, indelOffsetOnRead = 0;
    if ( indelCE.type == 'I' ) {
        // for reads that end in an insertion
        if ( myPosOnAlt + aRead.getReadLength() < endOfFirstBlock + indelCE.length ) {
            int partialInsertionLength = myPosOnAlt + aRead.getReadLength() - endOfFirstBlock;
            // if we also started inside the insertion, then we need to modify the length
            if ( !sawAlignmentStart )
                partialInsertionLength = aRead.getReadLength();
            readCigar.push_back( CigarOp('I', partialInsertionLength));
            aRead.setCigar(readCigar);
            return true;
        }
        
        // for reads that start in an insertion
        if ( !sawAlignmentStart && myPosOnAlt < endOfFirstBlock + indelCE.length ) {
            aRead.setAlignmentStart(leftmostIndex + endOfFirstBlock);
            readCigar.push_back(CigarOp('I', indelCE.length - (myPosOnAlt - endOfFirstBlock)));
            //indelOffsetOnRead = myPosOnAlt - endOfFirstBlock;
            sawAlignmentStart = true;
        } else if ( sawAlignmentStart ) {
            readCigar.push_back(indelCE);
            //indelOffsetOnRead = indelCE.getLength();
        }
    } else if ( indelCE.type == 'D' ) {
        if ( sawAlignmentStart )
            readCigar.push_back(indelCE);
        //indelOffsetOnRef = indelCE.getLength();
    }
    
    // for reads that start after the indel
    if ( !sawAlignmentStart ) {
        //aRead.setAlignmentStart(leftmostIndex + myPosOnAlt + indelOffsetOnRef - indelOffsetOnRead);
        //readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
        //aRead.setCigar(readCigar);
        aRead.clearCigar(); // reset to original alignment
        return true;
    }
    
    int readRemaining = aRead.getReadBases().size();
    for ( vector<CigarOp>::iterator ce = readCigar.begin(); ce != readCigar.end(); ce++ ) {
        if ( ce->type != 'D' )
            readRemaining -= ce->length;
    }
    if ( readRemaining > 0 )
        readCigar.push_back(CigarOp('M', readRemaining));
    aRead.setCigar(readCigar);
    
    return true;
}

bool LocalRealignment::alternateReducesEntropy(const vector<AlignedRead *> & reads, const string & reference, const int leftmostIndex) const {
    const size_t array_size = reference.size();
    int originalMismatchBases[array_size];
    int cleanedMismatchBases[array_size];
    int totalOriginalBases[array_size];
    int totalCleanedBases[array_size];
    
    // set to 1 to prevent dividing by zero
    for ( int i=0; i < array_size; i++ ) {
        originalMismatchBases[i] = 0;
        totalOriginalBases[i] = 0;
        cleanedMismatchBases[i] = 0;
        totalCleanedBases[i] = 0;
    }
    
    for (int i=0; i < reads.size(); i++) {
        AlignedRead &read = *reads[i];
        
        int count_m_eq_x = 0;

        vector<CigarOp> cigar = read.getRead()->getCigarData();
        for(vector<CigarOp>::const_iterator i = cigar.begin(); i != cigar.end(); i++)
            if(i->type == 'M' || i->type == '=' || i->type == 'X')
                count_m_eq_x++;

        //if ( read.getRead()->getAlignmentBlocks().size() > 1 )
        if(count_m_eq_x > 1)
            continue;
        
        int refIdx = read.getOriginalAlignmentStart() - leftmostIndex;
        const string readStr = read.getReadBases();
        const string quals = read.getBaseQualities();
        
        for (int j=0; j < readStr.size(); j++, refIdx++ ) {
            if ( refIdx < 0 || refIdx >= reference.size() ) {
                //System.out.println( "Read: "+read.getRead().getReadName() + "; length = " + readStr.length() );
                //System.out.println( "Ref left: "+ leftmostIndex +"; ref length=" + reference.length() + "; read alignment start: "+read.getOriginalAlignmentStart() );
                break;
            }
            totalOriginalBases[refIdx] += quals[j] - 33;
            if ( readStr[j] != reference[refIdx] )
                originalMismatchBases[refIdx] += quals[j] - 33;
        }
        
        // reset and now do the calculation based on the cleaning
        refIdx = read.getAlignmentStart() - leftmostIndex;
        int altIdx = 0;
        vector<CigarOp> c = read.getCigar();
        for (int j = 0 ; j < c.size() ; j++) {
            CigarOp ce = c[j];
            int elementLength = ce.length;
            switch ( ce.type ) {
                case 'M':
                    for (int k = 0 ; k < elementLength ; k++, refIdx++, altIdx++ ) {
                        if ( refIdx >= reference.size() )
                            break;
                        totalCleanedBases[refIdx] += quals[altIdx] - 33;
                        if ( readStr[altIdx] != reference[refIdx] )
                            cleanedMismatchBases[refIdx] += quals[altIdx] - 33;
                    }
                    break;
                case 'I':
                    altIdx += elementLength;
                    break;
                case 'D':
                    refIdx += elementLength;
                    break;
                case 'S':
                default:
                    break;
            }
        }
    }
    
    int originalMismatchColumns = 0, cleanedMismatchColumns = 0;
    stringstream sb;
    for ( int i=0; i < array_size; i++ ) {
        if ( cleanedMismatchBases[i] == originalMismatchBases[i] )
            continue;
        bool didMismatch = false, stillMismatches = false;
        if ( originalMismatchBases[i] > totalOriginalBases[i] * MISMATCH_THRESHOLD )  {
            didMismatch = true;
            originalMismatchColumns++;
            if ( totalCleanedBases[i] > 0 && ((double)cleanedMismatchBases[i] / (double)totalCleanedBases[i]) > ((double)originalMismatchBases[i] / (double)totalOriginalBases[i]) * (1.0 - MISMATCH_COLUMN_CLEANED_FRACTION) ) {
                stillMismatches = true;
                cleanedMismatchColumns++;
            }
        } else if ( cleanedMismatchBases[i] > totalCleanedBases[i] * MISMATCH_THRESHOLD ) {
            cleanedMismatchColumns++;
        }
        
#ifdef LR_SUPPORT_ADDITIONAL_OUTPUT_FILES
        if ( output_snps ) {
            if ( didMismatch ) {
                const SamSequenceDictionary & s = sequence_dictionary;
                const string & ref_name = s[reads[0]->getRead()->getRefID()].Name;

                sb << ref_name << ":" << (leftmostIndex + i);
                if ( stillMismatches )
                    sb << " SAME_SNP\n" << endl;
                else
                    sb << " NOT_SNP" << endl;
            }
        }
#endif
    }
    
    //if(verbose)
    //    cerr << "Original mismatch columns = " << originalMismatchColumns << "; cleaned mismatch columns = " << cleanedMismatchColumns << endl;
    
    const bool reduces = (originalMismatchColumns == 0 || cleanedMismatchColumns < originalMismatchColumns);
#ifdef LR_SUPPORT_ADDITIONAL_OUTPUT_FILES
    if ( reduces && output_snps ) {
        snpsOutput << sb.str();
    }
#endif
    return reduces;
}
    
vector<CigarOp> LocalRealignment::unclipCigar(const vector<CigarOp> & cigar) {
    vector<CigarOp> elements;
    elements.reserve(cigar.size());

    for (int i = 0; i < cigar.size(); i++ ) {
        CigarOp ce = cigar[i];
        if ( !isClipOperator(ce) )
            elements.push_back(ce);
    }
    return elements;
}

bool LocalRealignment::isClipOperator(const CigarOp op) {
    return op.type == 'S' || op.type == 'H' || op.type == 'P';
}

vector<CigarOp> LocalRealignment::reclipCigar(const vector<CigarOp> & cigar, OGERead * read) {
    vector<CigarOp> elements, cigarData;
    
    cigarData = read->getCigarData();
    
    int i = 0;
    int n = read->getCigarData().size();

    while ( i < n && isClipOperator(cigarData[i]) )
        elements.push_back(cigarData[i++]);
    
    //add all elements of cigar to elements
    elements.insert(elements.end(), cigar.begin(), cigar.end());    
    
    i++;

    while ( i < n && !isClipOperator(cigarData[i]) )
        i++;
    
    while ( i < n && isClipOperator(cigarData[i]) )
        elements.push_back(cigarData[i++]);
    
    return elements;
}

LocalRealignment::LocalRealignment() 
: LOD_THRESHOLD(5.0)
, manager(NULL)
, consensusModel(USE_READS)
, MISMATCH_THRESHOLD(.15)
, MAX_RECORDS_IN_MEMORY(150000)
, MAX_ISIZE_FOR_MOVEMENT(3000)
, MAX_POS_MOVE_ALLOWED(200)
, MAX_CONSENSUSES(30)
, MAX_READS_FOR_CONSENSUSES(120)
, MAX_READS(20000)
, NO_ORIGINAL_ALIGNMENT_TAGS(false)
#ifdef LR_SUPPORT_ADDITIONAL_OUTPUT_FILES
, write_out_indels(false)
, write_out_stats(false)
, write_out_snps(false)
#endif
, referenceReader(NULL)
, loc_parser(NULL)
, sawReadInCurrentInterval(false)
, loading_interval_data(NULL)
#ifdef LR_SUPPORT_ADDITIONAL_OUTPUT_FILES
, outputIndels(false)
, output_stats(false)
, output_snps(false)
#endif
{}

int LocalRealignment::runInternal()
{
    ogeNameThread("LRmain");
    sequence_dictionary = getHeader().getSequences();
    initialize();

    interval_it = intervalsFile.begin();

    loading_interval_data = new IntervalData(intervalsFile.empty() ? NULL : intervalsFile.front());
    loading_interval_data->readsToClean.initialize(loc_parser, &sequence_dictionary);

    while(true) {
        OGERead * al = getInputAlignment();
        
        if(!al)
            break;
        
        const ReadMetaDataTracker rmdt(loc_parser, al, std::map<int, RODMetaDataContainer>() );
        map_func( al, rmdt);
        
        //throttle ourselves, so that we don't overwhelm the downstream queue.
        while(manager->isReadAddQueueFull())
            usleep(50000);
    };
    
    onTraversalDone(*loading_interval_data, 0);
    
    return 0;
}