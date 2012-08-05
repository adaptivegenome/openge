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

/*
package org.broadinstitute.sting.gatk.walkers.indels;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.*;
import net.sf.samtools.util.RuntimeIOException;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StringUtil;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.NWaySAMFileWriter;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
*/


#include "local_realignment.h"
#include "api/algorithms/Sort.h"

using namespace BamTools;

#include <cassert>
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

string cigarToString(vector<CigarOp> cigar)
{
    stringstream ss;
    for(vector<CigarOp>::iterator i = cigar.begin(); i != cigar.end(); i++)
        ss << i->Length << i->Type;
    
    return ss.str();
}

// pull out the bases that aren't clipped out
void LocalRealignment::AlignedRead::getUnclippedBases() {
    readBases.clear();
    readBases.reserve(getReadLength());
    baseQuals.clear();
    baseQuals.reserve(getReadLength());
    const string & actualReadBases = read->QueryBases;
    const string & actualBaseQuals = read->Qualities;
    int fromIndex = 0, toIndex = 0;
    
    for ( vector<CigarOp>::iterator ce = read->CigarData.begin(); ce != read->CigarData.end(); ce++ ) {
        uint32_t & elementLength = ce->Length;
        switch ( ce->Type ) {
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
vector<CigarOp> LocalRealignment::AlignedRead::reclipCigar(const vector<CigarOp> & cigar) {
    return LocalRealignment::reclipCigar(cigar, read);
}

vector<CigarOp> & LocalRealignment::AlignedRead::getCigar() {
    return (newCigar.size() > 0 ? newCigar : read->CigarData);
}

// tentatively sets the new Cigar, but it needs to be confirmed later
void LocalRealignment::AlignedRead::setCigar(const vector<CigarOp> & cigar_in, bool fixClippedCigar) {

    bool reclip = false;
    if ( fixClippedCigar && getReadBases().length() < read->Length ) {
        reclip = true;
    }
    
    const vector<CigarOp> cigar = reclip ? reclipCigar(cigar_in) : cigar_in;

    // no change?
    if ( read->CigarData == cigar ) {
        clearCigar();
        return;
    }

    // no indel?
    bool has_i_or_d = false;
    for (vector<CigarOp>::const_iterator i = cigar.begin(); i != cigar.end(); i++)
        if(i->Type == 'I' || i->Type == 'D')
            has_i_or_d = true;

    if ( !has_i_or_d) {
        cerr << "Modifying a read with no associated indel; although this is possible, it is highly unlikely.  Perhaps this region should be double-checked: " << read->Name << " near read->getReferenceName():" << read->Position << endl;
        clearCigar();
        return;
    }
    newCigar = cigar;
}

// tentatively sets the new start, but it needs to be confirmed later
void LocalRealignment::AlignedRead::setAlignmentStart(int start) {
    newStart = start;
}

int LocalRealignment::AlignedRead::getAlignmentStart() const{
    return (newStart != -1 ? newStart : read->Position);
}

int LocalRealignment::AlignedRead::getOriginalAlignmentStart() const {
    return read->Position;
}

// constizes the changes made.
// returns true if this record actually changes, false otherwise
bool LocalRealignment::AlignedRead::constizeUpdate() {
    // if we haven't made any changes, don't do anything
    if ( newCigar.size() == 0 ) {
        return false;
    }
    if ( newStart == -1 )
        newStart = read->Position;
    else if ( abs(newStart - read->Position) > MAX_POS_MOVE_ALLOWED ) {
        cerr << "Attempting to realign read " << read->Name << " at " << read->Position << " more than " << MAX_POS_MOVE_ALLOWED << " bases to " << newStart << ".";
        return false;
    }
    
    // annotate the record with the original cigar (and optionally the alignment start)
    if ( !NO_ORIGINAL_ALIGNMENT_TAGS ) {
        read->AddTag(ORIGINAL_CIGAR_TAG, "Z", cigarToString(read->CigarData));
        if ( newStart != read->Position )
            read->AddTag(ORIGINAL_POSITION_TAG, "i", read->Position+1);
    }

    read->CigarData = newCigar;
    read->Position = newStart;
    
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
void LocalRealignment::ReadBin::add(BamAlignment * read) {
    
    GenomeLoc locForRead = loc_parser->createGenomeLoc(*read);
    if ( loc == NULL )
        loc = new GenomeLoc(locForRead);
    else if ( locForRead.getStop() > loc->getStop() )
        loc = new GenomeLoc(loc_parser->createGenomeLoc(loc->getContig(), loc->getStart(), locForRead.getStop()));

    reads.push_back(read);
}

string LocalRealignment::ReadBin::getReference(FastaReader & referenceReader) {
    // set up the reference if we haven't done so yet
    if ( reference.size() == 0 ) {
        // first, pad the reference to handle deletions in narrow windows (e.g. those with only 1 read)
        int padLeft = max(loc->getStart()-REFERENCE_PADDING, 0);    //LCB 1->0
        int padRight = min(loc->getStop()+REFERENCE_PADDING, atoi(referenceReader.getSequenceDictionary()[loc->getContig()].Length.c_str()));
        loc = new GenomeLoc(loc_parser->createGenomeLoc(loc->getContig(), padLeft, padRight));
        reference = referenceReader.readSequence(loc->getContig(), loc->getStart(), loc->getStop() - loc->getStart() + 1);

        //transform to upper case. STL sometimes makes things harder than they should be.
        std::transform(reference.begin(), reference.end(), reference.begin(), (int (*)(int))std::toupper); 
    }
    
    return reference;
}

void LocalRealignment::ReadBin::clear() {
    reads.clear();
    reference.clear();
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
    if(!referenceReader->Open(reference_filename)) {
        cerr << "Error opening reference file " << reference_filename << endl;
        exit(-1);
    }
    
    loc_parser = new GenomeLocParser(getHeader().Sequences);

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

    
    AlignedRead::NO_ORIGINAL_ALIGNMENT_TAGS = NO_ORIGINAL_ALIGNMENT_TAGS;
    AlignedRead::MAX_POS_MOVE_ALLOWED = MAX_POS_MOVE_ALLOWED;

    // set up the output writer
    manager = new ConstrainedMateFixingManager(this, MAX_ISIZE_FOR_MOVEMENT, MAX_POS_MOVE_ALLOWED, MAX_RECORDS_IN_MEMORY, loc_parser);
    
    if ( write_out_indels ) 
        indelOutput.open(OUT_INDELS.c_str());
    if ( write_out_stats )
        statsOutput.open(OUT_STATS.c_str());
    if ( write_out_snps )
        snpsOutput.open(OUT_SNPS.c_str());
}

void LocalRealignment::emit(BamAlignment * read) {
    
    // check to see whether the read was modified by looking at the temporary tag
    bool wasModified = readsActuallyCleaned.count(read) > 0;
    manager->addRead(read, wasModified);
}

void LocalRealignment::emitReadLists() {
    // pre-merge lists to sort them in preparation for constrained SAMFileWriter
    vector<BamAlignment *> to_insert = readsToClean.getReads();
    readsNotToClean.insert(readsNotToClean.end(), to_insert.begin(), to_insert.end());
    std::stable_sort( readsNotToClean.begin(), readsNotToClean.end(), Algorithms::Sort::ByPosition() );
    manager->addReads(readsNotToClean, readsActuallyCleaned);
    readsToClean.clear();
    readsNotToClean.clear();
    readsActuallyCleaned.clear();
}

int LocalRealignment::map_func( BamAlignment * read, const ReadMetaDataTracker & metaDataTracker, GenomeLoc * & currentInterval) {
    const int NO_ALIGNMENT_REFERENCE_INDEX = -1;
    if ( currentInterval == NULL ) {
        emit(read);
        return 0;
    }
    
    // edge case: when the last target interval abuts the end of the genome, we'll get one of the
    //   unmapped reads while the currentInterval still isn't null.  We need to trigger the cleaning
    //   at this point without trying to create a GenomeLoc.
    if ( read->RefID == NO_ALIGNMENT_REFERENCE_INDEX ) {
        cleanAndCallMap(NULL, *currentInterval);
        
        do {
            currentInterval = (++interval_it != intervalsFile.end()) ? *interval_it : NULL;
            
        } while ( currentInterval != NULL );
        
        sawReadInCurrentInterval = false;
        map_func(read, metaDataTracker, currentInterval);
        return 0;
    }
    GenomeLoc readLoc = loc_parser->createGenomeLoc(*read);
    // hack to get around unmapped reads having screwy locations
    if ( readLoc.getStop() == 0 )
        readLoc = loc_parser->createGenomeLoc(readLoc.getContig(), readLoc.getStart(), readLoc.getStart());
    
    if ( readLoc.isBefore(*currentInterval) ) {
        if ( !sawReadInCurrentInterval )
            emit(read);
        else
            readsNotToClean.push_back(read);
    }
    else if ( readLoc.overlapsP(*currentInterval) ) {
        sawReadInCurrentInterval = true;
        
        if ( doNotTryToClean(*read) ) {
            readsNotToClean.push_back(read);
        } else {
            readsToClean.add(read);
            
            // add the rods to the list of known variants
            populateKnownIndels(metaDataTracker);
        }
        
        if ( readsToClean.size() + readsNotToClean.size() >= MAX_READS ) {
            fprintf(stderr, "Not attempting realignment in interval %s because there are too many reads.", currentInterval->toString().c_str());
            emitReadLists();
            currentInterval = (++interval_it != intervalsFile.end()) ? *interval_it : NULL;
            sawReadInCurrentInterval = false;
        }
    }
    else {  // the read is past the current interval
        cleanAndCallMap(&readLoc, *currentInterval);
        do {
            currentInterval = (++interval_it != intervalsFile.end()) ? *interval_it : NULL;
            
        } while ( currentInterval != NULL && currentInterval->isBefore(readLoc) );
        sawReadInCurrentInterval = false;
        map_func(read, metaDataTracker, currentInterval);
    }
    
    return 0;
}

bool LocalRealignment::doNotTryToClean( const BamAlignment & read) {
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
    read.MapQuality == 0 ||
    read.Position == NO_ALIGNMENT_START ||
    ConstrainedMateFixingManager::iSizeTooBigToMove(read, MAX_ISIZE_FOR_MOVEMENT) ||
    is454;

    // TODO -- it would be nice if we could use indels from 454 reads as alternate consenses
}

void LocalRealignment::cleanAndCallMap( GenomeLoc * readLoc, const GenomeLoc & currentInterval) {
    if ( readsToClean.size() > 0 ) {
        GenomeLoc earliestPossibleMove = loc_parser->createGenomeLoc(*(readsToClean.getReads()[0]));
        if ( manager->canMoveReads(earliestPossibleMove) )
            clean(readsToClean, currentInterval);
    }
    knownIndelsToTry.clear();
    indelRodsSeen.clear();
    
    emitReadLists();
}

int LocalRealignment::reduceInit() {
    return 0;
}

int LocalRealignment::reduce(int value, int sum) {
    return sum + value;
}

void LocalRealignment::onTraversalDone(int result) {
    if ( readsToClean.size() > 0 ) {
        BamAlignment * read = readsToClean.getReads()[0];
        GenomeLoc earliestPossibleMove = loc_parser->createGenomeLoc(*read);
        if ( manager->canMoveReads(earliestPossibleMove) )
            clean(readsToClean, GenomeLoc("FinalClean", 0, 0, 0));
        emitReadLists();
    } else if ( readsNotToClean.size() > 0 ) {
        emitReadLists();                            
    }
    
    knownIndelsToTry.clear();
    indelRodsSeen.clear();
    
    if ( write_out_indels )
        indelOutput.close();
    if ( write_out_stats )
        statsOutput.close();
    if ( write_out_snps ) 
        snpsOutput.close();
    
    manager->close();
    
    for(vector<GenomeLoc *>::const_iterator interval_it = intervalsFile.begin(); interval_it != intervalsFile.end(); interval_it++)
        delete *interval_it;
}

void LocalRealignment::populateKnownIndels(const ReadMetaDataTracker & metaDataTracker) {
    map<int, set<GATKFeature> > contigOffsetMapping = metaDataTracker.getContigOffsetMapping();

    for ( map<int, set<GATKFeature> >::iterator rods_it  = contigOffsetMapping.begin(); rods_it != contigOffsetMapping.end(); rods_it++) {
        for(set<GATKFeature>::iterator rodIter = rods_it->second.begin(); rodIter != rods_it->second.end();rodIter++) {
            const GATKFeature & feature = *rodIter;
            const GATKFeature & rod = feature;    // LCB TODO verify the correctness of this function vs original source code
            if ( indelRodsSeen.count(rod) > 0 )
                continue;
            indelRodsSeen.insert(rod);
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
    for (int readIndex = 0 ; readIndex < readSeq.size() ; refIndex++, readIndex++ ) {
        if ( refIndex >= refSeq.size() ) {
            sum += MAX_QUAL;
            // optimization: once we pass the threshold, stop calculating
            if ( sum > quitAboveThisValue )
                return sum;
        } else {
            char refChr = refSeq[refIndex];
            char readChr = readSeq[readIndex];
            if ( !BaseUtils::isRegularBase(readChr) || !BaseUtils::isRegularBase(refChr) )
                continue; // do not count Ns/Xs/etc ?
            if ( readChr != refChr ) {
                sum += (int)quals[readIndex] -33;   // Our qualities are still stored in ASCII
                // optimization: once we pass the threshold, stop calculating
                if ( sum > quitAboveThisValue )
                    return sum;
            }
        }
    }
    return sum;
}

void LocalRealignment::clean(ReadBin readsToClean, const GenomeLoc & currentInterval) {
    
    vector<BamAlignment *> reads = readsToClean.getReads();
    if ( reads.size() == 0 )
        return;
    
    string reference = readsToClean.getReference(*referenceReader);
    int leftmostIndex = readsToClean.getLocation().getStart();
    
    vector<BamAlignment *> refReads;                 // reads that perfectly match ref
    vector<AlignedRead *> altReads;               // reads that don't perfectly match
    vector<AlignedRead *> altAlignmentsToTest;  // should we try to make an alt consensus from the read?
    set<Consensus *> altConsenses;               // list of alt consenses
    
    // if there are any known indels for this region, get them and create alternate consenses
    generateAlternateConsensesFromKnownIndels(altConsenses, leftmostIndex, reference);

    // decide which reads potentially need to be cleaned;
    // if there are reads with a single indel in them, add that indel to the list of alternate consenses
    long totalRawMismatchSum = determineReadsThatNeedCleaning(reads, refReads, altReads, altAlignmentsToTest, altConsenses, leftmostIndex, reference);
    
    if ( verbose ) cerr << "------\nChecking consenses...(" << altConsenses.size() << ")\n--------\n" << endl;
    
    Consensus * bestConsensus = NULL;

    for (set<Consensus *>::iterator iter = altConsenses.begin(); iter != altConsenses.end(); iter++) {
        Consensus &consensus = **iter;
        if(verbose)
            cerr << "Trying new consensus: " << cigarToString( consensus.cigar) /*<< " " << consensus.str*/ << endl;
        
        consensus.cigar.size();
        if ( false ) {
            cerr << "Checking consensus with alignment at " << consensus.positionOnReference << " cigar " << cigarToString(consensus.cigar) << endl;
            //cerr << consensus.str << endl;
            int z = 0;
            for ( ; z < consensus.positionOnReference; z++ ) 
                cerr << ".";
            for ( z=0 ; z < consensus.cigar[0].Length ; z++ ) 
                cerr << ".";
            if ( consensus.cigar[1].Type == 'I' ) 
                for ( z= 0; z < consensus.cigar[1].Length; z++ )  
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
        
        if(verbose)
            cerr << "Mismatch sum of new consensus: " << consensus.mismatchSum << endl;
        if ( bestConsensus == NULL || bestConsensus->mismatchSum > consensus.mismatchSum) {
            // we do not need this alt consensus, release memory right away!!
            if ( bestConsensus != NULL ) {
                bestConsensus->readIndexes.clear();
            }
            bestConsensus = &consensus;
            if(verbose)
                cerr << "New consensus " << cigarToString(bestConsensus->cigar) <<  " is now best consensus" << endl;
        } else {
            // we do not need this alt consensus, release memory right away!!
            consensus.readIndexes.clear();
        }
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
            if ( write_out_stats ) {
                statsOutput << currentInterval.toString();
                statsOutput << "\tFAIL (bad indel)\t"; // if improvement > LOD_THRESHOLD *BUT* entropy is not reduced (SNPs still exist)
                statsOutput << improvement << endl;
            }
        } else {
            assert (bestConsensus->cigar.size() < 100);
            //cerr << "CLEAN: " << cigarToString(bestConsensus->cigar) << " " << string(bestConsensus->str) << " " << bestConsensus->cigar.size() << endl;

            if ( indelOutput != NULL && bestConsensus->cigar.size() > 1 ) {
                // NOTE: indels are printed out in the format specified for the low-coverage pilot1
                //  indel calls (tab-delimited): chr position size type sequence
                stringstream str;
                str << getHeader().Sequences[reads[0]->RefID].Name;
                int position = bestConsensus->positionOnReference + bestConsensus->cigar[0].Length;
                str << "\t" << (leftmostIndex + position - 1);
                CigarOp ce = bestConsensus->cigar[1];
                str << "\t"  << ce.Length << "\t" << ce.Type << "\t";
                int length = ce.Length;
                if ( ce.Type == 'D' ) {
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
                statsOutput << currentInterval.toString() << "\tCLEAN"; // if improvement > LOD_THRESHOLD *AND* entropy is reduced
                if ( bestConsensus->cigar.size() > 1 )
                    statsOutput << " (found indel)\t" << improvement << endl;
            }
            
            // finish cleaning the appropriate reads
            for ( vector<pair<int, int> >::iterator indexPair = bestConsensus->readIndexes.begin(); indexPair !=  bestConsensus->readIndexes.end(); indexPair++ ) {
                AlignedRead & aRead = *altReads[indexPair->first];
                if ( aRead.constizeUpdate() ) {
                    // We need to update the mapping quality score of the cleaned reads;
                    // however we don't have enough info to use the proper MAQ scoring system.
                    // For now, we will just arbitrarily add 10 to the mapping quality. [EB, 6/7/2010].
                    // TODO -- we need a better solution here
                    BamAlignment * read = aRead.getRead();  //TODO LCB verify the below mapping quality, since BamTools and SAM format are a bit different.
                    if ( read->MapQuality != 255 ) // 255 == Unknown, so don't modify it
                        read->MapQuality = min(aRead.getRead()->MapQuality + 10, 254);
                    
                    // before we fix the attribute tags we first need to make sure we have enough of the reference sequence
                    int neededBasesToLeft = leftmostIndex - read->Position;
                    int neededBasesToRight = read->Position + read->Length - leftmostIndex - reference.size() + 1;
                    int neededBases = max(neededBasesToLeft, neededBasesToRight);
                    if ( neededBases > 0 ) {
                        int padLeft = max(leftmostIndex-neededBases, 1);
                        int padRight = min(leftmostIndex+reference.size()+neededBases, (unsigned long)atol(referenceReader->getSequenceDictionary()[currentInterval.getContig()].Length.c_str()));
                        reference = referenceReader->getSubsequenceAt(currentInterval.getContig(), padLeft, padRight);
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
                    readsActuallyCleaned.insert(read);
                }
            }
        }
        
        // END IF ( improvement >= LOD_THRESHOLD )
        
    } else if ( write_out_stats ) {
        statsOutput << currentInterval.toString() << "\tFAIL\t" << improvement << endl;
    }
}

void LocalRealignment::generateAlternateConsensesFromKnownIndels(set<Consensus *> & altConsensesToPopulate, const int leftmostIndex, const string reference) {
    for ( vector<VariantContext>::iterator knownIndelIt  = knownIndelsToTry.begin(); knownIndelIt != knownIndelsToTry.end(); knownIndelIt++ ) {
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

long LocalRealignment::determineReadsThatNeedCleaning( vector<BamAlignment *> & reads,
                                            vector<BamAlignment *> & refReadsToPopulate,
                                            vector<AlignedRead *> & altReadsToPopulate,
                                            vector<AlignedRead *> & altAlignmentsToTest,
                                            set<Consensus *> & altConsenses,
                                            int leftmostIndex,
                                            string & reference) {
    
    long totalRawMismatchSum = 0L;
    
    for ( int read_ctr = 0; read_ctr < reads.size(); read_ctr++ ) {
        BamAlignment * read = reads[read_ctr];
        
        // we can not deal with screwy records
        if ( read->CigarData.size() == 0 ) {
            refReadsToPopulate.push_back(read);
            continue;
        }
        
        AlignedRead * aRead = new AlignedRead(read);
        
        // first, move existing indels (for 1 indel reads only) to leftmost position within identical sequence
        int numBlocks = 0;
        for(vector<CigarOp>::iterator i = read->CigarData.begin(); i != read->CigarData.end(); i++)
            if(i->Type == 'M' || i->Type == '=' || i->Type == 'X')
                numBlocks++;

        if ( numBlocks == 2 ) {
            vector<CigarOp> newCigar = AlignmentUtils::leftAlignIndel(unclipCigar(read->CigarData), reference, read->QueryBases, read->Position-leftmostIndex, 0);
            aRead->setCigar(newCigar, false);
        }
        
        const int startOnRef = read->Position-leftmostIndex;
        const int rawMismatchScore = mismatchQualitySumIgnoreCigar(*aRead, reference, startOnRef, INT_MAX);
        
        // if this doesn't match perfectly to the reference, let's try to clean it
        if ( rawMismatchScore > 0 ) {
            if(verbose)
                cerr << "Adding " << read->Name << " with raw mismatch score " << rawMismatchScore << " to non-ref reads" << endl;
            
            if ( !read->IsDuplicate() )
                totalRawMismatchSum += rawMismatchScore;
            aRead->setMismatchScoreToReference(rawMismatchScore);
            aRead->setAlignerMismatchScore(AlignmentUtils::mismatchingQualities(aRead->getRead(), reference, startOnRef));
            
            // if it has an indel, let's see if that's the best consensus
            if ( consensusModel != KNOWNS_ONLY && numBlocks == 2 )  {
                Consensus * c = createAlternateConsensus(startOnRef, aRead->getCigar(), reference, aRead->getReadBases());
                if ( c != NULL ) {
                    altConsenses.insert(c);
                }
            } else {
                altAlignmentsToTest.push_back(aRead);
            }
            
            altReadsToPopulate.push_back(aRead);
        }
        // otherwise, we can emit it as is
        else {
            if(verbose)
                cerr << "Adding " << read->Name << " with raw mismatch score " << rawMismatchScore << " to ref reads" << endl;
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
LocalRealignment::Consensus * LocalRealignment::createAlternateConsensus(const int indexOnRef, const vector<CigarOp> & c, const string reference, const string readStr) {
    if ( indexOnRef < 0 )
        return NULL;
    
    // if there are no indels, we do not need this consensus, can abort early:
    if ( c.size() == 1 && c[0].Type == 'M' ) return NULL;
    
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
        int elementLength = ce.Length;
        switch( ce.Type ) {
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
LocalRealignment::Consensus * LocalRealignment::createAlternateConsensus(const int indexOnRef, const string & reference, const string & indelStr, VariantContext indel) {
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


pair<int, int> LocalRealignment::findBestOffset(const string & ref, AlignedRead read, const int leftmostIndex) {
    
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
/*
string LocalRealignment::cigarToString(const vector<CigarOp> & cigar)
{
    stringstream ss;
    for(vector<CigarOp>::const_iterator i = cigar.begin(); i != cigar.end(); i++)
        ss << i->Length << i->Type;
    
    return ss.str();
}
 */

bool LocalRealignment::updateRead(const vector<CigarOp> & altCigar, const int altPosOnRef, const int myPosOnAlt, AlignedRead & aRead, const int leftmostIndex) {
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
    if ( altCE1.Type == 'I'  ) {
        indelCE=altCE1;
        if ( altCE2.Type != 'M' ) {
            cerr << "When the first element of the alt consensus is I, the second one must be M. Actual: " << cigarToString(altCigar) << ".  Skipping this site..." << endl;
            return false;
        }
    }
    else {
        if ( altCE1.Type != 'M'  ) {
            cerr << "First element of the alt consensus cigar must be M or I. Actual: " << cigarToString(altCigar) << ".  Skipping this site..." << endl;
            return false;
        }
        if ( altCE2.Type == 'I'  || altCE2.Type == 'D' ) {
            indelCE=altCE2;
        } else {
            cerr << "When first element of the alt consensus is M, the second one must be I or D. Actual: " << cigarToString(altCigar) << ".  Skipping this site..." << endl;
            return false;
        }
        leadingMatchingBlockLength = altCE1.Length;
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
    if ( indelCE.Type == 'I' ) {
        // for reads that end in an insertion
        if ( myPosOnAlt + aRead.getReadLength() < endOfFirstBlock + indelCE.Length ) {
            int partialInsertionLength = myPosOnAlt + aRead.getReadLength() - endOfFirstBlock;
            // if we also started inside the insertion, then we need to modify the length
            if ( !sawAlignmentStart )
                partialInsertionLength = aRead.getReadLength();
            readCigar.push_back( CigarOp('I', partialInsertionLength));
            aRead.setCigar(readCigar);
            return true;
        }
        
        // for reads that start in an insertion
        if ( !sawAlignmentStart && myPosOnAlt < endOfFirstBlock + indelCE.Length ) {
            aRead.setAlignmentStart(leftmostIndex + endOfFirstBlock);
            readCigar.push_back(CigarOp('I', indelCE.Length - (myPosOnAlt - endOfFirstBlock)));
            //indelOffsetOnRead = myPosOnAlt - endOfFirstBlock;
            sawAlignmentStart = true;
        } else if ( sawAlignmentStart ) {
            readCigar.push_back(indelCE);
            //indelOffsetOnRead = indelCE.getLength();
        }
    } else if ( indelCE.Type == 'D' ) {
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
        if ( ce->Type != 'D' )
            readRemaining -= ce->Length;
    }
    if ( readRemaining > 0 )
        readCigar.push_back(CigarOp('M', readRemaining));
    aRead.setCigar(readCigar);
    
    return true;
}

bool LocalRealignment::alternateReducesEntropy(vector<AlignedRead *> & reads, const string & reference, const int leftmostIndex) {
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
        for(vector<CigarOp>::iterator i = read.getRead()->CigarData.begin(); i != read.getRead()->CigarData.end(); i++)
            if(i->Type == 'M' || i->Type == '=' || i->Type == 'X')
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
            totalOriginalBases[refIdx] += quals[j];
            if ( readStr[j] != reference[refIdx] )
                originalMismatchBases[refIdx] += quals[j];
        }
        
        // reset and now do the calculation based on the cleaning
        refIdx = read.getAlignmentStart() - leftmostIndex;
        int altIdx = 0;
        vector<CigarOp> c = read.getCigar();
        for (int j = 0 ; j < c.size() ; j++) {
            CigarOp ce = c[j];
            int elementLength = ce.Length;
            switch ( ce.Type ) {
                case 'M':
                    for (int k = 0 ; k < elementLength ; k++, refIdx++, altIdx++ ) {
                        if ( refIdx >= reference.size() )
                            break;
                        totalCleanedBases[refIdx] += quals[altIdx];
                        if ( readStr[altIdx] != reference[refIdx] )
                            cleanedMismatchBases[refIdx] += quals[altIdx];
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
        if ( snpsOutput != NULL ) {
            if ( didMismatch ) {
                const string & ref_name = getHeader().Sequences[reads[0]->getRead()->RefID].Name;
                sb << ref_name << ":" << (leftmostIndex + i);
                if ( stillMismatches )
                    sb << " SAME_SNP\n" << endl;
                else
                    sb << " NOT_SNP" << endl;
            }
        }
    }
    
    if(verbose)
        cerr << "Original mismatch columns = " << originalMismatchColumns << "; cleaned mismatch columns = " << cleanedMismatchColumns << endl;
    
    const bool reduces = (originalMismatchColumns == 0 || cleanedMismatchColumns < originalMismatchColumns);
    if ( reduces && snpsOutput != NULL ) {
        snpsOutput << sb.str();
    }
    return reduces;
}
    
vector<CigarOp> LocalRealignment::unclipCigar(const vector<CigarOp> & cigar) {
    vector<CigarOp> elements;
    elements.reserve(cigar.size());

    for (int i = 0; i < cigar.size(); i++ ) {
        CigarOp ce = cigar[i];
        if ( !isClipOperator(ce.Type) )
            elements.push_back(ce);
    }
    return elements;
}

bool LocalRealignment::isClipOperator(const CigarOp op) {
    return op.Type == 'S' || op.Type == 'H' || op.Type == 'P';
}

vector<CigarOp> LocalRealignment::reclipCigar(const vector<CigarOp> & cigar, BamAlignment * read) {
    vector<CigarOp> elements;
    
    int i = 0;
    int n = read->CigarData.size();
    while ( i < n && isClipOperator(read->CigarData[i].Type) )
        elements.push_back(read->CigarData[i++]);
    
    //add all elements of cigar to elements
    elements.insert(elements.end(), cigar.begin(), cigar.end());    
    
    i++;
    while ( i < n && !isClipOperator(read->CigarData[i].Type) )
        i++;
    
    while ( i < n && isClipOperator(read->CigarData[i].Type) )
        elements.push_back(read->CigarData[i++]);
    
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
, write_out_indels(true)
, write_out_stats(true)
, write_out_snps(true)
, referenceReader(NULL)
, loc_parser(NULL)
, sawReadInCurrentInterval(false)
, outputIndels(false)
, output_stats(false)
, output_snps(false)
{}

int LocalRealignment::runInternal()
{
    initialize();

    GenomeLoc * currentInterval = intervalsFile.empty() ? NULL : intervalsFile.front();
    
    readsToClean.initialize(loc_parser, getHeader());
    
    while(true)
    {
        BamAlignment * al = getInputAlignment();
        
        if(!al)
            break;
        
        const ReadMetaDataTracker rmdt(loc_parser, al, std::map<int, RODMetaDataContainer>() );
        map_func(al, rmdt, currentInterval);
    };
    
    onTraversalDone(0);
    
    return 0;
}