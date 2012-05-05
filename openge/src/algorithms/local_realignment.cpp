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
#include "api/BamReader.h"
#include "api/BamWriter.h"
using namespace BamTools;

//still have to port these classes
class GenomeLoc;
class IntervalBinding;
class GATKFeature;
class VariantContext;
class ReferenceContext;
class ReadMetaDataTracker;
class ConstrainedMateFixingManager;
class IndexedFastaSequenceFile;
template <class T> class RodBinding<T>;

#include <list>
#include <set>
#include <vector>
#include <string>
#include <sstream>
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
 *     <li>Running the realigner over those intervals (IndelRealigner)</li>
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
 *   -T IndelRealigner \
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
static const string PROGRAM_RECORD_NAME("OpenGE IndelRealigner");

class IndelRealigner //extends ReadWalker<int, int> 
{
public:
    
    typedef enum {
        /**
         * Uses only indels from a provided ROD of known indels.
         */
        KNOWNS_ONLY,
        /**
         * Additionally uses indels already present in the original alignments of the reads.
         */
        USE_READS,
        /**
         * Additionally uses 'Smith-Waterman' to generate alternate consenses.
         */
        USE_SW
    } ConsensusDeterminationModel_t;
    
    /**
     * Any number of VCF files representing known indels to be used for constructing alternate consenses.
     * Could be e.g. dbSNP and/or official 1000 Genomes indel calls.  Non-indel variants in these files will be ignored.
     */
    //@Input(fullName="knownAlleles", shortName = "known", doc="Input VCF file(s) with known indels", required=false)
    list<RodBinding<VariantContext> > known;
    
protected:
    /**
     * The interval list output from the RealignerTargetCreator tool using the same bam(s), reference, and known indel file(s).
     */
    //@Input(fullName="targetIntervals", shortName="targetIntervals", doc="intervals file output from RealignerTargetCreator", required=true)
    IntervalBinding<Feature *> intervalsFile = NULL;
    
    /**
     * This term is equivalent to "significance" - i.e. is the improvement significant enough to merit realignment? Note that this number
     * should be adjusted based on your particular data set. For low coverage and/or when looking for indels with low allele frequency,
     * this number should be smaller.
     */
    //@Argument(fullName="LODThresholdForCleaning", shortName="LOD", doc="LOD threshold above which the cleaner will clean", required=false)
    double LOD_THRESHOLD = 5.0;
    
    /**
     * The realigned bam file.
     */
    //@Output(required=false, doc="Output bam")
    BamWriter writer;
    ConstrainedMateFixingManager manager;
    BamWriter writerToUse;

public:
    /**
     * We recommend that users run with USE_READS when trying to realign high quality longer read data mapped with a gapped aligner;
     * Smith-Waterman is really only necessary when using an ungapped aligner (e.g. MAQ in the case of single-end read data).
     */
    //@Argument(fullName = "consensusDeterminationModel", shortName = "model", doc = "Determines how to compute the possible alternate consenses", required = false)
    ConsensusDeterminationModel_t consensusModel = USE_READS;
    
    
    // ADVANCED OPTIONS FOLLOW
    
    /**
     * For expert users only!  This is similar to the argument in the RealignerTargetCreator walker. The point here is that the realigner
     * will only proceed with the realignment (even above the given threshold) if it minimizes entropy among the reads (and doesn't simply
     * push the mismatch column to another position). This parameter is just a heuristic and should be adjusted based on your particular data set.
     */
protected:
    //@Argument(fullName="entropyThreshold", shortName="entropy", doc="percentage of mismatches at a locus to be considered having high entropy", required=false)
    double MISMATCH_THRESHOLD = 0.15;
    
    /**
     * For expert users only!  To minimize memory consumption you can lower this number (but then the tool may skip realignment on regions with too much coverage;
     * and if the number is too low, it may generate errors during realignment). Just make sure to give Java enough memory! 4Gb should be enough with the default value.
     */
    //@Advanced
    //@Argument(fullName="maxReadsInMemory", shortName="maxInMemory", doc="max reads allowed to be kept in memory at a time by the SAMFileWriter", required=false)
    int MAX_RECORDS_IN_MEMORY = 150000;
    
    /**
     * For expert users only!
     */
    //@Advanced
    //@Argument(fullName="maxIsizeForMovement", shortName="maxIsize", doc="maximum insert size of read pairs that we attempt to realign", required=false)
    int MAX_ISIZE_FOR_MOVEMENT = 3000;
    
    /**
     * For expert users only!
     */
    //@Advanced
    //@Argument(fullName="maxPositionalMoveAllowed", shortName="maxPosMove", doc="maximum positional move in basepairs that a read can be adjusted during realignment", required=false)
    static int MAX_POS_MOVE_ALLOWED = 200;
    
    /**
     * For expert users only!  If you need to find the optimal solution regardless of running time, use a higher number.
     */
    //@Advanced
    //@Argument(fullName="maxConsensuses", shortName="maxConsensuses", doc="max alternate consensuses to try (necessary to improve performance in deep coverage)", required=false)
    int MAX_CONSENSUSES = 30;
    
    /**
     * For expert users only!  If you need to find the optimal solution regardless of running time, use a higher number.
     */
    //@Advanced
    //@Argument(fullName="maxReadsForConsensuses", shortName="greedy", doc="max reads used for finding the alternate consensuses (necessary to improve performance in deep coverage)", required=false)
    int MAX_READS_FOR_CONSENSUSES = 120;
    
    /**
     * For expert users only!  If this value is exceeded at a given interval, realignment is not attempted and the reads are passed to the output file(s) as-is.
     * If you need to allow more reads (e.g. with very deep coverage) regardless of memory, use a higher number.
     */
    //@Advanced
    //@Argument(fullName="maxReadsForRealignment", shortName="maxReads", doc="max reads allowed at an interval for realignment", required=false)
    int MAX_READS = 20000;
    
    //@Advanced
    //@Argument(fullName="noOriginalAlignmentTags", shortName="noTags", required=false, doc="Don't output the original cigar or alignment start tags for each realigned read in the output bam")
    static bool NO_ORIGINAL_ALIGNMENT_TAGS = false;
    
    /**
     * Reads from all input files will be realigned together, but then each read will be saved in the output file corresponding to the input file that
     * the read came from. There are two ways to generate output bam file names: 1) if the value of this argument is a general string (e.g. '.cleaned.bam'),
     * then extensions (".bam" or ".sam") will be stripped from the input file names and the provided string value will be pasted on instead; 2) if the
     * value ends with a '.map' (e.g. input_output.map), then the two-column tab-separated file with the specified name must exist and list unique output
     * file name (2nd column) for each input file name (1st column).
     */
    //@Argument(fullName="nWayOut", shortName="nWayOut", required=false, doc="Generate one output file for each input (-I) bam file")
    string N_WAY_OUT;
    
    //@Hidden
    //@Argument(fullName="generate_nWayOut_md5s",doc="Generate md5sums for BAMs")
    bool generateMD5s = false;
    
    // DEBUGGING OPTIONS FOLLOW
    
    //@Hidden
    //@Argument(fullName="check_early",shortName="check_early",required=false,doc="Do early check of reads against existing consensuses")
    bool CHECKEARLY = false;
    
    //@Hidden
    //@Argument(fullName="noPGTag", shortName="noPG", required=false, doc="Don't output the usual PG tag in the realigned bam file header. FOR DEBUGGING PURPOSES ONLY.  This option is required in order to pass integration tests.")
    bool NO_PG_TAG = false;
    
    //@Hidden
    //@Argument(fullName="keepPGTags", shortName="keepPG", required=false, doc="Keep older PG tags left in the bam header by previous runs of this tool (by default, all these "+ "historical tags will be replaced by the latest tag generated in the current run).")
    bool KEEP_ALL_PG_RECORDS = false;
    
    //@Hidden
    //@Output(fullName="indelsFileForDebugging", shortName="indels", required=false, doc="Output file (text) for the indels found; FOR DEBUGGING PURPOSES ONLY")
    string OUT_INDELS = NULL;
    
    // @Hidden
    //@Output(fullName="statisticsFileForDebugging", shortName="stats", doc="print out statistics (what does or doesn't get cleaned); FOR DEBUGGING PURPOSES ONLY", required=false)
    string OUT_STATS = NULL;
    
    //@Hidden
    //@Output(fullName="SNPsFileForDebugging", shortName="snps", doc="print out whether mismatching columns do or don't get cleaned out; FOR DEBUGGING PURPOSES ONLY", required=false)
    string OUT_SNPS = NULL;
    
private:
    static vector<CigarOp> reclipCigar(vector<CigarOp> & cigar, BamAlignment * read);
    
    class AlignedRead {
    private:
        BamAlignment * read;
        string * readBases;
        string * baseQuals;
        vector<CigarOp> * newCigar = NULL;
        int newStart = -1;
        int mismatchScoreToReference = 0;
        long alignerMismatchScore = 0;
    public:
        AlignedRead(BamAlignment * read) 
        : read(read) 
        , mismatchScoreToReference(0)
        { }
        
        BamAlignment * getRead() {
            return read;
        }
        
        int getReadLength() {
            return readBases != NULL ? readBases->size() : read->Length;
        }
        
        string getReadBases() {
            if ( readBases == NULL )
                getUnclippedBases();
            return *readBases;
        }
        
        string getBaseQualities() {
            if ( baseQuals == NULL )
                getUnclippedBases();
            return *baseQuals;
        }
        
    private:
        // pull out the bases that aren't clipped out
        void getUnclippedBases() {
            readBases = new string('?', getReadLength());
            baseQuals = new string('?', getReadLength());
            const string actualReadBases = read->QueryBases;// LCB should be AlignedBases?
            const string actualBaseQuals = read->Qualities;
            int fromIndex = 0, toIndex = 0;
            
            for ( vector<CigarOp>::iterator ce = read->CigarData.begin(); ce != read->CigarData.end(); ce++ ) {
                uint32_t & elementLength = ce->Length;
                switch ( ce->Type ) {
                    case 'S':
                        fromIndex += elementLength;
                        break;
                    case 'M':
                    case 'I':
                        readBases->replace( fromIndex, elementLength, string(actualReadBases, fromIndex, elementLength));
                        baseQuals->replace( fromIndex, elementLength, string(actualBaseQuals, fromIndex, elementLength));
                        //System.arraycopy(actualReadBases, fromIndex, readBases, toIndex, elementLength);
                        //System.arraycopy(actualBaseQuals, fromIndex, baseQuals, toIndex, elementLength);
                        fromIndex += elementLength;
                        toIndex += elementLength;
                    default:
                        break;
                }
            }
            
            // if we got clipped, trim the array
            if ( fromIndex != toIndex ) {
                readBases->erase(toIndex);
                baseQuals->erase(toIndex);
            }
        }

        // pull out the bases that aren't clipped out
        vector<CigarOp> reclipCigar(vector<CigarOp> & cigar) {
            return IndelRealigner::reclipCigar(cigar, read);
        }
    public:
        vector<CigarOp> & getCigar() {
            return (newCigar != NULL ? *newCigar : read->CigarData);
        }
        
        // tentatively sets the new Cigar, but it needs to be confirmed later
        void setCigar(vector<CigarOp> * cigar, bool fixClippedCigar = true) {
            if ( cigar == NULL ) {
                delete newCigar;
                newCigar = NULL;
                return;
            }
            
            if ( fixClippedCigar && getReadBases().length() < read->Length ) {
                cigar = new vector<CigarOp>(reclipCigar(*cigar));
            }
            
            // no change?
            if ( read->CigarData == *cigar ) {
                newCigar = NULL;
                return;
            }

            // no indel?
            bool has_i_or_d = false;
            for (vector<CigarOp>::iterator i = cigar->begin(); i != cigar->end(); i++)
                if(i->Type == 'I' || i->Type == 'D')
                    has_i_or_d = true;

            if ( !has_i_or_d) {
                cerr << "Modifying a read with no associated indel; although this is possible, it is highly unlikely.  Perhaps this region should be double-checked: " << read->Name << " near read->getReferenceName():" << read->Position << endl; // FIXME LCB, can't easily access the reference name with current code organization. We should do something about this.
                //    newCigar = NULL;
                //    return;
            }
            
            newCigar = cigar;
        }
        
    public:
        // tentatively sets the new start, but it needs to be confirmed later
        void setAlignmentStart(int start) {
            newStart = start;
        }
        
        int getAlignmentStart() {
            return (newStart != -1 ? newStart : read->Position);
        }
        
        int getOriginalAlignmentStart() {
            return read->Position;
        }
        
        // constizes the changes made.
        // returns true if this record actually changes, false otherwise
        bool constizeUpdate() {
            // if we haven't made any changes, don't do anything
            if ( newCigar == NULL )
                return false;
            if ( newStart == -1 )
                newStart = read->Position;
            else if ( abs(newStart - read->Position) > MAX_POS_MOVE_ALLOWED ) {
                cerr << "Attempting to realign read " << read->Name << " at " << read->Position << " more than " << MAX_POS_MOVE_ALLOWED << " bases to " << newStart << ".";
                return false;
            }
            
            // annotate the record with the original cigar (and optionally the alignment start)
            if ( !NO_ORIGINAL_ALIGNMENT_TAGS ) {
                stringstream cigar_ss("");
                for(vector<CigarOp>::iterator i = read->CigarData.begin(); i != read->CigarData.end(); i++)
                    cigar_ss << i->Length << i->Type;
                read->AddTag(ORIGINAL_CIGAR_TAG, "Z", cigar_ss.str());
                if ( newStart != read->Position )
                    read->AddTag(ORIGINAL_POSITION_TAG, "i", read->Position);
            }
            
            read->CigarData = *newCigar;
            read->Position = newStart;
            
            return true;
        }
        
        void setMismatchScoreToReference(int score) {
            mismatchScoreToReference = score;
        }
        
        int getMismatchScoreToReference() {
            return mismatchScoreToReference;
        }
        
        void setAlignerMismatchScore(long score) {
            alignerMismatchScore = score;
        }
        
        long getAlignerMismatchScore() {
            return alignerMismatchScore;
        }
    };
    
    class Consensus {
        public:
        const string str;
        const vector<pair<int, int> > readIndexes;
        const int positionOnReference;
        int mismatchSum;
        vector<CigarOp> cigar;
        
        Consensus(string str, vector<CigarOp> cigar, int positionOnReference) 
        : str(str)
        , cigar(cigar)
        , positionOnReference(positionOnReference)
        , mismatchSum(0)
        {}
        
        /*
        public bool equals(Object o) {
            return ( this == o || (o instanceof Consensus && Arrays.equals(this.str,(((Consensus)o).str)) ) );
        } */
        
        bool operator==(const Consensus & c) {
            return ( this == &c || str == c.str ) ;
        }
    };
    
    class ReadBin {// implements HasGenomeLocation {
    private:
        const vector<BamAlignment *> reads;
        string * reference = NULL;
        GenomeLoc * loc = NULL;
    public:
        ReadBin() { }
        
        // Return false if we can't process this read bin because the reads are not correctly overlapping.
        // This can happen if e.g. there's a large known indel with no overlapping reads.
        void add(BamAlignment *) {
            
            GenomeLoc locForRead = getToolkit().getGenomeLocParser().createGenomeLoc(read);
            if ( loc == NULL )
                loc = locForRead;
            else if ( locForRead->getStop() > loc->getStop() )
                loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), loc.getStart(), locForRead.getStop());
            
            reads.push_back(read);
        }
        
        vector<BamAlignment *> getReads() { return reads; }
        
        string getReference(IndexedFastaSequenceFile & referenceReader) {
            // set up the reference if we haven't done so yet
            if ( reference == NULL ) {
                // first, pad the reference to handle deletions in narrow windows (e.g. those with only 1 read)
                int padLeft = max(loc.getStart()-REFERENCE_PADDING, 1);
                int padRight = min(loc.getStop()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(loc.getContig()).getSequenceLength());
                loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), padLeft, padRight);
                reference = referenceReader.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
                StringUtil.toUpperCase(reference);
            }
            
            return reference;
        }
        
        GenomeLoc getLocation() { return loc; }
        
        int size() { return reads.size(); }
        
        void clear() {
            reads.clear();
            reference = NULL;
            loc = NULL;
        }
        
    };

private:
    // fasta reference reader to supplement the edges of the reference sequence
    IndexedFastaSequenceFile referenceReader;
    
    // the intervals input by the user
    vector<GenomeLoc>::iterator * intervals = NULL;
    
    // the current interval in the list
    GenomeLoc * currentInterval = NULL;
    bool sawReadInCurrentInterval = false;
    
    // the reads and known indels that fall into the current interval
    ReadBin readsToClean;
    vector<BamAlignment *> readsNotToClean;
    vector<VariantContext *> knownIndelsToTry;
    set<GATKFeature> indelRodsSeen;
    set<BamAlignment *> readsActuallyCleaned;
    
    static const int MAX_QUAL = 99;
    
    // fraction of mismatches that need to no longer mismatch for a column to be considered cleaned
    static const double MISMATCH_COLUMN_CLEANED_FRACTION = 0.75;
    
    static const double SW_MATCH = 30.0;      // 1.0;
    static const double SW_MISMATCH = -10.0;  //-1.0/3.0;
    static const double SW_GAP = -10.0;       //-1.0-1.0/3.0;
    static const double SW_GAP_EXTEND = -2.0; //-1.0/.0;
    
    // reference base padding size
    // TODO -- make this a command-line argument if the need arises
    static const int REFERENCE_PADDING = 30;
    
    // other output files
    bool outputIndels = true, output_stats = true, output_snps = true;
    ostream indelOutput, statsOutput, snpsOutput;
    
    //###protected Map<SAMReaderID, ConstrainedMateFixingManager> nwayWriters = NULL;
    
    
    // debug info for lazy SW evaluation:
    long exactMatchesFound = 0; // how many reads exactly matched a consensus we already had
    long SWalignmentRuns = 0; // how many times (=for how many reads) we ran SW alignment
    long SWalignmentSuccess = 0; // how many SW alignments were "successful" (i.e. found a workable indel and resulted in non-null consensus)
    
    map<string,string> loadFileNameMap(string mapFile) {
        Map<String,String> fname_map = new HashMap<String,String>();
        
        try {
            
            XReadLines reader = new XReadLines(new File(mapFile),true);
            for ( String line : reader ) {
                if ( line.length() == 0 ) continue;
                
                String fields[] = line.split("\t");
                
                if ( fields.length != 2 )
                    throw new UserException.BadInput("Input-output map file must have exactly two columns. Offending line:\n"+line);
                if ( fields[0].length() == 0 || fields[1].length() == 0 )
                    throw new UserException.BadInput("Input-output map file can not have empty strings in either column. Offending line:\n"+line);
                
                if ( fname_map.containsKey(fields[0]) )
                    throw new UserException.BadInput("Input-output map file contains duplicate entries for input name "+fields[0]);
                if ( fname_map.containsValue(fields[1]) )
                    throw new UserException.BadInput("Input-output map file maps multiple entries onto single output name "+fields[1]);
                
                fname_map.put(fields[0],fields[1]);
            }
        } catch (IOException e) {
            throw new StingException("I/O Error while reading input-output map file "+N_WAY_OUT+": "+e.getMessage());
        }
        return fname_map;
    }
public:
    void initialize() {
        
        if ( N_WAY_OUT == NULL && writer == NULL ) {
            throw new UserException.CommandLineException("Either -o or -nWayOut must be specified");
        }
        if ( N_WAY_OUT != NULL && writer != NULL ) {
            throw new UserException.CommandLineException("-o and -nWayOut can not be used simultaneously");
        }
        if ( LOD_THRESHOLD < 0.0 )
            throw new RuntimeException("LOD threshold cannot be a negative number");
        if ( MISMATCH_THRESHOLD <= 0.0 || MISMATCH_THRESHOLD > 1.0 )
            throw new RuntimeException("Entropy threshold must be a fraction between 0 and 1");
        
        try {
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile,ex);
        }
        
        intervals = intervalsFile.getIntervals(getToolkit()).iterator();
        
        currentInterval = intervals.hasNext() ? intervals.next() : NULL;
        
        writerToUse = writer;
        
        if ( N_WAY_OUT != NULL ) {
            boolean createIndex =  true;
            
            if ( N_WAY_OUT.toUpperCase().endsWith(".MAP") ) {
                writerToUse = new NWaySAMFileWriter(getToolkit(),loadFileNameMap(N_WAY_OUT),
                                                    SAMFileHeader.SortOrder.coordinate,true, createIndex, generateMD5s,createProgramRecord(),KEEP_ALL_PG_RECORDS);
            } else {
                writerToUse = new NWaySAMFileWriter(getToolkit(),N_WAY_OUT,SAMFileHeader.SortOrder.coordinate,true,
                                                    createIndex, generateMD5s,createProgramRecord(),KEEP_ALL_PG_RECORDS);
            }
        }   else {
            
            // set up the output writer
            setupWriter(getToolkit().getSAMFileHeader());
        }
        manager = new ConstrainedMateFixingManager(writerToUse, getToolkit().getGenomeLocParser(), MAX_ISIZE_FOR_MOVEMENT, MAX_POS_MOVE_ALLOWED, MAX_RECORDS_IN_MEMORY);
        
        if ( OUT_INDELS != NULL ) {
            try {
                indelOutput = new FileWriter(new File(OUT_INDELS));
            } catch (Exception e) {
                logger.error("Failed to create output file "+ OUT_INDELS+". Indel output will be suppressed");
                logger.error(e.getMessage());
                indelOutput = NULL;
            }
        }
        if ( OUT_STATS != NULL ) {
            try {
                statsOutput = new FileWriter(new File(OUT_STATS));
            } catch (Exception e) {
                logger.error("Failed to create output file "+ OUT_STATS+". Cleaning stats output will be suppressed");
                logger.error(e.getMessage());
                statsOutput = NULL;
            }
        }
        if ( OUT_SNPS != NULL ) {
            try {
                snpsOutput = new FileWriter(new File(OUT_SNPS));
            } catch (Exception e) {
                logger.error("Failed to create output file "+ OUT_SNPS+". Cleaning snps output will be suppressed");
                logger.error(e.getMessage());
                snpsOutput = NULL;
            }
        }
    }

private:
    void setupWriter(const SamHeader & header) {
        
        if ( !NO_PG_TAG ) {
            const SAMProgramRecord programRecord = createProgramRecord();
            
            List<SAMProgramRecord> oldRecords = header.getProgramRecords();
            List<SAMProgramRecord> newRecords = new vector<SAMProgramRecord>(oldRecords.size()+1);
            for ( SAMProgramRecord record : oldRecords ) {
                if ( !record.getId().startsWith(PROGRAM_RECORD_NAME) || KEEP_ALL_PG_RECORDS )
                    newRecords.add(record);
            }
            newRecords.add(programRecord);
            header.setProgramRecords(newRecords);
        }
        
        writer.writeHeader(header);
        writer.setPresorted(true);
    }
    
    
    SamProgram createProgramRecord() {
        if ( NO_PG_TAG ) return NULL;
        
        const SAMProgramRecord programRecord = new SAMProgramRecord(PROGRAM_RECORD_NAME);
        const ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
        try {
            const String version = headerInfo.getString("org.broadinstitute.sting.gatk.version");
            programRecord.setProgramVersion(version);
        } catch (MissingResourceException e) {}
        programRecord.setCommandLine(getToolkit().createApproximateCommandLineArgumentString(getToolkit(), this));
        return programRecord;
    }
    
    void emit(const BamAlignment * read) {
        
        // check to see whether the read was modified by looking at the temporary tag
        boolean wasModified = readsActuallyCleaned.contains(read);
        
        try {
            manager.addRead(read, wasModified);
        } catch (RuntimeIOException e) {
            throw new UserException.ErrorWritingBamFile(e.getMessage());
        }
    }
    
    void emitReadLists() {
        // pre-merge lists to sort them in preparation for constrained SAMFileWriter
        readsNotToClean.addAll(readsToClean.getReads());
        ReadUtils.sortReadsByCoordinate(readsNotToClean);
        manager.addReads(readsNotToClean, readsActuallyCleaned);
        readsToClean.clear();
        readsNotToClean.clear();
        readsActuallyCleaned.clear();
    }
public:
    int map(ReferenceContext ref, BamAlignment * read, ReadMetaDataTracker metaDataTracker) {
        if ( currentInterval == NULL ) {
            emit(read);
            return 0;
        }
        
        // edge case: when the last target interval abuts the end of the genome, we'll get one of the
        //   unmapped reads while the currentInterval still isn't null.  We need to trigger the cleaning
        //   at this point without trying to create a GenomeLoc.
        if ( read.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX ) {
            cleanAndCallMap(ref, read, metaDataTracker, NULL);
            return 0;
        }
        
        GenomeLoc readLoc = getToolkit().getGenomeLocParser().createGenomeLoc(read);
        // hack to get around unmapped reads having screwy locations
        if ( readLoc.getStop() == 0 )
            readLoc = getToolkit().getGenomeLocParser().createGenomeLoc(readLoc.getContig(), readLoc.getStart(), readLoc.getStart());
        
        if ( readLoc.isBefore(currentInterval) ) {
            if ( !sawReadInCurrentInterval )
                emit(read);
            else
                readsNotToClean.add(read);
        }
        else if ( readLoc.overlapsP(currentInterval) ) {
            sawReadInCurrentInterval = true;
            
            if ( doNotTryToClean(read) ) {
                readsNotToClean.add(read);
            } else {
                readsToClean.add(read);
                
                // add the rods to the list of known variants
                populateKnownIndels(metaDataTracker, ref);
            }
            
            if ( readsToClean.size() + readsNotToClean.size() >= MAX_READS ) {
                logger.info("Not attempting realignment in interval " + currentInterval + " because there are too many reads.");
                abortCleanForCurrentInterval();
            }
        }
        else {  // the read is past the current interval
            cleanAndCallMap(ref, read, metaDataTracker, readLoc);
        }
        
        return 0;
    }
    
private:
    void abortCleanForCurrentInterval() {
        emitReadLists();
        currentInterval = intervals.hasNext() ? intervals.next() : NULL;
        sawReadInCurrentInterval = false;
    }
    
    bool doNotTryToClean(const BamAlignment * read) {
        return read.getReadUnmappedFlag() ||
        read.getNotPrimaryAlignmentFlag() ||
        read.getReadFailsVendorQualityCheckFlag() ||
        read.getMappingQuality() == 0 ||
        read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START ||
        ConstrainedMateFixingManager.iSizeTooBigToMove(read, MAX_ISIZE_FOR_MOVEMENT) ||
        ReadUtils.is454Read(read);
        // TODO -- it would be nice if we could use indels from 454 reads as alternate consenses
    }
    
    void cleanAndCallMap(ReferenceContext ref, const BamAlignment * read, ReadMetaDataTracker metaDataTracker, GenomeLoc readLoc) {
        if ( readsToClean.size() > 0 ) {
            GenomeLoc earliestPossibleMove = getToolkit().getGenomeLocParser().createGenomeLoc(readsToClean.getReads().get(0));
            if ( manager.canMoveReads(earliestPossibleMove) )
                clean(readsToClean);
        }
        knownIndelsToTry.clear();
        indelRodsSeen.clear();
        
        emitReadLists();
        try {
            do {
                currentInterval = intervals.hasNext() ? intervals.next() : NULL;
                
            } while ( currentInterval != NULL && (readLoc == NULL || currentInterval.isBefore(readLoc)) );
        } catch (ReviewedStingException e) {
            throw new UserException.MissortedFile(new File(intervalsFile.getSource()), " *** Are you sure that your interval file is sorted? If not, you must use the --targetIntervalsAreNotSorted argument. ***", e);
        }
        sawReadInCurrentInterval = false;
        
        // call back into map now that the state has been updated
        map(ref, read, metaDataTracker);
    }

public:
    int reduceInit() {
        return 0;
    }
    
    int reduce(int value, int sum) {
        return sum + value;
    }
    
    void onTraversalDone(int result) {
        if ( readsToClean.size() > 0 ) {
            GenomeLoc earliestPossibleMove = getToolkit().getGenomeLocParser().createGenomeLoc(readsToClean.getReads().get(0));
            if ( manager.canMoveReads(earliestPossibleMove) )
                clean(readsToClean);
            emitReadLists();
        } else if ( readsNotToClean.size() > 0 ) {
            emitReadLists();                            
        }
        
        knownIndelsToTry.clear();
        indelRodsSeen.clear();
        
        if ( OUT_INDELS != NULL ) {
            try {
                indelOutput.close();
            } catch (Exception e) {
                logger.error("Failed to close "+OUT_INDELS+" gracefully. Data may be corrupt.");
            }
        }
        if ( OUT_STATS != NULL ) {
            try {
                statsOutput.close();
            } catch (Exception e) {
                logger.error("Failed to close "+OUT_STATS+" gracefully. Data may be corrupt.");
            }
        }
        if ( OUT_SNPS != NULL ) {
            try {
                snpsOutput.close();
            } catch (Exception e) {
                logger.error("Failed to close "+OUT_SNPS+" gracefully. Data may be corrupt.");
            }
        }
        
        manager.close();
        if ( N_WAY_OUT != NULL ) writerToUse.close();
        
        if ( CHECKEARLY ) {
            logger.info("SW alignments runs: "+SWalignmentRuns);
            logger.info("SW alignments successfull: "+SWalignmentSuccess + " ("+SWalignmentSuccess/SWalignmentRuns+"% of SW runs)");
            logger.info("SW alignments skipped (perfect match): "+exactMatchesFound);
            logger.info("Total reads SW worked for: "+(SWalignmentSuccess + exactMatchesFound)+
                        " ("+(SWalignmentSuccess+exactMatchesFound)/(SWalignmentRuns+exactMatchesFound)+"% of all reads requiring SW)");
        }
    }
    
private:
    void populateKnownIndels(ReadMetaDataTracker metaDataTracker, ReferenceContext ref) {
        for ( Collection<GATKFeature> rods : metaDataTracker.getContigOffsetMapping().values() ) {
            Iterator<GATKFeature> rodIter = rods.iterator();
            while ( rodIter.hasNext() ) {
                Object rod = rodIter.next().getUnderlyingObject();
                if ( indelRodsSeen.contains(rod) )
                    continue;
                indelRodsSeen.add(rod);
                if ( rod instanceof VariantContext )
                    knownIndelsToTry.add((VariantContext)rod);
            }
        }
    }
    
    static int mismatchQualitySumIgnoreCigar(const AlignedRead aRead, const string refSeq, int refIndex, int quitAboveThisValue) {
        const byte[] readSeq = aRead.getReadBases();
        const byte[] quals = aRead.getBaseQualities();
        int sum = 0;
        for (int readIndex = 0 ; readIndex < readSeq.length ; refIndex++, readIndex++ ) {
            if ( refIndex >= refSeq.length ) {
                sum += MAX_QUAL;
                // optimization: once we pass the threshold, stop calculating
                if ( sum > quitAboveThisValue )
                    return sum;
            } else {
                byte refChr = refSeq[refIndex];
                byte readChr = readSeq[readIndex];
                if ( !BaseUtils.isRegularBase(readChr) || !BaseUtils.isRegularBase(refChr) )
                    continue; // do not count Ns/Xs/etc ?
                if ( readChr != refChr ) {
                    sum += (int)quals[readIndex];
                    // optimization: once we pass the threshold, stop calculating
                    if ( sum > quitAboveThisValue )
                        return sum;
                }
            }
        }
        return sum;
    }
    
    void clean(ReadBin readsToClean) {
        
        const List<GATKSAMRecord> reads = readsToClean.getReads();
        if ( reads.size() == 0 )
            return;
        
        byte[] reference = readsToClean.getReference(referenceReader);
        int leftmostIndex = readsToClean.getLocation().getStart();
        
        const vector<GATKSAMRecord> refReads = new vector<GATKSAMRecord>();                 // reads that perfectly match ref
        const vector<AlignedRead> altReads = new vector<AlignedRead>();               // reads that don't perfectly match
        const LinkedList<AlignedRead> altAlignmentsToTest = new LinkedList<AlignedRead>();  // should we try to make an alt consensus from the read?
        const Set<Consensus> altConsenses = new LinkedHashSet<Consensus>();               // list of alt consenses
        
        // if there are any known indels for this region, get them and create alternate consenses
        generateAlternateConsensesFromKnownIndels(altConsenses, leftmostIndex, reference);
        
        // decide which reads potentially need to be cleaned;
        // if there are reads with a single indel in them, add that indel to the list of alternate consenses
        long totalRawMismatchSum = determineReadsThatNeedCleaning(reads, refReads, altReads, altAlignmentsToTest, altConsenses, leftmostIndex, reference);
        
        // use 'Smith-Waterman' to create alternate consenses from reads that mismatch the reference, using totalRawMismatchSum as the random seed
        if ( consensusModel == ConsensusDeterminationModel.USE_SW )
            generateAlternateConsensesFromReads(altAlignmentsToTest, altConsenses, reference, leftmostIndex);
        
        // if ( debugOn ) System.out.println("------\nChecking consenses...\n--------\n");
        
        Consensus bestConsensus = NULL;
        Iterator<Consensus> iter = altConsenses.iterator();
        
        while ( iter.hasNext() ) {
            Consensus consensus = iter.next();
            //logger.debug("Trying new consensus: " + consensus.cigar + " " + new String(consensus.str));
            
            //            if ( DEBUG ) {
            //                System.out.println("Checking consensus with alignment at "+consensus.positionOnReference+" cigar "+consensus.cigar);
            //                System.out.println(new String(consensus.str));
            //                int z = 0;
            //                for ( ; z < consensus.positionOnReference; z++ )  System.out.print('.');
            //                for ( z=0 ; z < consensus.cigar.getCigarElement(0).getLength() ; z++ ) System.out.print('.');
            //                if ( consensus.cigar.getCigarElement(1).getOperator() == CigarOperator.I ) for ( z= 0; z < consensus.cigar.getCigarElement(1).getLength(); z++ )  System.out.print('I');
            //                System.out.println();
            //            }
            
            // if ( debugOn ) System.out.println("Consensus: "+consensus.str);
            
            for ( int j = 0; j < altReads.size(); j++ ) {
                AlignedRead toTest = altReads.get(j);
                pair<int, int> altAlignment = findBestOffset(consensus.str, toTest, leftmostIndex);
                
                // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
                int myScore = altAlignment.second;
                
                if ( myScore > toTest.getAlignerMismatchScore() || myScore >= toTest.getMismatchScoreToReference() )
                    myScore = toTest.getMismatchScoreToReference();
                // keep track of reads that align better to the alternate consensus.
                // By pushing alignments with equal scores to the alternate, it means we'll over-call (het -> hom non ref) but are less likely to under-call (het -> ref, het non ref -> het)
                else
                    consensus.readIndexes.add(new pair<int, int>(j, altAlignment.first));
                
                //logger.debug(consensus.cigar +  " vs. " + toTest.getRead().getReadName() + "-" + toTest.getRead().getReadString() + " => " + myScore + " vs. " + toTest.getMismatchScoreToReference());
                if ( !toTest.getRead().getDuplicateReadFlag() )
                    consensus.mismatchSum += myScore;
                
                // optimization: once the mismatch sum is higher than the best consensus, quit since this one can't win
                //  THIS MUST BE DISABLED IF WE DECIDE TO ALLOW MORE THAN ONE ALTERNATE CONSENSUS!
                if ( bestConsensus != NULL && consensus.mismatchSum > bestConsensus.mismatchSum )
                    break;
            }
            
            //logger.debug("Mismatch sum of new consensus: " + consensus.mismatchSum);
            if ( bestConsensus == NULL || bestConsensus.mismatchSum > consensus.mismatchSum) {
                // we do not need this alt consensus, release memory right away!!
                if ( bestConsensus != NULL )
                    bestConsensus.readIndexes.clear();
                bestConsensus = consensus;
                //logger.debug("New consensus " + bestConsensus.cigar +  " is now best consensus");
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
        const double improvement = (bestConsensus == NULL ? -1 : ((double)(totalRawMismatchSum - bestConsensus.mismatchSum))/10.0);
        if ( improvement >= LOD_THRESHOLD ) {
            
            bestConsensus.cigar = AlignmentUtils.leftAlignIndel(bestConsensus.cigar, reference, bestConsensus.str, bestConsensus.positionOnReference, bestConsensus.positionOnReference);
            
            // start cleaning the appropriate reads
            for ( pair<int, int> indexPair : bestConsensus.readIndexes ) {
                AlignedRead aRead = altReads.get(indexPair.first);
                if ( !updateRead(bestConsensus.cigar, bestConsensus.positionOnReference, indexPair.second, aRead, leftmostIndex) )
                    return;
            }
            if ( consensusModel != ConsensusDeterminationModel.KNOWNS_ONLY && !alternateReducesEntropy(altReads, reference, leftmostIndex) ) {
                if ( statsOutput != NULL ) {
                    try {
                        statsOutput.write(currentInterval.toString());
                        statsOutput.write("\tFAIL (bad indel)\t"); // if improvement > LOD_THRESHOLD *BUT* entropy is not reduced (SNPs still exist)
                        statsOutput.write(Double.toString(improvement));
                        statsOutput.write("\n");
                        statsOutput.flush();
                    } catch (Exception e) {
                        throw new UserException.CouldNotCreateOutputFile("statsOutput", "Failed to write stats output file", e);
                    }
                }
            } else {
                //logger.debug("CLEAN: " + bestConsensus.cigar + " " + bestConsensus.str.toString() + " " + bestConsensus.cigar.numCigarElements() );
                if ( indelOutput != NULL && bestConsensus.cigar.numCigarElements() > 1 ) {
                    // NOTE: indels are printed out in the format specified for the low-coverage pilot1
                    //  indel calls (tab-delimited): chr position size type sequence
                    StringBuilder str = new StringBuilder();
                    str.append(reads.get(0).getReferenceName());
                    int position = bestConsensus.positionOnReference + bestConsensus.cigar.getCigarElement(0).getLength();
                    str.append("\t" + (leftmostIndex + position - 1));
                    CigarElement ce = bestConsensus.cigar.getCigarElement(1);
                    str.append("\t" + ce.getLength() + "\t" + ce.getOperator() + "\t");
                    int length = ce.getLength();
                    if ( ce.getOperator() == CigarOperator.D ) {
                        for ( int i = 0; i < length; i++)
                            str.append((char)reference[position+i]);
                    } else {
                        for ( int i = 0; i < length; i++)
                            str.append((char)bestConsensus.str[position+i]);
                    }
                    str.append("\t" + (((double)(totalRawMismatchSum - bestConsensus.mismatchSum))/10.0) + "\n");
                    try {
                        indelOutput.write(str.toString());
                        indelOutput.flush();
                    } catch (Exception e) {
                        throw new UserException.CouldNotCreateOutputFile("indelOutput", "Failed to write indel output file", e);
                    }
                }
                if ( statsOutput != NULL ) {
                    try {
                        statsOutput.write(currentInterval.toString());
                        statsOutput.write("\tCLEAN"); // if improvement > LOD_THRESHOLD *AND* entropy is reduced
                        if ( bestConsensus.cigar.numCigarElements() > 1 )
                            statsOutput.write(" (found indel)");
                        statsOutput.write("\t");
                        statsOutput.write(Double.toString(improvement));
                        statsOutput.write("\n");
                        statsOutput.flush();
                    } catch (Exception e) {
                        throw new UserException.CouldNotCreateOutputFile("statsOutput", "Failed to write stats output file", e);
                    }
                }
                
                // finish cleaning the appropriate reads
                for ( pair<int, int> indexPair : bestConsensus.readIndexes ) {
                    const AlignedRead aRead = altReads.get(indexPair.first);
                    if ( aRead.constizeUpdate() ) {
                        // We need to update the mapping quality score of the cleaned reads;
                        // however we don't have enough info to use the proper MAQ scoring system.
                        // For now, we will just arbitrarily add 10 to the mapping quality. [EB, 6/7/2010].
                        // TODO -- we need a better solution here
                        GATKSAMRecord read = aRead.getRead();
                        if ( read.getMappingQuality() != 255 ) // 255 == Unknown, so don't modify it
                            read.setMappingQuality(Math.min(aRead.getRead().getMappingQuality() + 10, 254));
                        
                        // before we fix the attribute tags we first need to make sure we have enough of the reference sequence
                        int neededBasesToLeft = leftmostIndex - read.getAlignmentStart();
                        int neededBasesToRight = read.getAlignmentEnd() - leftmostIndex - reference.length + 1;
                        int neededBases = Math.max(neededBasesToLeft, neededBasesToRight);
                        if ( neededBases > 0 ) {
                            int padLeft = Math.max(leftmostIndex-neededBases, 1);
                            int padRight = Math.min(leftmostIndex+reference.length+neededBases, referenceReader.getSequenceDictionary().getSequence(currentInterval.getContig()).getSequenceLength());
                            reference = referenceReader.getSubsequenceAt(currentInterval.getContig(), padLeft, padRight).getBases();
                            leftmostIndex = padLeft;
                        }
                        
                        // now, fix the attribute tags
                        // TODO -- get rid of this try block when Picard does the right thing for reads aligned off the end of the reference
                        try {
                            if ( read.getAttribute(SAMTag.NM.name()) != NULL )
                                read.setAttribute(SAMTag.NM.name(), SequenceUtil.calculateSamNmTag(read, reference, leftmostIndex-1));
                            if ( read.getAttribute(SAMTag.UQ.name()) != NULL )
                                read.setAttribute(SAMTag.UQ.name(), SequenceUtil.sumQualitiesOfMismatches(read, reference, leftmostIndex-1));
                        } catch (Exception e) {
                            // ignore it
                        }
                        // TODO -- this is only temporary until Tim adds code to recalculate this value
                        if ( read.getAttribute(SAMTag.MD.name()) != NULL )
                            read.setAttribute(SAMTag.MD.name(), NULL);
                        
                        // mark that it was actually cleaned
                        readsActuallyCleaned.add(read);
                    }
                }
            }
            
            // END IF ( improvement >= LOD_THRESHOLD )
            
        } else if ( statsOutput != NULL ) {
            try {
                statsOutput.write(String.format("%s\tFAIL\t%.1f%n",
                                                currentInterval.toString(), improvement));
                statsOutput.flush();
            } catch (Exception e) {
                throw new UserException.CouldNotCreateOutputFile("statsOutput", "Failed to write stats output file", e);
            }
        }
    }
    
    void generateAlternateConsensesFromKnownIndels(const set<Consensus> altConsensesToPopulate, const int leftmostIndex, const string reference) {
        for ( VariantContext knownIndel : knownIndelsToTry ) {
            if ( knownIndel == NULL || !knownIndel.isIndel() || knownIndel.isComplexIndel() )
                continue;
            byte[] indelStr = knownIndel.isSimpleInsertion() ? knownIndel.getAlternateAllele(0).getBases() : Utils.dupBytes((byte)'-', knownIndel.getReference().length());
            int start = knownIndel.getStart() - leftmostIndex + 1;
            Consensus c = createAlternateConsensus(start, reference, indelStr, knownIndel);
            if ( c != NULL )
                altConsensesToPopulate.add(c);
        }
    }
    
    long determineReadsThatNeedCleaning(const vector<BamAlignment *> reads,
                                                const vector<BamAlignment *> refReadsToPopulate,
                                                const vector<AlignedRead> & altReadsToPopulate,
                                                const list<AlignedRead> & altAlignmentsToTest,
                                                const set<Consensus> & altConsenses,
                                                const int leftmostIndex,
                                                const string & reference) {
        
        long totalRawMismatchSum = 0L;
        
        for ( const BamAlignment * read : reads ) {
            
            // we can not deal with screwy records
            if ( read.getCigar().numCigarElements() == 0 ) {
                refReadsToPopulate.add(read);
                continue;
            }
            
            const AlignedRead aRead = new AlignedRead(read);
            
            // first, move existing indels (for 1 indel reads only) to leftmost position within identical sequence
            int numBlocks = AlignmentUtils.getNumAlignmentBlocks(read);
            if ( numBlocks == 2 ) {
                Cigar newCigar = AlignmentUtils.leftAlignIndel(unclipCigar(read.getCigar()), reference, read.getReadBases(), read.getAlignmentStart()-leftmostIndex, 0);
                aRead.setCigar(newCigar, false);
            }
            
            const int startOnRef = read.getAlignmentStart()-leftmostIndex;
            const int rawMismatchScore = mismatchQualitySumIgnoreCigar(aRead, reference, startOnRef, int.MAX_VALUE);
            
            // if this doesn't match perfectly to the reference, let's try to clean it
            if ( rawMismatchScore > 0 ) {
                altReadsToPopulate.add(aRead);
                //logger.debug("Adding " + read.getReadName() + " with raw mismatch score " + rawMismatchScore + " to non-ref reads");
                
                if ( !read.getDuplicateReadFlag() )
                    totalRawMismatchSum += rawMismatchScore;
                aRead.setMismatchScoreToReference(rawMismatchScore);
                aRead.setAlignerMismatchScore(AlignmentUtils.mismatchingQualities(aRead.getRead(), reference, startOnRef));
                
                // if it has an indel, let's see if that's the best consensus
                if ( consensusModel != ConsensusDeterminationModel.KNOWNS_ONLY && numBlocks == 2 )  {
                    Consensus c = createAlternateConsensus(startOnRef, aRead.getCigar(), reference, aRead.getReadBases());
                    if ( c != NULL )
                        altConsenses.add(c);
                } else {
                    altAlignmentsToTest.add(aRead);
                }
            }
            // otherwise, we can emit it as is
            else {
                //logger.debug("Adding " + read.getReadName() + " with raw mismatch score " + rawMismatchScore + " to ref reads");
                refReadsToPopulate.add(read);
            }
        }
        
        return totalRawMismatchSum;
    }
    
    void generateAlternateConsensesFromReads(const list<AlignedRead> & altAlignmentsToTest,
                                                     const set<Consensus> & altConsensesToPopulate,
                                                     const string & reference,
                                                     const int leftmostIndex) {
        
        // if we are under the limit, use all reads to generate alternate consenses
        if ( altAlignmentsToTest.size() <= MAX_READS_FOR_CONSENSUSES ) {
            for ( AlignedRead aRead : altAlignmentsToTest ) {
                if ( CHECKEARLY ) createAndAddAlternateConsensus1(aRead, altConsensesToPopulate, reference,leftmostIndex);
                else createAndAddAlternateConsensus(aRead.getReadBases(), altConsensesToPopulate, reference);
            }
        }
        // otherwise, choose reads for alternate consenses randomly
        else {
            int readsSeen = 0;
            while ( readsSeen++ < MAX_READS_FOR_CONSENSUSES && altConsensesToPopulate.size() <= MAX_CONSENSUSES) {
                int index = GenomeAnalysisEngine.getRandomGenerator().nextInt(altAlignmentsToTest.size());
                AlignedRead aRead = altAlignmentsToTest.remove(index);
                if ( CHECKEARLY ) createAndAddAlternateConsensus1(aRead, altConsensesToPopulate, reference,leftmostIndex);
                else createAndAddAlternateConsensus(aRead.getReadBases(), altConsensesToPopulate, reference);
            }
        }
    }
    
    void createAndAddAlternateConsensus(const string & read, const set<Consensus> altConsensesToPopulate, const string & reference) {
        
        // do a pairwise alignment against the reference
        SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, read, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND);
        Consensus c = createAlternateConsensus(swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar(), reference, read);
        if ( c != NULL )
            altConsensesToPopulate.add(c);
    }
    
    void createAndAddAlternateConsensus1(AlignedRead read, const set<Consensus> altConsensesToPopulate,
                                                 const string & reference, const int leftmostIndex) {
        
        for ( Consensus known : altConsensesToPopulate ) {
            pair<int, int> altAlignment = findBestOffset(known.str, read, leftmostIndex);
            // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
            int myScore = altAlignment.second;
            if ( myScore == 0 ) {exactMatchesFound++; return; }// read matches perfectly to a known alt consensus - no need to run SW, we already know the answer
        }
        // do a pairwise alignment against the reference
        SWalignmentRuns++;
        SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, read.getReadBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND);
        Consensus c = createAlternateConsensus(swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar(), reference, read.getReadBases());
        if ( c != NULL ) {
            altConsensesToPopulate.add(c);
            SWalignmentSuccess++;
        }
    }
    
    // create a Consensus from cigar/read strings which originate somewhere on the reference
    Consensus createAlternateConsensus(const int indexOnRef, const vector<CigarOp> & c, const string reference, const string readStr) {
        if ( indexOnRef < 0 )
            return NULL;
        
        // if there are no indels, we do not need this consensus, can abort early:
        if ( c.numCigarElements() == 1 && c.getCigarElement(0).getOperator() == CigarOperator.M ) return NULL;
        
        // create the new consensus
        vector<CigarElement> elements = new vector<CigarElement>(c.numCigarElements()-1);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < indexOnRef; i++)
            sb.append((char)reference[i]);
        
        int indelCount = 0;
        int altIdx = 0;
        int refIdx = indexOnRef;
        boolean ok_flag = true;
        for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
            CigarElement ce = c.getCigarElement(i);
            int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
                case D:
                    refIdx += elementLength;
                    indelCount++;
                    elements.add(ce);
                    break;
                case M:
                    altIdx += elementLength;
                case N:
                    if ( reference.length < refIdx + elementLength )
                        ok_flag = false;
                    else  {
                        for (int j = 0; j < elementLength; j++)
                            sb.append((char)reference[refIdx+j]);
                    }
                    refIdx += elementLength;
                    elements.add(new CigarElement(elementLength, CigarOperator.M));
                    break;
                case I:
                    for (int j = 0; j < elementLength; j++) {
                        if ( ! BaseUtils.isRegularBase(readStr[altIdx+j]) ) {
                            // Insertions with N's in them cause real problems sometimes; it's better to drop them altogether
                            ok_flag=false;
                            break;
                        }
                        sb.append((char)readStr[altIdx + j]);
                    }
                    altIdx += elementLength;
                    indelCount++;
                    elements.add(ce);
                    break;
                case S:
                default:
                    break;
            }
        }
        // make sure that there is at most only a single indel and it aligns appropriately!
        if ( !ok_flag || indelCount != 1 || reference.length < refIdx )
            return NULL;
        
        for (int i = refIdx; i < reference.length; i++)
            sb.append((char)reference[i]);
        byte[] altConsensus =  StringUtil.stringToBytes(sb.toString()); // alternative consensus sequence we just built from the current read
        
        return new Consensus(altConsensus, new Cigar(elements), indexOnRef);
    }
    
    // create a Consensus from just the indel string that falls on the reference
    Consensus createAlternateConsensus(const int indexOnRef, const string & reference, const string & indelStr, const VariantContext indel) {
        if ( indexOnRef < 0 || indexOnRef >= reference.length )
            return NULL;
        
        // create the new consensus
        StringBuilder sb = new StringBuilder();
        Cigar cigar = new Cigar();
        int refIdx;
        
        for (refIdx = 0; refIdx < indexOnRef; refIdx++)
            sb.append((char)reference[refIdx]);
        if ( indexOnRef > 0 )
            cigar.add(new CigarElement(indexOnRef, CigarOperator.M));
        
        if ( indel.isSimpleDeletion() ) {
            refIdx += indelStr.length;
            cigar.add(new CigarElement(indelStr.length, CigarOperator.D));
        }
        else if ( indel.isSimpleInsertion() ) {
            for ( byte b : indelStr )
                sb.append((char)b);
            cigar.add(new CigarElement(indelStr.length, CigarOperator.I));
        } else {
            throw new IllegalStateException("Creating an alternate consensus from a complex indel is not allows");
        }
        
        if ( reference.length - refIdx > 0 )
            cigar.add(new CigarElement(reference.length - refIdx, CigarOperator.M));        
        for (; refIdx < reference.length; refIdx++)
            sb.append((char)reference[refIdx]);
        byte[] altConsensus =  StringUtil.stringToBytes(sb.toString()); // alternative consensus sequence we just built from the current read
        
        return new Consensus(altConsensus, cigar, 0);
    }
    
    pair<int, int> findBestOffset(const string & ref, const AlignedRead read, const int leftmostIndex) {
        
        // optimization: try the most likely alignment first (to get a low score to beat)
        int originalAlignment = read.getOriginalAlignmentStart() - leftmostIndex;
        int bestScore = mismatchQualitySumIgnoreCigar(read, ref, originalAlignment, INT_MAX);
        int bestIndex = originalAlignment;
        
        // optimization: we can't get better than 0, so we can quit now
        if ( bestScore == 0 )
            return new pair<int, int>(bestIndex, 0);
        
        // optimization: the correct alignment shouldn't be too far from the original one (or else the read wouldn't have aligned in the first place)
        for ( int i = 0; i < originalAlignment; i++ ) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if ( bestScore == 0 )
                return new pair<int, int>(bestIndex, 0);
        }
        
        const int maxPossibleStart = ref.length - read.getReadLength();
        for ( int i = originalAlignment + 1; i <= maxPossibleStart; i++ ) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if ( bestScore == 0 )
                return new pair<int, int>(bestIndex, 0);
        }
        
        return new pair<int, int>(bestIndex, bestScore);
    }

    bool updateRead(const vector<CigarOp> & altCigar, const int altPosOnRef, const int myPosOnAlt, const AlignedRead & aRead, const int leftmostIndex) {
        vector<CigarOp> readCigar;
        
        // special case: there is no indel
        if ( altCigar.getCigarElements().size() == 1 ) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
            readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            aRead.setCigar(readCigar);
            return true;
        }
        
        CigarElement altCE1 = altCigar.getCigarElement(0);
        CigarElement altCE2 = altCigar.getCigarElement(1);
        
        int leadingMatchingBlockLength = 0; // length of the leading M element or 0 if the leading element is I
        
        CigarElement indelCE;
        if ( altCE1.getOperator() == CigarOperator.I  ) {
            indelCE=altCE1;
            if ( altCE2.getOperator() != CigarOperator.M  ) {
                logger.warn("When the first element of the alt consensus is I, the second one must be M. Actual: " + altCigar.toString() + ".  Skipping this site...");
                return false;
            }
        }
        else {
            if ( altCE1.getOperator() != CigarOperator.M  ) {
                logger.warn("First element of the alt consensus cigar must be M or I. Actual: " + altCigar.toString() + ".  Skipping this site...");
                return false;
            }
            if ( altCE2.getOperator() == CigarOperator.I  || altCE2.getOperator() == CigarOperator.D ) {
                indelCE=altCE2;
            } else {
                logger.warn("When first element of the alt consensus is M, the second one must be I or D. Actual: " + altCigar.toString() + ".  Skipping this site...");
                return false;
            }
            leadingMatchingBlockLength = altCE1.getLength();
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
                aRead.setCigar(NULL); // reset to original alignment
                return true;
            }
            readCigar.add(new CigarElement(endOfFirstBlock - myPosOnAlt, CigarOperator.M));
        }
        
        // forward along the indel
        //int indelOffsetOnRef = 0, indelOffsetOnRead = 0;
        if ( indelCE.getOperator() == CigarOperator.I ) {
            // for reads that end in an insertion
            if ( myPosOnAlt + aRead.getReadLength() < endOfFirstBlock + indelCE.getLength() ) {
                int partialInsertionLength = myPosOnAlt + aRead.getReadLength() - endOfFirstBlock;
                // if we also started inside the insertion, then we need to modify the length
                if ( !sawAlignmentStart )
                    partialInsertionLength = aRead.getReadLength();
                readCigar.add(new CigarElement(partialInsertionLength, CigarOperator.I));
                aRead.setCigar(readCigar);
                return true;
            }
            
            // for reads that start in an insertion
            if ( !sawAlignmentStart && myPosOnAlt < endOfFirstBlock + indelCE.getLength() ) {
                aRead.setAlignmentStart(leftmostIndex + endOfFirstBlock);
                readCigar.add(new CigarElement(indelCE.getLength() - (myPosOnAlt - endOfFirstBlock), CigarOperator.I));
                //indelOffsetOnRead = myPosOnAlt - endOfFirstBlock;
                sawAlignmentStart = true;
            } else if ( sawAlignmentStart ) {
                readCigar.add(indelCE);
                //indelOffsetOnRead = indelCE.getLength();
            }
        } else if ( indelCE.getOperator() == CigarOperator.D ) {
            if ( sawAlignmentStart )
                readCigar.add(indelCE);
            //indelOffsetOnRef = indelCE.getLength();
        }
        
        // for reads that start after the indel
        if ( !sawAlignmentStart ) {
            //aRead.setAlignmentStart(leftmostIndex + myPosOnAlt + indelOffsetOnRef - indelOffsetOnRead);
            //readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            //aRead.setCigar(readCigar);
            aRead.setCigar(NULL); // reset to original alignment
            return true;
        }
        
        int readRemaining = aRead.getReadBases().length;
        for ( CigarElement ce : readCigar.getCigarElements() ) {
            if ( ce.getOperator() != CigarOperator.D )
                readRemaining -= ce.getLength();
        }
        if ( readRemaining > 0 )
            readCigar.add(new CigarElement(readRemaining, CigarOperator.M));
        aRead.setCigar(readCigar);
        
        return true;
    }
    
    bool alternateReducesEntropy(const list<AlignedRead> & reads, const string & reference, const int leftmostIndex) {
        const int[] originalMismatchBases = new int[reference.length];
        const int[] cleanedMismatchBases = new int[reference.length];
        const int[] totalOriginalBases = new int[reference.length];
        const int[] totalCleanedBases = new int[reference.length];
        
        // set to 1 to prevent dividing by zero
        for ( int i=0; i < reference.length; i++ )
            originalMismatchBases[i] = totalOriginalBases[i] = cleanedMismatchBases[i] = totalCleanedBases[i] = 0;
        
        for (int i=0; i < reads.size(); i++) {
            const AlignedRead read = reads.get(i);
            if ( read.getRead().getAlignmentBlocks().size() > 1 )
                continue;
            
            int refIdx = read.getOriginalAlignmentStart() - leftmostIndex;
            const byte[] readStr = read.getReadBases();
            const byte[] quals = read.getBaseQualities();
            
            for (int j=0; j < readStr.length; j++, refIdx++ ) {
                if ( refIdx < 0 || refIdx >= reference.length ) {
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
            Cigar c = read.getCigar();
            for (int j = 0 ; j < c.numCigarElements() ; j++) {
                CigarElement ce = c.getCigarElement(j);
                int elementLength = ce.getLength();
                switch ( ce.getOperator() ) {
                    case M:
                        for (int k = 0 ; k < elementLength ; k++, refIdx++, altIdx++ ) {
                            if ( refIdx >= reference.length )
                                break;
                            totalCleanedBases[refIdx] += quals[altIdx];
                            if ( readStr[altIdx] != reference[refIdx] )
                                cleanedMismatchBases[refIdx] += quals[altIdx];
                        }
                        break;
                    case I:
                        altIdx += elementLength;
                        break;
                    case D:
                        refIdx += elementLength;
                        break;
                    case S:
                    default:
                        break;
                }
            }
        }
        
        int originalMismatchColumns = 0, cleanedMismatchColumns = 0;
        StringBuilder sb = new StringBuilder();
        for ( int i=0; i < reference.length; i++ ) {
            if ( cleanedMismatchBases[i] == originalMismatchBases[i] )
                continue;
            boolean didMismatch = false, stillMismatches = false;
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
                    sb.append(reads.get(0).getRead().getReferenceName() + ":");
                    sb.append((leftmostIndex + i));
                    if ( stillMismatches )
                        sb.append(" SAME_SNP\n");
                    else
                        sb.append(" NOT_SNP\n");
                }
            }
        }
        
        //logger.debug("Original mismatch columns = " + originalMismatchColumns + "; cleaned mismatch columns = " + cleanedMismatchColumns);
        
        const boolean reduces = (originalMismatchColumns == 0 || cleanedMismatchColumns < originalMismatchColumns);
        if ( reduces && snpsOutput != NULL ) {
            try {
                snpsOutput.write(sb.toString());
                snpsOutput.flush();
            } catch (Exception e) {
                throw new UserException.CouldNotCreateOutputFile("snpsOutput", "Failed to write SNPs output file", e);
            }
        }
        return reduces;
    }
    
    static vector<CigarOp> unclipCigar(vector<CigarOp> cigar) {
        vector<CigarOp> elements;
        elements.reserve(cigar.size());

        for ( CigarElement ce : cigar.getCigarElements() ) {
            if ( !isClipOperator(ce.getOperator()) )
                elements.add(ce);
        }
        return elements;
    }
    
    static bool isClipOperator(CigarOp op) {
        return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
    }
    
    static vector<CigarOp> reclipCigar(vector<CigarOp> & cigar, BamAlignment * read) {
        vector<CigarOp> elements;
        
        int i = 0;
        int n = read.getCigar().numCigarElements();
        while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
            elements.add(read.getCigar().getCigarElement(i++));
        
        elements.addAll(cigar.getCigarElements());
        
        i++;
        while ( i < n && !isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
            i++;
        
        while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
            elements.add(read.getCigar().getCigarElement(i++));
        
        return elements;
    
}


int LocalRealignment::runInternal()
{
    
}