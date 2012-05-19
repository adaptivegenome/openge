#ifndef OGE_ALGO_LOCAL_REALIGNMENT_H
#define OGE_ALGO_LOCAL_REALIGNMENT_H

#include "algorithm_module.h"

#include <string>

#include "../util/fasta_reader.h"

#include "../util/gatk/AlignmentUtils.h"
#include "../util/gatk/GATKFeature.h"
#include "../util/gatk/GenomeLoc.h"
#include "../util/gatk/GenomeLocParser.h"
#include "../util/gatk/ConstrainedMateFixingManager.h"
#include "../util/gatk/BaseUtils.h"
#include "../util/gatk/GATKFeature.h"
#include "../util/gatk/VariantContext.h"
#include "../util/gatk/ReadMetaDataTracker.h"
#include "../util/gatk/SequenceUtil.h"


class LocalRealignment : public AlgorithmModule
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
    //vector<RodBinding<VariantContext> > known;
    
    std::string reference_filename;
    std::string output_filename;
    std::string intervals_filename;
    
protected:
    /**
     * The interval list output from the RealignerTargetCreator tool using the same bam(s), reference, and known indel file(s).
     */
    //@Input(fullName="targetIntervals", shortName="targetIntervals", doc="intervals file output from RealignerTargetCreator", required=true)
    std::vector<GenomeLoc *> intervalsFile;
    
    /**
     * This term is equivalent to "significance" - i.e. is the improvement significant enough to merit realignment? Note that this number
     * should be adjusted based on your particular data set. For low coverage and/or when looking for indels with low allele frequency,
     * this number should be smaller.
     */
    //@Argument(fullName="LODThresholdForCleaning", shortName="LOD", doc="LOD threshold above which the cleaner will clean", required=false)
    double LOD_THRESHOLD;

    ConstrainedMateFixingManager * manager;
    
public:
    /**
     * We recommend that users run with USE_READS when trying to realign high quality longer read data mapped with a gapped aligner;
     * Smith-Waterman is really only necessary when using an ungapped aligner (e.g. MAQ in the case of single-end read data).
     */
    //@Argument(fullName = "consensusDeterminationModel", shortName = "model", doc = "Determines how to compute the possible alternate consenses", required = false)
    ConsensusDeterminationModel_t consensusModel;
    
    
    // ADVANCED OPTIONS FOLLOW
    
    /**
     * For expert users only!  This is similar to the argument in the RealignerTargetCreator walker. The point here is that the realigner
     * will only proceed with the realignment (even above the given threshold) if it minimizes entropy among the reads (and doesn't simply
     * push the mismatch column to another position). This parameter is just a heuristic and should be adjusted based on your particular data set.
     */
protected:
    //@Argument(fullName="entropyThreshold", shortName="entropy", doc="percentage of mismatches at a locus to be considered having high entropy", required=false)
    double MISMATCH_THRESHOLD;
    
    /**
     * For expert users only!  To minimize memory consumption you can lower this number (but then the tool may skip realignment on regions with too much coverage;
     * and if the number is too low, it may generate errors during realignment). Just make sure to give Java enough memory! 4Gb should be enough with the default value.
     */
    //@Advanced
    //@Argument(fullName="maxReadsInMemory", shortName="maxInMemory", doc="max reads allowed to be kept in memory at a time by the SAMFileWriter", required=false)
    int MAX_RECORDS_IN_MEMORY;
    
    /**
     * For expert users only!
     */
    //@Advanced
    //@Argument(fullName="maxIsizeForMovement", shortName="maxIsize", doc="maximum insert size of read pairs that we attempt to realign", required=false)
    int MAX_ISIZE_FOR_MOVEMENT;
    
    /**
     * For expert users only!
     */
    //@Advanced
    //@Argument(fullName="maxPositionalMoveAllowed", shortName="maxPosMove", doc="maximum positional move in basepairs that a read can be adjusted during realignment", required=false)
    int MAX_POS_MOVE_ALLOWED;
    
    /**
     * For expert users only!  If you need to find the optimal solution regardless of running time, use a higher number.
     */
    //@Advanced
    //@Argument(fullName="maxConsensuses", shortName="maxConsensuses", doc="max alternate consensuses to try (necessary to improve performance in deep coverage)", required=false)
    int MAX_CONSENSUSES;
    
    /**
     * For expert users only!  If you need to find the optimal solution regardless of running time, use a higher number.
     */
    //@Advanced
    //@Argument(fullName="maxReadsForConsensuses", shortName="greedy", doc="max reads used for finding the alternate consensuses (necessary to improve performance in deep coverage)", required=false)
    int MAX_READS_FOR_CONSENSUSES;
    
    /**
     * For expert users only!  If this value is exceeded at a given interval, realignment is not attempted and the reads are passed to the output file(s) as-is.
     * If you need to allow more reads (e.g. with very deep coverage) regardless of memory, use a higher number.
     */
    //@Advanced
    //@Argument(fullName="maxReadsForRealignment", shortName="maxReads", doc="max reads allowed at an interval for realignment", required=false)
    int MAX_READS;
    
    //@Advanced
    //@Argument(fullName="noOriginalAlignmentTags", shortName="noTags", required=false, doc="Don't output the original cigar or alignment start tags for each realigned read in the output bam")
    bool NO_ORIGINAL_ALIGNMENT_TAGS;
    
    //@Hidden
    //@Argument(fullName="generate_nWayOut_md5s",doc="Generate md5sums for BAMs")
    bool generateMD5s;
    
    // DEBUGGING OPTIONS FOLLOW
    
    //@Hidden
    //@Argument(fullName="check_early",shortName="check_early",required=false,doc="Do early check of reads against existing consensuses")
    bool CHECKEARLY;
    
    //@Hidden
    //@Argument(fullName="keepPGTags", shortName="keepPG", required=false, doc="Keep older PG tags left in the bam header by previous runs of this tool (by default, all these "+ "historical tags will be replaced by the latest tag generated in the current run).")
    bool KEEP_ALL_PG_RECORDS;
    
    //@Hidden
    //@Output(fullName="indelsFileForDebugging", shortName="indels", required=false, doc="Output file (text) for the indels found; FOR DEBUGGING PURPOSES ONLY")
    std::string OUT_INDELS;
    bool write_out_indels;
    
    // @Hidden
    //@Output(fullName="statisticsFileForDebugging", shortName="stats", doc="print out statistics (what does or doesn't get cleaned); FOR DEBUGGING PURPOSES ONLY", required=false)
    std::string OUT_STATS;
    bool write_out_stats;
    
    //@Hidden
    //@Output(fullName="SNPsFileForDebugging", shortName="snps", doc="print out whether mismatching columns do or don't get cleaned out; FOR DEBUGGING PURPOSES ONLY", required=false)
    std::string OUT_SNPS;
    bool write_out_snps;
    
private:
    static std::vector<BamTools::CigarOp> unclipCigar(const std::vector<BamTools::CigarOp> & cigar);    
    static bool isClipOperator(const BamTools::CigarOp op);
    static std::vector<BamTools::CigarOp> reclipCigar(std::vector<BamTools::CigarOp> & cigar, BamTools::BamAlignment * read);
    
    class AlignedRead {
    private:
        BamTools::BamAlignment * read;
        std::string * readBases;
        std::string * baseQuals;
        std::vector<BamTools::CigarOp> * newCigar;
        int newStart;
        int mismatchScoreToReference;
        long alignerMismatchScore;
    public:
        static int MAX_POS_MOVE_ALLOWED;
        static int NO_ORIGINAL_ALIGNMENT_TAGS;
        AlignedRead(BamTools::BamAlignment * read) 
        : read(read) 
        , newCigar(NULL)
        , newStart(-1)
        , mismatchScoreToReference(0)
        , alignerMismatchScore(0)
        { }
        
        BamTools::BamAlignment * getRead() {
            return read;
        }
        
        int getReadLength() {
            return readBases != NULL ? readBases->size() : read->Length;
        }
        
        std::string getReadBases();
        
        std::string getBaseQualities() ;
        
    private:
        // pull out the bases that aren't clipped out
        void getUnclippedBases();
        
        // pull out the bases that aren't clipped out
        std::vector<BamTools::CigarOp> reclipCigar(std::vector<BamTools::CigarOp> & cigar);
    public:
        std::vector<BamTools::CigarOp> & getCigar();
        // tentatively sets the new Cigar, but it needs to be confirmed later
        void setCigar(std::vector<BamTools::CigarOp> * cigar, bool fixClippedCigar = true);
        
    public:
        void setAlignmentStart(int start);
        int getAlignmentStart();
        int getOriginalAlignmentStart();        
        bool constizeUpdate();

        void setMismatchScoreToReference(int score);
        int getMismatchScoreToReference();        
        void setAlignerMismatchScore(long score);
        long getAlignerMismatchScore();
    };
    
#pragma mark Consensus
    class Consensus {
    public:
        std::string str;
        std::vector<std::pair<int, int> > readIndexes;
        int positionOnReference;
        int mismatchSum;
        std::vector<BamTools::CigarOp> cigar;
        
        Consensus(std::string str, std::vector<BamTools::CigarOp> cigar, int positionOnReference) 
        : str(str)
        , positionOnReference(positionOnReference)
        , mismatchSum(0)
        , cigar(cigar)
        {}

        bool operator==(const Consensus & c) {
            return ( this == &c || str == c.str ) ;
        }
    };
    
#pragma mark ReadBin
    class ReadBin {
    private:
        std::vector<BamTools::BamAlignment *> reads;
        std::string * reference;
        GenomeLoc * loc;
        GenomeLocParser * loc_parser;
        BamTools::SamSequenceDictionary sequences;
    public:
        ReadBin() 
        : reference(NULL)
        , loc(NULL)
        , loc_parser(NULL)
        { }
        
        void initialize(GenomeLocParser * loc_parser, const BamTools::SamHeader & header) {
            this->loc_parser = loc_parser;
            sequences = header.Sequences;
        }
        
        // Return false if we can't process this read bin because the reads are not correctly overlapping.
        // This can happen if e.g. there's a large known indel with no overlapping reads.
        void add(BamTools::BamAlignment * read);
        
        std::vector<BamTools::BamAlignment *> getReads() { return reads; }
        
        std::string * getReference(FastaReader & referenceReader);
        GenomeLoc getLocation() { return *loc; }
        
        int size() { return reads.size(); }
        
        void clear();
    };
    
private:
    // fasta reference reader to supplement the edges of the reference sequence
    FastaReader * referenceReader;
    GenomeLocParser * loc_parser;
    
    // the intervals input by the user
    std::vector<GenomeLoc *>::iterator interval_it;
    
    // the current interval in the list
    GenomeLoc * currentInterval;
    bool sawReadInCurrentInterval;
    
    // the reads and known indels that fall into the current interval
    ReadBin readsToClean;
    std::vector<BamTools::BamAlignment *> readsNotToClean;
    std::vector<VariantContext *> knownIndelsToTry;
    std::set<GATKFeature> indelRodsSeen;
    std::set<BamTools::BamAlignment *> readsActuallyCleaned;
    
    static const int MAX_QUAL;
    
    // fraction of mismatches that need to no longer mismatch for a column to be considered cleaned
    static const double MISMATCH_COLUMN_CLEANED_FRACTION;
    
    static const double SW_MATCH;      // 1.0;
    static const double SW_MISMATCH;  //-1.0/3.0;
    static const double SW_GAP;       //-1.0-1.0/3.0;
    static const double SW_GAP_EXTEND; //-1.0/.0;
    
    // reference base padding size
    // TODO -- make this a command-line argument if the need arises
    static const int REFERENCE_PADDING;
    
    // other output files
    bool outputIndels, output_stats, output_snps;
    std::ofstream indelOutput, statsOutput, snpsOutput;
    
    //###protected Map<SAMReaderID, ConstrainedMateFixingManager> nwayWriters = NULL;
    
    
    // debug info for lazy SW evaluation:
    long exactMatchesFound; // how many reads exactly matched a consensus we already had
    long SWalignmentRuns; // how many times (=for how many reads) we ran SW alignment
    long SWalignmentSuccess; // how many SW alignments were "successful" (i.e. found a workable indel and resulted in non-null consensus)
    
public:
    void initialize();
    void writeRead(BamTools::BamAlignment * read) { putOutputAlignment(read); }

private:
    void setupWriter( BamTools::SamHeader header);
    BamTools::SamProgram createProgramRecord();
    void emit(BamTools::BamAlignment * read);
    void emitReadLists();
    
public:
    int map_func(BamTools::BamAlignment * read, ReadMetaDataTracker metaDataTracker);

private:
    void abortCleanForCurrentInterval();
    bool doNotTryToClean( BamTools::BamAlignment * read);
    void cleanAndCallMap(BamTools::BamAlignment * read, ReadMetaDataTracker metaDataTracker, GenomeLoc * readLoc);
    
public:
    int reduceInit();
    int reduce(int value, int sum);
    void onTraversalDone(int result);

private:
    void populateKnownIndels(ReadMetaDataTracker metaDataTracker) ;
    
    static int mismatchQualitySumIgnoreCigar(AlignedRead aRead, const std::string refSeq, int refIndex, int quitAboveThisValue);
    
    void clean(ReadBin readsToClean) ;
    void generateAlternateConsensesFromKnownIndels(std::set<Consensus *> & altConsensesToPopulate, const int leftmostIndex, const std::string reference);
    long determineReadsThatNeedCleaning( std::vector<BamTools::BamAlignment *> & reads,
                                        std::vector<BamTools::BamAlignment *> & refReadsToPopulate,
                                        std::vector<AlignedRead> & altReadsToPopulate,
                                        std::vector<AlignedRead> & altAlignmentsToTest,
                                        std::set<Consensus *> & altConsenses,
                                        int leftmostIndex,
                                        std::string & reference) ;
    void generateAlternateConsensesFromReads( std::vector<AlignedRead> & altAlignmentsToTest,
                                             std::set<Consensus *> & altConsensesToPopulate,
                                             const std::string & reference,
                                             const int leftmostIndex);
    void createAndAddAlternateConsensus(const std::string & read, std::set<Consensus *> & altConsensesToPopulate, const std::string & reference);
    void createAndAddAlternateConsensus1(AlignedRead & read, std::set<Consensus *> & altConsensesToPopulate,
                                         const std::string & reference, const int leftmostIndex);
    Consensus * createAlternateConsensus(const int indexOnRef, const std::vector<BamTools::CigarOp> & c, const std::string reference, const std::string readStr);
    Consensus * createAlternateConsensus(const int indexOnRef, const std::string & reference, const std::string & indelStr, VariantContext indel);
    std::pair<int, int> findBestOffset(const std::string & ref, AlignedRead read, const int leftmostIndex) ;
    std::string cigarToString(const std::vector<BamTools::CigarOp> & cigar);
    bool updateRead(const std::vector<BamTools::CigarOp> & altCigar, const int altPosOnRef, const int myPosOnAlt, AlignedRead & aRead, const int leftmostIndex);
    bool alternateReducesEntropy(std::vector<AlignedRead> & reads, const std::string & reference, const int leftmostIndex) ;

protected:
public:
    bool verbose;
    LocalRealignment();
    
    std::string getReferenceFilename() { return reference_filename; }
    void setReferenceFilename(const std::string & filename) { this->reference_filename = filename; }
    std::string getIntervalsFilename() { return intervals_filename; }
    void setIntervalsFilename(const std::string & filename) { this->intervals_filename = filename; }

protected:
    int runInternal();
};

#endif
