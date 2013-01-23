#ifndef OGE_REPEATSEQ_H
#define OGE_REPEATSEQ_H
/*********************************************************************
 *
 * repeatseq.cpp:
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
#include "../util/fasta_reader.h"
#include "../util/oge_read.h"
#include "../util/thread_pool.h"
#include <string>
#include <sstream>
#include <vector>
#include <fstream>

class Repeatseq : public AlgorithmModule {
public:
    typedef enum {
        ERROR_NORMAL, ERROR_SOMATIC_A, ERROR_SOMATIC_B
    } error_model_t;

protected:
    int MAX_READ_SIZE;
	int LR_CHARS_TO_PRINT;
	int mode;
	bool multi;
	bool properlyPaired;
	bool makeRepeatseqFile;
	bool makeCallsFile;
    bool makeVcfFile;
	int readLengthMin;
	int readLengthMax;
	int consLeftFlank;
	int consRightFlank;
	int MapQuality;
    std::string paramString;
    std::string fasta_filename;
    std::string intervals_filename;
    std::string output_filename_base;
    std::string somatic_input;

    std::ofstream vcfFile, oFile, callsFile;

    FastaReader fasta_reader;
    
    error_model_t error_model;
public:
	Repeatseq()
    : MAX_READ_SIZE(200)
    , LR_CHARS_TO_PRINT(8)
    , mode(2)
    , multi(false)
    , properlyPaired(false)
    , makeRepeatseqFile(false)
    , makeCallsFile(false)
    , makeVcfFile(true)
    , readLengthMin(0)
    , readLengthMax(0)
    , consLeftFlank(3)
    , consRightFlank(3)
    , MapQuality(0)
    , error_model(ERROR_NORMAL)
    , last_deleted_job(0)
    { }
    
    void setLength(int len_min, int len_max) { readLengthMin = len_min; readLengthMax = len_min; }
    void setHaploid(bool haploid) { if(haploid) mode = 1; else mode = 2; }
    void setProperlyPaired(bool paired) { properlyPaired = paired; }
    void setMulti(bool multi) { this->multi = multi; }
    void setLeftFlank(int flank) { consLeftFlank = flank; }
    void setRightFlank(int flank) { consRightFlank = flank; }
    void setLRCharsToPrint(int chars) { LR_CHARS_TO_PRINT = chars; }
    void setParamString(const std::string param_string) { paramString = param_string; }
    void setMinimumQuality(const int quality) { MapQuality = quality; }
    void setMakeRepeatseqFile(bool make) { makeRepeatseqFile = make; }
    void setMakeCallsFile(bool make) { makeCallsFile = make; }
    void setMakeVcfFile(bool make) { makeVcfFile = make; }
    void setFastaFilename(const std::string & filename) { fasta_filename = filename; }
    void setIntervalsFilename(const std::string & filename ) { intervals_filename = filename; }
    void setOutputFilename(const std::string & filename) { output_filename_base = filename; }
    void setSomaticInput(const std::string & s) { somatic_input = s; }
    void setErrorModel(error_model_t e) { error_model = e; }
    
    ///////////////////
protected:
    //class for parsing region argument:
    class Region {
    public:
        std::string startSeq;
        int startPos;
        int stopPos;
        
        Region(std::string& region);
        std::string toString() const;
        int length(void) const;
    };
    struct RegionStringComparator {
        const BamSequenceRecords & d;
        RegionStringComparator(const BamSequenceRecords & d);
        bool operator() (std::string a, std::string b) const;
    };
    struct Sequences {
        bool insertions;
        std::string preSeq;
        std::string alignedSeq;
        std::string postSeq;
        
        Sequences();
        Sequences(std::string, std::string, std::string, bool);
    };
    
    struct STRING_GT {
        std::string print;
        Sequences reads;
        int GT;
        bool paired;
        bool reverse;
        int MapQ;
        int minFlank;
        double avgBQ;
        
        STRING_GT(std::string, Sequences, int, bool, int, int, bool, double);
        STRING_GT();
        bool operator<(const STRING_GT &other) const; //for sorting
    };

    //to be used for printing GT: information to header line:
    struct GT {
        int readlength;
        int occurrences;
        int reverse;
        int avgMinFlank;
        double avgBQ;
        
        GT(int rl, int oc, int rev, int minF, double avg);
        
        //ordering should be in reverse order- largest to smallest for display.
        inline bool operator<(const GT & b) const { return (occurrences > b.occurrences); }
    };

    virtual int runInternal();
    void flushWrites();
    std::string getReference(const Region & target) const;
    std::vector<STRING_GT> makeToPrint(const std::vector<OGERead *> reads, const Region & target, const std::string & leftReference, const std::string & centerReference, const std::string & rightReference, int & numStars, int & depth) const;
    std::vector<GT> makeGenotypeInformation(std::vector<STRING_GT> & toPrint) const;
    inline std::vector<int> printGenoPerc(std::vector<GT> vectorGT, int ref_length, int unit_size, double &confidence, int mode, std::stringstream & ofile_out) const ;
    inline std::vector<int> somaticConfidence(std::vector<GT> & vectorGT, const std::vector<GT> & vectorGT_reference, const Region & target, error_model_t model, double &confidence, std::stringstream & ofile_out) const ;
    void print_output(const std::string &, std::stringstream &vcf_buffer,  std::stringstream &o_buffer, std::stringstream &calls_buffer, const std::vector<OGERead *> & reads, const std::vector<OGERead *> & reads_somatic) const;
    
    class RepeatseqJob : public ThreadJob {
    public:
        RepeatseqJob(Repeatseq * rs, int job_id, const std::string & region)
        : repeatseq(rs)
        , job_id(job_id)
        , region(region)
        {}

        Repeatseq * repeatseq;
        int job_id;
        std::vector<OGERead *> reads, somatic_reads;
        const std::string & region;
        std::stringstream vcfFile, oFile, callsFile;
        
        virtual void runJob();
    };
    std::vector<RepeatseqJob *> jobs;
    Spinlock write_flush_lock;
    int last_deleted_job;
};

#endif