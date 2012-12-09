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

class Repeatseq : public AlgorithmModule {
    int MAX_READ_SIZE;
	int LR_CHARS_TO_PRINT;
	int mode;
	bool multi;
	bool properlyPaired;
	bool makeRepeatseqFile;
	bool makeCallsFile;
	int readLengthMin;
	int readLengthMax;
	int consLeftFlank;
	int consRightFlank;
	int MapQuality;
    std::string paramString;
    std::string fasta_filename;
    std::string intervals_filename;
    std::string output_filename_base;

    std::ofstream vcfFile, oFile, callsFile;

    FastaReader fasta_reader;

public:
	Repeatseq()
    : MAX_READ_SIZE(200)
    , LR_CHARS_TO_PRINT(8)
    , mode(2)
    , multi(false)
    , properlyPaired(false)
    , makeRepeatseqFile(false)
    , makeCallsFile(false)
    , readLengthMin(0)
    , readLengthMax(0)
    , consLeftFlank(3)
    , consRightFlank(3)
    , MapQuality(0)
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
    void setFastaFilename(const std::string & filename) { fasta_filename = filename; }
    void setIntervalsFilename(const std::string & filename ) { intervals_filename = filename; }
    void setOutputFilename(const std::string & filename) { output_filename_base = filename; }
    
    ///////////////////
protected:
    virtual int runInternal();
    void flushWrites();

    void print_output(const std::string &, std::stringstream &vcf_buffer,  std::stringstream &o_buffer, std::stringstream &calls_buffer, const std::vector<OGERead *> & reads) const;
    
    class RepeatseqJob : public ThreadJob {
    public:
        RepeatseqJob(Repeatseq * rs, int job_id, const std::string & region)
        : repeatseq(rs)
        , job_id(job_id)
        , region(region)
        {}

        Repeatseq * repeatseq;
        int job_id;
        std::vector<OGERead *> reads;
        const std::string & region;
        std::stringstream vcfFile, oFile, callsFile;
        
        virtual void runJob();
    };
    std::vector<RepeatseqJob *> jobs;
    Spinlock write_flush_lock;
    int last_deleted_job;
};

#endif