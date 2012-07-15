#ifndef OGE_ALGO_COVERAGE_H
#define OGE_ALGO_COVERAGE_H

/*********************************************************************
 *
 * measure_coverage.h:  Measure the coverage of a file.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 22 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include "algorithm_module.h"
#include "api/BamAlignment.h"

#include <string>
#include <vector>
#include <map>

class MeasureCoverage : public AlgorithmModule
{
public:
    MeasureCoverage();
    const std::map<std::string, std::vector<int> > & getCoverageMap() const;
    void setVerifyCorrectMapping(bool verify) { verify_mapping = verify; }
    void setOutputFile(const std::string & filename) { out_filename = filename; }
    void setPrintZeroCoverageBases(bool print_zero_cover_bases) { this->print_zero_cover_bases = print_zero_cover_bases; }
    void setBinSize(int bin_size) { this->binsize = bin_size; }
protected:
    virtual int runInternal();
protected:
    std::map<std::string, std::vector<int> > coverage_map;
    std::map<std::string, std::vector<int> > correctness_map;

    bool verify_mapping;
    bool print_zero_cover_bases;
    std::string out_filename;
    int binsize;
};

#endif