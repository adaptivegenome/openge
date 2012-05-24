#ifndef OGE_ALGO_FILTER_H
#define OGE_ALGO_FILTER_H

#include "algorithm_module.h"
#include "api/BamAlignment.h"

#include <string>

class Filter : public AlgorithmModule
{
public:
    Filter();
    void setRegion(std::string region) { region_string = region; has_region = true; }
    std::string getRegion() { return region_string; }
    void setCountLimit(int ct) { count_limit = ct;}
    size_t getCountLimit() { return count_limit;}
protected:
    virtual int runInternal();
protected:
    bool ParseRegionString(const std::string& regionString, BamTools::BamRegion& region);
    std::string region_string;
    bool has_region;
    size_t count_limit;
};

#endif