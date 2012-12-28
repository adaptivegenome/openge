#ifndef OGE_ALGO_SPLIT_CHROMO_H
#define OGE_ALGO_SPLIT_CHROMO_H

/*********************************************************************
 *
 * split_by_chromosome.h:  Split a stream of reads into multiple streams
 *                         by chromosome.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 7 June 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include "algorithm_module.h"

#include <string>

class SplitByChromosome : public AlgorithmModule
{
public:
    SplitByChromosome();
    void setSplitCount(int count) { number_of_splits = count; }
    int getSplitCount() { return number_of_splits; }
protected:
    virtual int runInternal();
    int number_of_splits;
};

#endif