#ifndef OGE_ALGO_SORTED_MERGE_H
#define OGE_ALGO_SORTED_MERGE_H

/*********************************************************************
 *
 * sorted_merge.h:  Merge multiple read streams, keeping them 
 *                  sorted.
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
 *********************************************************************
 * 
 * To merge, a SortedMergeInputProxy is created for each source,
 * and these are polled in turn by the SortedMerge parent class.
 * No work is actually done by the SortedMergeInputProxy class.
 * 
 *********************************************************************/

#include "algorithm_module.h"
#include "api/BamAlignment.h"

#include <string>

class SortedMerge : public AlgorithmModule
{
    class SortedMergeInputProxy : public AlgorithmModule
    {
    public:
        SortedMergeInputProxy(SortedMerge * parent);
        SortedMerge * merge_class;
        void mergeDone();   //< call when the last read comes in from this source.
        using AlgorithmModule::getInputAlignment;
    protected:
        virtual int runInternal();

        pthread_mutex_t merge_complete;
    };
    
public:
    SortedMerge();
    ~SortedMerge();
    void addSource(AlgorithmModule * source);
protected:
    virtual int runInternal();
    std::vector<SortedMergeInputProxy *> input_proxies;
};

#endif