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
 * When all reads have been read, all the proxies are notified with
 * a condition variable.
 * 
 *********************************************************************/

#include "algorithm_module.h"

#include <string>

class SortedMerge : public AlgorithmModule
{
    class SortedMergeInputProxy : public AlgorithmModule
    {
    public:
        SortedMergeInputProxy(SortedMerge * parent);
        SortedMerge * merge_class;
        using AlgorithmModule::getInputAlignment;
    protected:
        virtual int runInternal();
    };
    
    class SortedMergeElement{
    public:
        OGERead * read;
        SortedMergeInputProxy * source;
        bool operator<(const SortedMergeElement & t) const;
        SortedMergeElement(OGERead * read, SortedMergeInputProxy * source)
        : read(read)
        , source(source)
        {}
    };
    
public:
    ~SortedMerge();
    void addSource(AlgorithmModule * source);
protected:
    virtual int runInternal();
    std::vector<SortedMergeInputProxy *> input_proxies;
    SynchronizedFlag done;
    condition_variable done_signal;
    mutex done_signal_mutex;
};

#endif