#ifndef OGE_ALGO_STATISTICS_H
#define OGE_ALGO_STATISTICS_H
/*********************************************************************
 *
 * statistics.cpp: Calculate statistics on a stream of reads.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 24 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************
 *
 * This file is based largely on bamtools_stats.cpp. Original authors
 * Derek Barnett, Erik Garrison, Marth Lab, Department of Biology, 
 * Boston College
 *********************************************************************/

#include "algorithm_module.h"

#include <vector>

class Statistics : public AlgorithmModule
{
public:
    Statistics();
    void showInsertSizeSummary(bool show) { m_showInsertSizeSummary = show; }
    void showReadLengthSummary(bool show) { m_showLengthSummary = show; }
protected:
    virtual int runInternal();
protected:
    unsigned int m_numReads;
    unsigned int m_numPaired;
    unsigned int m_numProperPair;
    unsigned int m_numMapped;
    unsigned int m_numBothMatesMapped;
    unsigned int m_numForwardStrand;
    unsigned int m_numReverseStrand;
    unsigned int m_numFirstMate;
    unsigned int m_numSecondMate;
    unsigned int m_numSingletons;
    unsigned int m_numFailedQC;
    unsigned int m_numDuplicates;
    std::vector<int> m_insertSizes;
    bool m_showInsertSizeSummary;
    bool m_showLengthSummary;
};

#endif