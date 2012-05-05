#ifndef OGE_ALGO_STATISTICS_H
#define OGE_ALGO_STATISTICS_H

#include "algorithm_module.h"

#include <vector>

class Statistics : public AlgorithmModule
{
public:
    Statistics();
    void showInsertSizeSummary(bool show) { m_showInsertSizeSummary = show; }
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
};

#endif