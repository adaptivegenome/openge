#ifndef OGE_ALGO_LOCAL_REALIGNMENT_H
#define OGE_ALGO_LOCAL_REALIGNMENT_H

#include "algorithm_module.h"

class LocalRealignment : public AlgorithmModule
{
public:
    bool verbose;
    LocalRealignment();

protected:
    int runInternal();
};

#endif