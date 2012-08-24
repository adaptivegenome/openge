#ifndef OGE_ALGO_HISTORY_H
#define OGE_ALGO_HISTORY_H
/*********************************************************************
 *
 * print_history.h: Print the history of commands from a file's header.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 2 August 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include "algorithm_module.h"

class PrintHistory : public AlgorithmModule
{
protected:
    virtual int runInternal();
};

#endif