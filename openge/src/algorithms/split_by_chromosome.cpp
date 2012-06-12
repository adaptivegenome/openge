/*********************************************************************
 *
 * split_by_chromosome.cpp:  Split a stream of reads into multiple streams
 *                           by chromosome.
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

#include "split_by_chromosome.h"

#include <algorithm>
#include <numeric>

using namespace BamTools;
using namespace std;

SplitByChromosome::SplitByChromosome()
: number_of_splits(0)
{
}

int SplitByChromosome::runInternal()
{    
    ogeNameThread("am_split_chromo");

    BamAlignment * read;
    number_of_splits = sinks.size();
    while(true) {
        
        read = getInputAlignment();
        
        if(!read) break;

        //this is essentially a modified implementation of AlgorithmModule::putOutputAlignment()
        write_count++;
        
        int chain = read->RefID % number_of_splits;

        if(read->RefID < 0)
            chain = 0;

        sinks[chain]->putInputAlignment(read);
    }
    
    if(verbose)
        for(int i = 0; i < sinks.size(); i++)
            cerr << "Chain " << i << " wrote " << sinks[i]->getReadCount() << endl;

    return 0;
}