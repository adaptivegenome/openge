// ***************************************************************************
// Largely based on bamtool's bamtools_stats.cpp:
// (c) 2010 Derek Barnett, Erik Garrison, Lee C. Baker
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 29 April 2012
// ---------------------------------------------------------------------------
// Prints general alignment statistics for BAM file(s).
// ***************************************************************************

#include "statistics.h"

#include <numeric>

using namespace BamTools;
using namespace std;

Statistics::Statistics()
: m_numReads(0)
, m_numPaired(0)
, m_numProperPair(0)
, m_numMapped(0)
, m_numBothMatesMapped(0)
, m_numForwardStrand(0)
, m_numReverseStrand(0)
, m_numFirstMate(0)
, m_numSecondMate(0)
, m_numSingletons(0)
, m_numFailedQC(0)
, m_numDuplicates(0)
, m_showInsertSizeSummary(false)
{ 
    m_insertSizes.reserve(100000);
}

// median is of type double because in the case of even number of data elements,
// we need to return the average of middle 2 elements
bool CalculateMedian(vector<int>& data, double& median) { 
    
    // skip if data empty
    if ( data.empty() ) return false;
    
    // find middle element
    size_t middleIndex = data.size() / 2;
    vector<int>::iterator target = data.begin() + middleIndex;
    nth_element(data.begin(), target, data.end());
    
    // odd number of elements
    if ( (data.size() % 2) != 0) {
        median = (double)(*target);
        return true;
    }
    
    // even number of elements
    else {
        double rightTarget = (double)(*target);
        vector<int>::iterator leftTarget = target - 1;
        nth_element(data.begin(), leftTarget, data.end());
        median = (double)((rightTarget+*leftTarget)/2.0);
        return true;
    }
}

int Statistics::runInternal()
{
    BamAlignment * pal;
    while(NULL != (pal = getInputAlignment())) {
        BamAlignment & al = *pal;

        // increment total alignment counter
        ++m_numReads;
        
        // incrememt counters for pairing-independent flags
        if ( al.IsDuplicate() ) ++m_numDuplicates;
        if ( al.IsFailedQC()  ) ++m_numFailedQC;
        if ( al.IsMapped()    ) ++m_numMapped;
        
        // increment strand counters
        if ( al.IsReverseStrand() ) 
            ++m_numReverseStrand;
        else 
            ++m_numForwardStrand;
        
        // if alignment is paired-end
        if ( al.IsPaired() ) {
            
            // increment PE counter
            ++m_numPaired;
            
            // increment first mate/second mate counters
            if ( al.IsFirstMate()  ) ++m_numFirstMate;
            if ( al.IsSecondMate() ) ++m_numSecondMate;
            
            // if alignment is mapped, check mate status
            if ( al.IsMapped() ) {
                // if mate mapped
                if ( al.IsMateMapped() ) 
                    ++m_numBothMatesMapped;
                // else singleton
                else 
                    ++m_numSingletons;
            }
            
            // check for explicit proper pair flag
            if ( al.IsProperPair() ) ++m_numProperPair;
            
            // store insert size for first mate 
            if ( m_showInsertSizeSummary && al.IsFirstMate() && (al.InsertSize != 0) ) {
                int insertSize = abs(al.InsertSize);
                m_insertSizes.push_back( insertSize );
            }
        }

        putOutputAlignment(pal);
    }
          
    cout << endl;
    cout << "**********************************************" << endl;
    cout << "Stats for BAM file(s): " << endl;
    cout << "**********************************************" << endl;
    cout << endl;
    cout << "Total reads:       " << m_numReads << endl;
    cout << "Mapped reads:      " << m_numMapped << "\t(" << ((float)m_numMapped/m_numReads)*100 << "%)" << endl;
    cout << "Forward strand:    " << m_numForwardStrand << "\t(" << ((float)m_numForwardStrand/m_numReads)*100 << "%)" << endl;
    cout << "Reverse strand:    " << m_numReverseStrand << "\t(" << ((float)m_numReverseStrand/m_numReads)*100 << "%)" << endl;
    cout << "Failed QC:         " << m_numFailedQC << "\t(" << ((float)m_numFailedQC/m_numReads)*100 << "%)" << endl;
    cout << "Duplicates:        " << m_numDuplicates << "\t(" << ((float)m_numDuplicates/m_numReads)*100 << "%)" << endl;
    cout << "Paired-end reads:  " << m_numPaired << "\t(" << ((float)m_numPaired/m_numReads)*100 << "%)" << endl;

    if ( m_numPaired != 0 ) {
      cout << "'Proper-pairs':    " << m_numProperPair << "\t(" << ((float)m_numProperPair/m_numPaired)*100 << "%)" << endl;
      cout << "Both pairs mapped: " << m_numBothMatesMapped << "\t(" << ((float)m_numBothMatesMapped/m_numPaired)*100 << "%)" << endl;
      cout << "Read 1:            " << m_numFirstMate << endl;
      cout << "Read 2:            " << m_numSecondMate << endl;
      cout << "Singletons:        " << m_numSingletons << "\t(" << ((float)m_numSingletons/m_numPaired)*100 << "%)" << endl;
    }

    if ( m_showInsertSizeSummary ) {
      
      double avgInsertSize = 0.0;
      if ( !m_insertSizes.empty() ) {
          avgInsertSize = ( accumulate(m_insertSizes.begin(), m_insertSizes.end(), 0.0) / (double)m_insertSizes.size() );
          cout << "Average insert size (absolute value): " << avgInsertSize << endl;
      }
      
      double medianInsertSize = 0.0;
      if ( CalculateMedian(m_insertSizes, medianInsertSize) )
          cout << "Median insert size (absolute value): " << medianInsertSize << endl;
    }
    cout << endl;
    return 0;
}