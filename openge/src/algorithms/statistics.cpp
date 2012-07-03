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

#include "statistics.h"

#include <algorithm>
#include <numeric>

#include <iomanip>

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
    
    map<int, int> read_len_ct;

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
        
        if(read_len_ct.end() != read_len_ct.find(al.Length))
            read_len_ct[al.Length]++;
        else
            read_len_ct[al.Length] = 1;

        putOutputAlignment(pal);
    }

    const int precision = 1;
    const int field_width = 5;
    cout << "Total reads:       " << setw(10) << m_numReads << endl;
    cout << "Mapped reads:      " << setw(10) << m_numMapped << " (" << setprecision(precision) << setw(field_width) << fixed << ((float)m_numMapped/m_numReads)*100 << "%)" << endl;
    cout << "Forward strand:    " << setw(10) << m_numForwardStrand << " (" << setprecision(precision) << setw(field_width) << fixed << ((float)m_numForwardStrand/m_numReads)*100 << "%)" << endl;
    cout << "Reverse strand:    " << setw(10) << m_numReverseStrand << " (" << setprecision(precision) << setw(field_width) << fixed << ((float)m_numReverseStrand/m_numReads)*100 << "%)" << endl;
    cout << "Failed QC:         " << setw(10) << m_numFailedQC << " (" << setprecision(precision) << setw(field_width) << fixed << ((float)m_numFailedQC/m_numReads)*100 << "%)" << endl;
    cout << "Duplicates:        " << setw(10) << m_numDuplicates << " (" << setprecision(precision) << setw(field_width) << fixed << ((float)m_numDuplicates/m_numReads)*100 << "%)" << endl;
    cout << "Paired-end reads:  " << setw(10) << m_numPaired << " (" << setprecision(precision) << setw(field_width) << fixed << ((float)m_numPaired/m_numReads)*100 << "%)" << endl;

    if ( m_numPaired != 0 ) {
      cout << "'Proper-pairs':    " << setw(10) << m_numProperPair << " (" << setprecision(precision) << setw(field_width) << fixed << ((float)m_numProperPair/m_numPaired)*100 << "%)" << endl;
      cout << "Both pairs mapped: " << setw(10) << m_numBothMatesMapped << " (" << setprecision(precision) << setw(field_width) << fixed << ((float)m_numBothMatesMapped/m_numPaired)*100 << "%)" << endl;
      cout << "Read 1:            " << setw(10) << m_numFirstMate << endl;
      cout << "Read 2:            " << setw(10) << m_numSecondMate << endl;
      cout << "Singletons:        " << setw(10) << m_numSingletons << " (" << setprecision(precision) << setw(field_width) << fixed << ((float)m_numSingletons/m_numPaired)*100 << "%)" << endl;
    }
    
    if (m_showLengthSummary) {
        cerr << "Read lengths:" << endl;

        for(map<int,int>::const_iterator it = read_len_ct.begin(); it != read_len_ct.end(); it++)
        {
            float pct = 100. * (double) it->second / (double) getReadCount();
            
            cout << " " << setw(5) << it->first << "bp:          " << setw(10) << it->second << " (" << setprecision(precision) << setw(field_width) << fixed << pct << "%)" << endl;
        }
    }

    cout << resetiosflags(ios_base::adjustfield) << resetiosflags(ios_base::floatfield);

    if ( m_showInsertSizeSummary ) {
        cout << "Insert size (absolute value):" << endl;
        double avgInsertSize = 0.0;
        if ( !m_insertSizes.empty() ) {
            avgInsertSize = ( accumulate(m_insertSizes.begin(), m_insertSizes.end(), 0.0) / (double)m_insertSizes.size() );
            cout << "    Mean:          " <<  setprecision(1) << setw(10) << fixed << avgInsertSize << endl;
        }
        
        double medianInsertSize = 0.0;
        if ( CalculateMedian(m_insertSizes, medianInsertSize) )
            cout << "    Median:        " <<  setprecision(1) << setw(10) << fixed << medianInsertSize << endl;
    }

    cout << resetiosflags(ios_base::adjustfield) << resetiosflags(ios_base::floatfield);
    return 0;
}