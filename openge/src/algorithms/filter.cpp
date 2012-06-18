/*********************************************************************
 *
 * filter.cpp:  Algorithm module that filters a stream of reads.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 22 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************
 *
 * Filter algorithm module class. Filters a stream of reads by the 
 * region of the genome, or by count, or both. Some code based on 
 * Bamtools' merge code, from bamtools_merge.cpp.
 *
 *********************************************************************/

#include "filter.h"

#include <algorithm>
#include <numeric>

using namespace BamTools;
using namespace std;


// this has been copied from bamtools utilities, since it isn't in the API. Original file is bamtools_utilities.cpp.
// Like the rest of Bamtools, it is under the BSD license.
bool Filter::ParseRegionString(const string& regionString, BamRegion& region)
{
    // -------------------------------
    // parse region string
    
    // check first for empty string
    if ( regionString.empty() ) 
        return false;   
    
    // non-empty string, look for a colom
    size_t foundFirstColon = regionString.find(':');
    
    // store chrom strings, and numeric positions
    string chrom;
    int startPos;
    int stopPos;
    
    // no colon found
    // going to use entire contents of requested chromosome 
    // just store entire region string as startChrom name
    // use BamReader methods to check if its valid for current BAM file
    if ( foundFirstColon == string::npos ) {
        chrom = regionString;
        startPos   = 0;
        stopPos    = -1;
    }
    
    // colon found, so we at least have some sort of startPos requested
    else {
        
        // store start chrom from beginning to first colon
        chrom = regionString.substr(0,foundFirstColon);
        
        // look for ".." after the colon
        size_t foundRangeDots = regionString.find("..", foundFirstColon+1);
        
        // no dots found
        // so we have a startPos but no range
        // store contents before colon as startChrom, after as startPos
        if ( foundRangeDots == string::npos ) {
            startPos   = atoi( regionString.substr(foundFirstColon+1).c_str() ); 
            stopPos    = -1;
        } 
        
        // ".." found, so we have some sort of range selected
        else {
            
            // store startPos between first colon and range dots ".."
            startPos = atoi( regionString.substr(foundFirstColon+1, foundRangeDots-foundFirstColon-1).c_str() );
            
            // look for second colon
            size_t foundSecondColon = regionString.find(':', foundRangeDots+1);
            
            // no second colon found
            // so we have a "standard" chrom:start..stop input format (on single chrom)
            if ( foundSecondColon == string::npos ) {
                stopPos    = atoi( regionString.substr(foundRangeDots+2).c_str() );
            } else {
                return false;
            }
        }
    }
    
    // -------------------------------
    // validate reference IDs & genomic positions
    
    const RefVector references = getReferences();
    
    int RefID = -1;
    for(int i = 0; i < references.size(); i++) {
        if(references[i].RefName == chrom)
            RefID = i;
    }
    
    // if startRefID not found, return false
    if ( RefID == -1 ) {
        cerr << "Can't find chromosome'" << chrom << "'" << endl;
        return false;
    }
    
    // startPos cannot be greater than or equal to reference length
    const RefData& startReference = references.at(RefID);
    if ( startPos >= startReference.RefLength ) {
        cerr << "Start position (" << startPos << ") after end of the reference sequence (" << startReference.RefLength << ")" << endl;
        return false;
    }
    
    // stopPosition cannot be larger than reference length
    const RefData& stopReference = references.at(RefID);
    if ( stopPos > stopReference.RefLength ) {
        cerr << "Start position (" << stopPos << ") after end of the reference sequence (" << stopReference.RefLength << ")" << endl;
        return false;
    }

    // if no stopPosition specified, set to reference end
    if ( stopPos == -1 ) stopPos = stopReference.RefLength;
    
    // -------------------------------
    // set up Region struct & return
    
    region.LeftRefID     = RefID;
    region.LeftPosition  = startPos;
    region.RightRefID    = RefID;;
    region.RightPosition = stopPos;
    return true;
}

Filter::Filter()
: has_region(false)
, count_limit(INT_MAX)
{}

int Filter::runInternal()
{
    size_t count = 0;

    // if no region specified, store entire contents of file(s)
    if ( !has_region ) {
        BamAlignment * al = NULL;
        while (NULL != (al = getInputAlignment()) && count < count_limit ) {
            putOutputAlignment(al);
            count++;
        }
        
        //flush rest of queue
        while (NULL != (al = getInputAlignment()) ) {
            delete al;
        }
    }
    
    // otherwise attempt to use region as constraint
    else {
        if(isVerbose())
            cerr << "Filtering to region " << region_string << endl;
        
        // if region string parses OK
        BamRegion region;
        if (ParseRegionString(region_string, region) ) {
            BamAlignment * al;
            while ( NULL != (al = getInputAlignment())  && count < count_limit) {
                if ( (al->RefID >= region.LeftRefID)  && ( (al->Position + al->Length) >= region.LeftPosition ) &&
                    (al->RefID <= region.RightRefID) && ( al->Position <= region.RightPosition) ) {
                    putOutputAlignment(al);
                    count++;
                } else {
                    delete al;
                }
            }
        }
        
        // error parsing REGION string
        else {
            cerr << "ERROR: could not parse region'" << region_string << "'" << endl;
            cerr << "Check that region description is in valid format (see documentation) and that the coordinates are valid"
            << endl;
            exit(-1);
            return -1;
        }
    }
    
    if(isVerbose())
        cerr << count << " alignments processed." << endl;
    
    return 0;
}