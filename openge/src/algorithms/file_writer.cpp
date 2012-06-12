/*********************************************************************
 *
 * file_writer.cpp:  Algorithm module writes a BAM or SAM file.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 31 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************
 *
 * File writer algorithm module. Writes a stream of reads to a BAM
 * file. Eventually this will be extended to support SAM or CRAM
 * formats.
 *
 *********************************************************************/

#include "file_writer.h"

#include "../util/sam_writer.h"

#include "api/BamWriter.h"
using namespace BamTools;
using namespace std;

#ifdef __linux__
#include <sys/prctl.h>
#endif


//from http://stackoverflow.com/questions/874134/find-if-string-endswith-another-string-in-c
bool hasEnding (std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

int FileWriter::runInternal()
{
#ifdef __linux__
    prctl(PR_SET_NAME,"am_FileWriter",0,0,0);
#endif
    
    if(filename.size() > 3 && (hasEnding(filename, "sam") || hasEnding(filename, "SAM"))) {
        SamWriter writer;
        
        if(!writer.Open(filename, getHeader().ToString(), getReferences())) {
            cerr << "Error opening BAM file to write." << endl;
            return -1;
        }

        BamAlignment * al;
        
        while(true)
        {
            al = getInputAlignment();
            if(!al)
                break;
            
            writer.SaveAlignment(*al);
            putOutputAlignment(al);
        }
        
        writer.Close();
        
    } else {
        BamWriter writer;
        
        if(!writer.Open(filename, getHeader(), getReferences())) {
            cerr << "Error opening BAM file to write." << endl;
            return -1;
        }
        
        writer.SetCompressionLevel(compression_level);
        
        BamAlignment * al;
        
        while(true)
        {
            al = getInputAlignment();
            if(!al)
                break;

            writer.SaveAlignment(*al);
            putOutputAlignment(al);
        }
        
        writer.Close();
    }
    
    if(isVerbose())
        cerr << "Wrote " << write_count << " reads to " << filename << endl;
    
    return 0;
}