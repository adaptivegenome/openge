/*********************************************************************
 *
 * file_reader.cpp:  Algorithm module that filters a stream of reads.
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
 * File reader algorithm module. Provides an abstracted interface to 
 * convert a BAM or SAM file into a stream of reads for other 
 * algorithm module subclasses.
 *
 *********************************************************************/

#include "file_reader.h"

#include "api/BamMultiReader.h"
#include "../util/read_stream_reader.h"

using namespace BamTools;
using namespace std;

int FileReader::runInternal()
{
    ogeNameThread("am_FileReader");

    {
        MultiReader reader;

        if(!reader.open(filenames)) {
            cerr << "Error opening one or more files." << endl;
            reader.close();
            exit(-1);
        }

        header = reader.getHeader();
        open = true;
        
        BamAlignment * al;
        
        while(true)
        {
            al = reader.read();
            
            if(!al)
                break;
            
            putOutputAlignment(al);
        }
        
        reader.close();
    }

    return 0;
}

void FileReader::addFile(std::string filename) 
{
    filenames.push_back(filename); 
}

void FileReader::addFiles(std::vector<std::string> filenames) 
{
    this->filenames.insert(this->filenames.end(), filenames.begin(), filenames.end()); 
}