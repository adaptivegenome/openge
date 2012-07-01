/*********************************************************************
 *
 * fastq.cpp: A simple FASTQ writer.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 16 March 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include <iostream>
#include <sstream>
#include <cassert>
#include "fastq_writer.h"

using namespace std;
using namespace BamTools;

FastqWriter::FastqWriter() 
: open(false)
{
    
}

bool FastqWriter::Open(const string& filename) {
    this->filename = filename;
    file.open(filename.c_str(), ios::out);
    
    if(file.fail()) {
        cerr << "Failed to open FASTQ output file " << filename << endl;
        return false;
    }

    open = true;
    
    return true;
}

bool FastqWriter::Close() {
    file.close();
    
    open = false;
    
    return true;
}

bool FastqWriter::SaveAlignment(BamTools::BamAlignment & a) {
    assert(open);
    
    
    file << "@" << a.Name << endl;
    file << a.QueryBases << endl;
    file << "+" << a.Name << endl;
    file << a.Qualities << endl;
    
    return true;
}
