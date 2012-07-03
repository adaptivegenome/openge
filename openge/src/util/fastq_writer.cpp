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
: output_stream(&cout)
, open(false)
{
    
}

bool FastqWriter::Open(const string& filename, const std::string& samHeaderText, const BamTools::RefVector& referenceSequences) {
    this->filename = filename;
    
    if(filename != "stdout") {
        file.open(filename.c_str(), ios::out);
        
        if(file.fail()) {
            cerr << "Failed to open FASTQ output file " << filename << endl;
            return false;
        }
        output_stream = &file;
    }

    open = true;
    
    return true;
}

bool FastqWriter::Close() {
    if(file.is_open())
        file.close();
    
    open = false;
    
    return true;
}

bool FastqWriter::SaveAlignment(BamTools::BamAlignment & a) {
    assert(open);
    
    a.BuildCharData();
    
    *output_stream << "@" << a.Name << endl;
    *output_stream << a.QueryBases << endl;
    *output_stream << "+" << a.Name << endl;
    *output_stream << a.Qualities << endl;
    
    return true;
}
