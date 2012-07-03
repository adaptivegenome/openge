/*********************************************************************
 *
 * fastq_writer.h: A simple SAM writer.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 1 July 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#ifndef OGE_FASTQWRITER_H
#define OGE_FASTQWRITER_H
#include <string>
#include <iostream>
#include <api/BamAlignment.h>
#include <api/SamHeader.h>

class FastqWriter
{
public:
    FastqWriter();
    bool Open(const std::string& filename);
    bool Close();
    bool SaveAlignment( BamTools::BamAlignment & al);
protected:
    std::ofstream file;
    std::ostream * output_stream;
    BamTools::SamHeader header;
    BamTools::RefVector references;
    std::string filename;

    bool open;
};

#endif
