/*********************************************************************
 *
 * sam_writer.h: A simple SAM writer.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 6 June 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#ifndef OGE_SAMWRITER_H
#define OGE_SAMWRITER_H
#include <string>
#include <iostream>
#include <fstream>
#include "bam_header.h"
#include "read_stream_writer.h"

// SamReader is capable of sequentially reading a SAM file. It doesn't support
// most of the features that BamReader does, only enough to support converting SAM
// files to the BAM format.
class SamWriter : public ReadStreamWriter
{
public:
    SamWriter();
    bool open(const std::string& filename,
              const BamHeader & samHeader);
    void close();
    bool write( const OGERead & al);
    bool is_open() const { return m_open; }
protected:
    std::ofstream file;
    std::ostream * output_stream;
    BamHeader header;
    std::string filename;

    bool m_open;
};

#endif
