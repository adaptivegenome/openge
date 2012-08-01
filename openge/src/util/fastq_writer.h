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

#include "file_io.h"
#include <string>
#include <iostream>
#include <api/BamAlignment.h>
#include <api/SamHeader.h>

#include <map>
#include <string>

class FastqWriter : public FileWriterClass
{
public:
    FastqWriter();
    bool Open(const std::string& filename, const BamTools::SamHeader & samHeader);
    bool Close();
    bool SaveAlignment( BamTools::BamAlignment & al);
protected:
    typedef struct {
        std::string qual, seq;
    } fastq_record_t;

    std::ofstream fwd_file, rev_file, orphan_file;
    std::ostream * fwd_stream, * rev_stream, * orphan_stream;
    BamTools::SamHeader header;
    BamTools::RefVector references;
    std::string filename;
    std::map<std::string, fastq_record_t> potential_pairs;

    bool open;
};

#endif
