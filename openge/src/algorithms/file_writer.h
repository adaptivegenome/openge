#ifndef OGE_ALGO_FILE_WRITER_H
#define OGE_ALGO_FILE_WRITER_H

/*********************************************************************
 *
 * file_writer.h:  Algorithm module writes a BAM or SAM file.
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

#include "algorithm_module.h"

#include <vector>
#include <string>

class FileWriter : public AlgorithmModule
{
    
protected:
    virtual int runInternal();
    std::string filename;
    int compression_level;
    
public:
    FileWriter() : compression_level(6) {}
    void setFilename(std::string filename) { this->filename = filename; }
    void setCompressionLevel(int level) { compression_level = level;}
    size_t getCount() { return write_count; }
    
// eventually we should also be able to set the format.
// void setFormat(???);
};
#endif