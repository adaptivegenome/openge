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
#include "../util/file_io.h"

class FileWriter : public AlgorithmModule
{
    
protected:
    virtual int runInternal();
    std::string filename;
    int compression_level;
    file_format_t file_format, default_file_format;
    std::string command_line_options;
public:
    FileWriter() : compression_level(6), file_format(FORMAT_UNKNOWN), default_file_format(FORMAT_BAM) {}
    void setFilename(std::string filename) { this->filename = filename; }
    void setCompressionLevel(int level) { compression_level = level;}
    size_t getCount() { return write_count; }
    void setFormat(file_format_t format) { file_format = format; }   //force this format to be used
    void setFormat(const std::string & format_name);
    void setDefaultFormat(file_format_t format) { default_file_format = format; }    //format to be used if auto detection doesn't work
    void addProgramLine(const std::string & command_options) { this->command_line_options = command_options; }
    file_format_t getFileFormat();
};
#endif