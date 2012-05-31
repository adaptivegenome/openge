#ifndef OGE_ALGO_FILE_READER_H
#define OGE_ALGO_FILE_READER_H

/*********************************************************************
 *
 * file_reader.h:  Algorithm module subclass that filters a stream of reads.
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

#include "algorithm_module.h"

#include <vector>
#include <string>

class FileReader : public AlgorithmModule
{
public:
    typedef enum
    {
        FORMAT_BAM, FORMAT_SAM, FORMAT_CRAM, FORMAT_UNKNOWN
    } file_format_t;
    
protected:
    virtual int runInternal();
    std::vector<std::string> filenames;
    bool open;
    BamTools::SamHeader header;
    BamTools::RefVector references;
    file_format_t format;
    bool format_specified;
    bool load_string_data;
    
    virtual BamTools::RefVector getReferences() { while (!open) usleep(10000); return references; }
    virtual BamTools::SamHeader getHeader() { while (!open) usleep(10000); return header; }
    
public:
    FileReader() : open(false), format(FORMAT_BAM), format_specified(false), load_string_data(true) {}
    void addFile(std::string filename);
    void addFiles(std::vector<std::string> filenames);
    void setFormat(file_format_t f) { format = f; format_specified = true;}
    file_format_t getFormat() { return format; }
    size_t getCount() { return write_count; }
    void setLoadStringData(bool string_data) { load_string_data = string_data; }
    bool getLoadStringData() { return load_string_data; }
    
    file_format_t deduceFileFormat();
};

#endif