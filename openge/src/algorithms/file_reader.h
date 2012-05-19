#ifndef OGE_ALGO_FILE_READER_H
#define OGE_ALGO_FILE_READER_H

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
    size_t count;
    bool load_string_data;
    
    virtual BamTools::RefVector getReferences() { while (!open) usleep(10000); return references; }
    virtual BamTools::SamHeader getHeader() { while (!open) usleep(10000); return header; }
    
public:
    FileReader() : open(false), format(FORMAT_BAM), format_specified(false), count(0), load_string_data(true) {}
    void addFile(std::string filename);
    void addFiles(std::vector<std::string> filenames);
    void setFormat(file_format_t f) { format = f; format_specified = true;}
    file_format_t getFormat() { return format; }
    size_t getCount() { return count; }
    void setLoadStringData(bool string_data) { load_string_data = string_data; }
    bool getLoadStringData() { return load_string_data; }
    
    file_format_t deduceFileFormat();
};

#endif