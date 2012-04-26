
#ifndef FILE_WRITER_H
#define FILE_WRITER_H

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
    
// eventually we should also be able to set the format.
// void setFormat(???);
};
#endif