#ifndef FILE_READER_H
#define FILE_READER_H

#include "algorithm_module.h"

#include <vector>
#include <string>

class FileReader : public AlgorithmModule
{
    
protected:
    virtual int runInternal();
    std::vector<std::string> filenames;
    bool open;
    BamTools::SamHeader header;
    BamTools::RefVector references;
    
    virtual BamTools::RefVector getReferences() { while (!open) usleep(10000); return references; }
    virtual BamTools::SamHeader getHeader() { while (!open) usleep(10000); return header; }
    
public:
    FileReader() : open(false) {}
    void addFile(std::string filename);
    void addFiles(std::vector<std::string> filenames);
};

#endif