#ifndef OGE_FILE_IO_H
#define OGE_FILE_IO_H

/*********************************************************************
 *
 * file_io.h:  Constants shared by file reading and writing.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 4 July 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include <api/BamAlignment.h>

#include <string>
#include <algorithm>

typedef enum
{
    FORMAT_BAM, FORMAT_SAM, FORMAT_CRAM, FORMAT_FASTQ, FORMAT_UNKNOWN
} file_format_t;

inline file_format_t detectFileFormatFromFilename(std::string filename) {
    std::transform(filename.begin(), filename.end(), filename.begin(), ::tolower);
    
    if(filename.size() >= 3 && filename.substr(filename.size() - 3, 3) == "bam")
        return FORMAT_BAM;
    if(filename.size() >= 3 && filename.substr(filename.size() - 3, 3) == "sam")
        return FORMAT_SAM;
    if(filename.size() >= 4 && filename.substr(filename.size() - 4, 4) == "cram")
        return FORMAT_CRAM;
    if(filename.size() >= 5 && filename.substr(filename.size() - 5, 5) == "fastq")
        return FORMAT_FASTQ;
    
    return FORMAT_UNKNOWN;
}

class FileWriterClass {
public:
    virtual bool Open(const std::string& filename, const std::string& samHeaderText, const BamTools::RefVector& referenceSequences) = 0;
    virtual bool SaveAlignment(BamTools::BamAlignment & a) = 0;
};

class FileReaderClass {
    
};

#endif