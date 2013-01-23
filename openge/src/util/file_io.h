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

#include "oge_read.h"

#include <string>
#include <algorithm>

typedef enum
{
    FORMAT_BAM, FORMAT_RAWBAM, FORMAT_SAM, FORMAT_CRAM, FORMAT_FASTQ, FORMAT_UNKNOWN
} file_format_t;

inline file_format_t detectFileFormatFromFilename(std::string filename) {
    std::transform(filename.begin(), filename.end(), filename.begin(), ::tolower);
    
    if(filename.size() >= 6 && filename.substr(filename.size() - 6, 6) == "rawbam")
        return FORMAT_RAWBAM;
    if(filename.size() >= 6 && filename.substr(filename.size() - 6, 6) == "bamraw")
        return FORMAT_RAWBAM;
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

#endif