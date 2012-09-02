/*********************************************************************
 *
 * read_stream_reader.h: Interface for reading a BAM/SAM/CRAM file.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 28 August 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial
 * Purpose License. A copy of this license has been provided in
 * the openge/ directory.
 *
 *********************************************************************/

#include "read_stream_reader.h"

#include "sam_reader.h"
#include "bgzf_input_stream.h"
#include "bam_deserializer.h"

bool MultiReader::open(const std::vector<std::string> & filenames) {
    for(std::vector<std::string>::const_iterator i = filenames.begin(); i != filenames.end(); i++) {
        ReadStreamReader * reader = NULL;
        
        switch (ReadStreamReader::detectFileFormat(*i)) {
            case FORMAT_BAM: reader = new BamDeserializer<BgzfInputStream>(); break;
            case FORMAT_RAWBAM: reader = new BamDeserializer<std::ifstream>(); break;
            case FORMAT_SAM: reader = new ::SamReader(); break;
            default:
                std::cerr << "File " << *i << " is of an unknown format. Aborting." << std::endl;
                exit(-1);
                break;
        }
        readers.push_back(reader);
        bool ret = readers.back()->open(*i);
        if(!ret) {
            std::cerr << "Failed to open " << *i << std::endl;
            close();
            return false;
        }
    }
    
    // first, get one read from each queue
    // make sure and deal with the case where one chain will never have any reads. TODO LCB
    
    for(std::vector<ReadStreamReader *>::iterator i = readers.begin(); i != readers.end(); i++)
    {
        BamTools::BamAlignment * read = (*i)->read();
        
        if(!read)
            continue;
        
        reads.insert(SortedMergeElement(read, (*i)));
    }
    return true;
}