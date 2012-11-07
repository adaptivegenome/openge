#ifndef OGE_READ_STREAM_READER_H
#define OGE_READ_STREAM_READER_H

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

#include <api/SamHeader.h>
#include <api/algorithms/Sort.h>
#include "oge_read.h"

#include <iostream>

#include <set>

class ReadStreamReader {
public:
    typedef enum
    {
        FORMAT_BAM, FORMAT_RAWBAM, FORMAT_SAM, FORMAT_CRAM, FORMAT_UNKNOWN
    } file_format_t;

    virtual bool open(const std::string & filename) = 0;
    virtual const BamTools::SamHeader & getHeader() const = 0;
    virtual void close() = 0;
    virtual OGERead * read() = 0;
    
    static inline file_format_t detectFileFormat(std::string filename);
};

ReadStreamReader::file_format_t ReadStreamReader::detectFileFormat(std::string filename) {
    file_format_t this_file_format = FORMAT_UNKNOWN;
    FILE * fp = NULL;
    if(filename == "stdin")
        fp = stdin;
    else {
        fp = fopen(filename.c_str(), "rb");
        
        if(!fp) {
            perror("Read stream reader can't open a file");
            std::cerr << "Couldn't open file " << filename << std::endl;
            return FORMAT_UNKNOWN;
        }
    }
    
    unsigned char data[2];
    if(2 != fread(data, 1,2, fp)) {
        perror("Read stream reader can't read file");
        std::cerr << "Couldn't read from file " << filename << std::endl;
        return FORMAT_UNKNOWN;
    }
    
    if(filename != "stdin")
        fclose(fp);
    else {
        int ret = ungetc(data[1], fp);
        assert(ret != EOF);
        ret = ungetc(data[0], fp);
        assert(ret != EOF);
    }
    
    if(data[0] == '@')
        this_file_format = FORMAT_SAM;
    else if(data[0] == 31 && data[1] == 139)
        this_file_format = FORMAT_BAM;
    else if(data[0] == 'B' && data[1] == 'A')
        this_file_format = FORMAT_RAWBAM;
    else
        this_file_format = FORMAT_UNKNOWN;
    
    return this_file_format;
}

class MultiReader : public ReadStreamReader {

    class SortedMergeElement{
    public:
        OGERead * read;
        ReadStreamReader * source;
        bool operator<(const SortedMergeElement & t) const {
            BamTools::Algorithms::Sort::ByPosition cmp = BamTools::Algorithms::Sort::ByPosition();
            return cmp(this->read, t.read);
        }

        SortedMergeElement(OGERead * read, ReadStreamReader * source)
        : read(read)
        , source(source)
        {}
    };

    std::vector<ReadStreamReader *> readers;
    std::multiset<SortedMergeElement> reads;
public:
    virtual bool open(const std::string & filename) {
        std::vector<std::string> fn;
        fn.push_back(filename);
        return open(fn);
    }

    virtual bool open(const std::vector<std::string> & filenames);

    virtual const BamTools::SamHeader & getHeader() const {
        BamTools::SamSequenceDictionary s = readers.front()->getHeader().Sequences;
        for( std::vector<ReadStreamReader *>::const_iterator i = readers.begin(); i != readers.end(); i++) {
            if(s != (*i)->getHeader().Sequences) {
                std::cerr << "Warning; sequence headers vary between files. Data may be corrupt." << std::endl;
            }
        }
        return readers.front()->getHeader();
    }

    virtual void close() {
        for( std::vector<ReadStreamReader *>::iterator i = readers.begin(); i != readers.end(); i++) {
            (*i)->close();
        }
    }

    virtual OGERead * read() {
        OGERead * ret = NULL;
        //now handle the steady state situation. When sources are done, We
        // won't have a read any more in the reads pqueue.
        while(!ret && !reads.empty()) {
            SortedMergeElement el = *reads.begin();
            reads.erase(reads.begin());
            
            ret = el.read;
            
            el.read = el.source->read();
            if(!el.read)
                continue;
            
            reads.insert(el);
        }

        return ret;
    }
};

#endif