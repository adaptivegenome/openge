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
#include <api/BamAlignment.h>
#include <api/algorithms/Sort.h>

class ReadStreamReader {
public:
    virtual bool open(const std::string & filename) = 0;
    virtual const BamTools::SamHeader & getHeader() const = 0;
    virtual void close() = 0;
    virtual BamTools::BamAlignment * read() = 0;
};

template <class T>
class MultiReader : public ReadStreamReader {

    class SortedMergeElement{
    public:
        BamTools::BamAlignment * read;
        T * source;
        bool operator<(const SortedMergeElement & t) const {
            BamTools::Algorithms::Sort::ByPosition cmp = BamTools::Algorithms::Sort::ByPosition();
            return cmp(this->read, t.read);
        }

        SortedMergeElement(BamTools::BamAlignment * read, T * source)
        : read(read)
        , source(source)
        {}
    };

    std::vector<T *> readers;
    std::multiset<SortedMergeElement> reads;
public:
    virtual bool open(const std::string & filename) {
        std::vector<std::string> fn;
        fn.push_back(filename);
        return open(fn);
    }

    virtual bool open(const std::vector<std::string> & filenames) {
        for(std::vector<std::string>::const_iterator i = filenames.begin(); i != filenames.end(); i++) {
            readers.push_back(new T());
            bool ret = readers.back()->open(*i);
            if(!ret) {
                close();
                return false;
            }
        }

        // first, get one read from each queue
        // make sure and deal with the case where one chain will never have any reads. TODO LCB
        
        for(typename std::vector<T *>::iterator i = readers.begin(); i != readers.end(); i++)
        {
            BamTools::BamAlignment * read = (*i)->read();
            
            if(!read)
                continue;

            reads.insert(SortedMergeElement(read, (*i)));
        }
        return true;
    }

    virtual const BamTools::SamHeader & getHeader() const {
        BamTools::SamSequenceDictionary s = readers.front()->getHeader().Sequences;
        for(typename std::vector<T *>::const_iterator i = readers.begin(); i != readers.end(); i++) {
            if(s != (*i)->getHeader().Sequences) {
                std::cerr << "Warning; sequence headers vary between files. Data may be corrupt." << std::endl;
            }
        }
        return readers.front()->getHeader();
    }

    virtual void close() {
        for(typename std::vector<T *>::iterator i = readers.begin(); i != readers.end(); i++) {
            (*i)->close();
        }
    }

    virtual BamTools::BamAlignment * read() {
        BamTools::BamAlignment * ret = NULL;
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