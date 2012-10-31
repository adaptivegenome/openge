#ifndef OGE_BGZF_INPUT_STREAM_H
#define OGE_BGZF_INPUT_STREAM_H

/*********************************************************************
 *
 * bgzf_input_stream.h: Use BGZF compression on a stream.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 29 August 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial
 * Purpose License. A copy of this license has been provided in
 * the openge/ directory.
 *
 *********************************************************************/
#include <fstream>
#include <map>
#include <vector>
#include "thread_pool.h"

#include <iostream>

class BgzfInputStream
{
    
    class BgzfCacheElement : public ThreadJob {
    public:
        std::vector<char> compressed_data;
        std::vector<char> uncompressed_data;
        size_t file_position;
        bool loaded;
        
        size_t next_block_address() { return file_position + compressed_data.size(); }
        virtual void runJob();
        BgzfCacheElement()
        : loaded(false)
        {}
    };
public:
    bool open(std::string filename);
    void read(char * data, size_t len);
    void close();
    bool is_open() { return *input_stream == std::cin || input_stream_real.is_open(); }
    bool fail() { return cache.empty() && input_stream->fail(); }
    bool eof() { return cache.empty() && input_stream->eof(); }
protected:
    size_t current_block;
    size_t current_offset;
    std::istream * input_stream;
    std::ifstream input_stream_real;
    std::map<size_t, BgzfCacheElement *> cache;
    int cached_blocks_read, cache_misses;

    bool requestNextBlock();
};

#endif