/*********************************************************************
 *
 * bgzf_input_stream.cpp: Use BGZF compression on a stream.
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

#include "api/api_global.h"
#include "bgzf_input_stream.h"
#include <zlib.h>

#include <iostream>
#include <cassert>
#include <cstring>

using namespace std;

const int NUM_BLOCKS_IN_CACHE = 200;    //tune this parameter to balance RAM usage vs. speed. one block is 64k of RAM
const int BGZF_BLOCK_SIZE = 65536;

void BgzfInputStream::BgzfCacheElement::runJob() {
    ogeNameThread("ogeBgzfInflate");
    if(compressed_data[0] != 31 || compressed_data[1] != (char)139) {
        cerr << "Error- BGZF block has invalid start block. Is this file corrupted?" << endl;
        exit(-1);
    }

    if(compressed_data[2] != 8 || compressed_data[3] != 4) {
        cerr << "Error- BGZF block has unexpected flags. Is this file corrupted?" << endl;
        exit(-1);
    }
    
    //uint32_t timestamp = *((uint32_t *) &compressed_data[4]);
    //uint8_t extra_flags = *((uint8_t *) &compressed_data[8]);
    //uint8_t os = *((uint8_t *) &compressed_data[9]);
    uint16_t extra_length = *((uint16_t *) &compressed_data[10]);

    char si1 = compressed_data[12];
    char si2 = compressed_data[13];
    //uint16_t slen = *((uint16_t *) &compressed_data[14]);
    uint32_t bsize = *((uint16_t *) &compressed_data[16]) + 1;

    if(extra_length != 6 || si1 != 66 || si2 != 67) {
        cerr << "Error- BGZF GZ extra field is incorrect. Is this file corrupted?" << endl;
        exit(-1);
    }
    
    char * data = &compressed_data[18];
    
    size_t uncompressed_position = bsize - 4;
    size_t uncompressed_size = *((uint32_t *) &compressed_data[uncompressed_position]);
    uncompressed_data.resize(BGZF_BLOCK_SIZE);

    {
        // setup zlib stream object
        z_stream zs;
        zs.zalloc    = NULL;
        zs.zfree     = NULL;
        zs.next_in   = (unsigned char *)data;
        zs.avail_in  = bsize - 16;
        zs.next_out  =  (unsigned char *) &uncompressed_data[0];
        zs.avail_out = uncompressed_data.size();

        // initialize
        int status = inflateInit2(&zs, -15);
        if ( status != Z_OK ) {
            cerr << "Error- Zlib initialization failed. Aborting." << endl;
            exit(-1);
        }
        
        // decompress
        status = inflate(&zs, Z_FINISH);
        if ( status != Z_STREAM_END ) {
            inflateEnd(&zs);
            cerr << "Error- Zlib inflate failed. Aborting." << endl;
            exit(-1);
        }
        
        // finalize
        status = inflateEnd(&zs);
        if ( status != Z_OK ) {
            inflateEnd(&zs);
            cerr << "Error- Zlib inflateEnd failed. Aborting." << endl;
            exit(-1);
        }

        assert(zs.total_out == uncompressed_size);
        uncompressed_data.resize(zs.total_out);
    }
    
    loaded = true;
}

bool BgzfInputStream::requestNextBlock() {

    if(input_stream->eof())
        return false;

    //find last element from file, so we can get the address for the next one.
    BgzfCacheElement * last_el = cache.empty() ? NULL : cache.begin()->second;
    BgzfCacheElement * new_block = new BgzfCacheElement();

    for(map<size_t, BgzfCacheElement *>::const_iterator i = cache.begin(); i != cache.end(); i++) {
        if(i->second->file_position > last_el->file_position)
            last_el = i->second;
    }

    size_t start_offset = input_stream->tellg();
    new_block->file_position = start_offset;
    
    //read the block header to get actual block size
    new_block->compressed_data.resize(18);
    input_stream->read(&new_block->compressed_data[0], 18);
    
    //if we fail at reads, and have read at least part of a block
    if(input_stream->fail()) {
        delete new_block;
        if(input_stream->gcount() != 0) {
            cerr << "BGZF file seems to be truncated. Is this file corrupt?" << endl;
            exit(-1);
        }
        return false;
    }

    uint32_t size = *((uint16_t*) &new_block->compressed_data[16]) + 1;

    new_block->file_position = last_el ? last_el->file_position + last_el->compressed_data.size() : 0;
    new_block->compressed_data.resize(size);
    input_stream->read(&new_block->compressed_data[18], size-18);
    
    if( input_stream->fail()) {
        delete new_block;
        cerr << "BGZF file seems to be truncated. Is this file corrupt?" << endl;
        exit(-1);
        return false;
    }

    cache[new_block->file_position] = new_block;
    ThreadPool::sharedPool()->addJob(new_block);
    return true;
}

bool BgzfInputStream::open(string filename) {
    if(filename == "stdin") {
        input_stream = &cin;
    } else {
        ifstream_buffer.resize(4000000);
        input_stream_real.rdbuf()->pubsetbuf(&ifstream_buffer[0], ifstream_buffer.size());
        input_stream = &input_stream_real;
        input_stream_real.open(filename.c_str());
    
        if(input_stream->fail()) {
            cerr << "BGZF open() failed for " << filename << endl;
            return false;
        }
    }

    for(int i = 0; i < NUM_BLOCKS_IN_CACHE; i++)
        requestNextBlock();
    
    current_offset = 0;
    current_block = 0;
    cached_blocks_read = 0;
    cache_misses = 0;

    return true;
}

void BgzfInputStream::read(char * data, size_t len) {
    size_t copied_data_len = 0;
    while(copied_data_len != len) {
        BgzfCacheElement * current_block_element = cache[current_block];
        cached_blocks_read++;
        
        if(!current_block_element->loaded)
            cache_misses++;

        //wait for element to be loaded if it isn't
        while(!current_block_element->loaded)
            usleep(10000);
        
        //handle moving to next cache block if needed.
        if(current_offset == current_block_element->uncompressed_data.size()) {
            requestNextBlock();
            
            cache.erase(current_block);

            //set up new block
            current_block += current_block_element->compressed_data.size();
            current_offset = 0;
            
            // remove existing block
            delete current_block_element;
            
            if(!cache.count(current_block)) {
                return;
            }
            current_block_element = cache[current_block];

            if(!current_block_element->loaded)
                cache_misses++;
            
            //wait for element to be loaded if it isn't
            while(!current_block_element->loaded)
                usleep(10000);
        }
        
        size_t copy_len = min(len - copied_data_len, current_block_element->uncompressed_data.size() - current_offset);
        
        memcpy(&data[copied_data_len], &current_block_element->uncompressed_data[current_offset], copy_len);
        current_offset += copy_len;
        
        copied_data_len += copy_len;
    }
}

void BgzfInputStream::close() {
    if(input_stream_real.is_open())
        input_stream_real.close();
    
    // use this line to see if NUM_BLOCKS_IN_CACHE is too small.
    //if(cache_misses)
    //    cerr << cache_misses << " cache misses in reading file (of " << cached_blocks_read << " blocks read)." << endl;

    for(map<size_t, BgzfCacheElement *>::const_iterator i = cache.begin(); i != cache.end(); i++)
        delete i->second;
}
