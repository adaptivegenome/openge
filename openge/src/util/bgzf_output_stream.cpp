/*********************************************************************
 *
 * bgzf_output_stream.cpp: Use BGZF compression on a stream.
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
 *********************************************************************
 *
 * This module uses both threadpools and single write thread to
 * perform work. When data is received (via write()), compress jobs
 * are queued up to the thread pool if necessary. When a compress 
 * job finishes, it signals the single write thread to write any
 * new data to a file. This ensures that only one thread accesses
 * the file at any time.
 *
 ********************************************************************/

#include "bgzf_output_stream.h"
#include <iostream>
#include <zlib.h>
using namespace std;

//////////////////
// BgzfBlock class

void BgzfOutputStream::BgzfBlock::runJob() {
    if(!compress()) {
        cerr << "Bgzf block compression failed. Aborting." << endl;
        exit(-1);
    }
}

bool BgzfOutputStream::BgzfBlock::compress() {
    
    int current_compression_level = stream->compression_level;
    while(true) {
        // compress
        z_stream zs = {0};
        zs.zalloc    = NULL;
        zs.zfree     = NULL;
        zs.next_in   = (Bytef*)&uncompressed_data[0];
        zs.avail_in  = uncompressed_size;
        zs.next_out  = (Bytef*)&compressed_data[18];
        zs.avail_out = BGZF_BLOCK_SIZE - 18 - 8;
        
        // initialize the zlib compression algorithm
        int init_status = deflateInit2(&zs,
                                  current_compression_level,
                                  Z_DEFLATED,
                                  -15,  //window
                                  8,    //memory level
                                  Z_DEFAULT_STRATEGY);
        
        if ( init_status != Z_OK ) {
            cerr << "BGZF writer: zlib deflateInit2 failed" << endl;
            exit(-1);
        }
        // compress the data
        int deflate_status = deflate(&zs, Z_FINISH);
        
        
        // finalize the compression routine
        int end_status = deflateEnd(&zs);
        if ( end_status != Z_OK ){
            cerr << "BGZF writer: zlib deflateEnd failed. Aborting." << endl;
            exit(-1);
        }
        
        if(deflate_status == Z_STREAM_END) {
            compressed_size = zs.total_out + 18 + 8;

            break;
        } else {
            
            // there was not enough space available in buffer
            // try to reduce the input length & re-start loop
            if ( deflate_status == Z_OK ) {
                current_compression_level++;
                if ( current_compression_level > Z_BEST_COMPRESSION ){
                    cerr << "BGZF writer: input reduction failed" << endl;
                    exit(-1);
                }
                
                continue;
            }
            cerr << "BGZF writer: zlib deflate failed. Aborting." << endl;
            exit(-1);
        }
    }
    

    // fill in GZ data fields in buffer
    // some day we should avoid doing this manually and use deflateSetHeader()
    compressed_data[0] = 31;  //magic
    compressed_data[1] = 139; //magic
    compressed_data[2] = 8;   //deflate
    compressed_data[3] = 4;   //flags
    *((uint32_t *)&compressed_data[4]) = 0;    //mod time
    compressed_data[8] = 0;   //xflags
    compressed_data[9] = 255; //unknown OS
    compressed_data[10] = 6;  //xtra fields length
    compressed_data[11] = 0;  //xtra fields length
    compressed_data[12] = 66;// Subfield id 1
    compressed_data[13] = 67;// Subfield id 2
    compressed_data[14] = 2;// Subfield length (lower)
    compressed_data[15] = 0;// Subfield length (upper)
    *((uint16_t *)&compressed_data[16]) = compressed_size-1;    //block size minus one
    
    //data is in compressed_data[16..BGZF_BLOCK_SIZE-8]
    
    //now GZ footer..
    unsigned int data_end = compressed_size - 8;
    *((uint32_t *)&compressed_data[data_end]) = crc32(crc32(0, NULL, 0), (Bytef *)&uncompressed_data[0], uncompressed_size);
    *((uint32_t *)&compressed_data[data_end+4]) = uncompressed_size;
    
    return true;
}

bool BgzfOutputStream::BgzfBlock::isFull() {
    unsigned int full_size = stream->compression_level == 0 ? BGZF_BLOCK_SIZE - 64 : BGZF_BLOCK_SIZE;
    return full_size == uncompressed_size;
}

unsigned int BgzfOutputStream::BgzfBlock::addData(const char * data, unsigned int length) {
    unsigned int full_size = stream->compression_level == 0 ? BGZF_BLOCK_SIZE - 64 : BGZF_BLOCK_SIZE;
    
    unsigned int copy_size = min(length, full_size - uncompressed_size);
    
    memcpy(&uncompressed_data[uncompressed_size], data, copy_size);
    
    uncompressed_size += copy_size;
    
    return copy_size;
}

bool BgzfOutputStream::BgzfBlock::write(std::ofstream & out) {
    out.write(compressed_data, compressed_size);
    return !out.fail();
}

/////////////////////////
// BgzfOutputStream class

bool BgzfOutputStream::open(std::string filename) {
    output_stream.open(filename.c_str());
    
    if(output_stream.fail())
        return false;

    current_block = new BgzfBlock(this);
    
    return true;
}

void BgzfOutputStream::write(const char * data, size_t len) {
    while(len) {
        int written = current_block->addData(data, len);
        len -= written;
        data += written;
        
        if(current_block->isFull()) {
            current_block->compress();
            current_block->write(output_stream);
            delete current_block;
            current_block = new BgzfBlock(this);
        }
    }
}

void BgzfOutputStream::close() {
    current_block->compress();
    current_block->write(output_stream);
    delete current_block;
    
    //write empty block
    BgzfBlock empty(this);
    empty.compress();
    empty.write(output_stream);
    
    output_stream.close();
}
