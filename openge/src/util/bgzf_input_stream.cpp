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

#include "bgzf_input_stream.h"
#include <zlib.h>

#include <iostream>
#include <cassert>
#include <cstring>
#include <stdint.h>

using namespace std;

///////////////////////////
// BgzfBlock implementation

void BgzfInputStream::BgzfBlock::runJob() {
    decompress();
}

unsigned int BgzfInputStream::BgzfBlock::read() {
    stream->input_stream->read(compressed_data, 18);
    
    if(stream->input_stream->eof()) {
        stream->eof_seen.set();
        return 0;
    }
    
    if(stream->input_stream->fail() && !stream->input_stream->eof()) {
        cerr << "Error reading BGZF block header. Aborting." << endl;
        exit(-1);
    }
    uint16_t length = *((uint16_t *)&compressed_data[16]) + 1;

    stream->input_stream->read(&compressed_data[18], length - 18);
    
    if(stream->input_stream->eof())
        stream->eof_seen.set();
    
    if(stream->input_stream->fail() && !stream->input_stream->eof()) {
        cerr << "Error reading BGZF block. Aborting." << endl;
        exit(-1);
    }
    
    compressed_size = stream->input_stream->gcount();
    
    assert(compressed_data[0] == 31 && compressed_data[1] == (char)139);
    
    return compressed_size;
}

bool BgzfInputStream::BgzfBlock::decompress() {
    
    decompression_start.lock();
    if(decompression_started.isSet()) {
        decompression_start.unlock();
        return true;
    }
    decompression_started.set();
    decompression_start.unlock();
    
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
    uncompressed_size = *((uint32_t *) &compressed_data[uncompressed_position]);
    
    {
        // setup zlib stream object
        z_stream zs;
        zs.zalloc    = NULL;
        zs.zfree     = NULL;
        zs.next_in   = (unsigned char *)data;
        zs.avail_in  = bsize - 16;
        zs.next_out  =  (unsigned char *) &uncompressed_data[0];
        zs.avail_out = 65536;
        
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
    }
    
    return true;
}

unsigned int BgzfInputStream::BgzfBlock::readData(void * dest, unsigned int max_size) {
    assert(isDone());
    
    unsigned int actual_read_len = min(max_size, uncompressed_size - read_size);
    
    memcpy(dest, &uncompressed_data[read_size], actual_read_len);
    read_size += actual_read_len;
    
    return actual_read_len;
}

/////////////////////////////////
// BgzfInputStream implementation

bool BgzfInputStream::open(string filename) {
    if(filename == "stdin") {
        input_stream = &cin;
    } else {
        input_stream = &input_stream_real;
        input_stream_real.open(filename.c_str());
        
        if(input_stream->fail()) {
            cerr << "BGZF open() failed for " << filename << endl;
            return false;
        }
    }
    
    int ret = pthread_create(&read_thread, NULL, block_readproc, this);
    if(0 != ret) {
        cerr << "Error creating BGZF read thread. Quitting. error = " << ret << endl;
        exit(-1);
    }
    
    return true;
}

void * BgzfInputStream::block_readproc(void * data) {
    BgzfInputStream * stream = (BgzfInputStream *) data;
    
    while(!stream->eof_seen.isSet()) {
        stream->read_signal_lock.lock();
        while(stream->block_queue.size() >= 100) {
            if(stream->eof_seen.isSet())
                return NULL;
            stream->read_signal_cv.wait(stream->read_signal_lock);
        }
        stream->read_signal_lock.unlock();
        
        while(stream->block_queue.size() < 100 && !stream->eof_seen.isSet()) {
            BgzfBlock * block = new BgzfBlock(stream);
            int read = block->read();
            if(read) {
                stream->block_queue.push(block);
                if(OGEParallelismSettings::isMultithreadingEnabled())
                    ThreadPool::sharedPool()->addJob(block);
                else
                    block->runJob();
            }
        }
    }
    
    return NULL;
}

bool BgzfInputStream::read(char * data, size_t len) {
    unsigned int read_len = 0;
    while(read_len != len) {
        while(block_queue.empty() || !block_queue.front()->isDecompressed()) {
            if(block_queue.empty() && eof_seen.isSet()) {
                return false;
            }
            if (!block_queue.empty() && !block_queue.front()->isDecompressed())
                block_queue.front()->decompress();
            usleep(50000);
        }
        
        BgzfBlock * block = block_queue.front();
        assert(block->isDecompressed());
        
        unsigned int actual_read_length = block->readData(&((char *)data)[read_len], len - read_len);
        
        read_len += actual_read_length;
        
        if(!block->dataRemaining()) {
            delete block;
            block_queue.pop();
            if(!eof_seen.isSet()) {
                //request another block
                read_signal_lock.lock();
                read_signal_cv.notify_one();
                read_signal_lock.unlock();
            }
        }
    }
    
    return true;
}
void BgzfInputStream::close() {
    
    eof_seen.set();
    
    read_signal_lock.lock();
    read_signal_cv.notify_one();
    read_signal_lock.unlock();
    
    int ret = pthread_join(read_thread, NULL);
    if(0 != ret) {
        cerr << "Error joining BGZF read thread (error " << ret << ")." << endl;
    }
    
    if(input_stream_real.is_open())
        input_stream_real.close();
}
