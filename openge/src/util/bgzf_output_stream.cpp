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
 *********************************************************************/

#include "bgzf_output_stream.h"
#include "api/api_global.h"
#include "thread_pool.h"
#include <zlib.h>
#include <iostream>
#include <stdint.h>
#include <cstring>

using namespace std;

const uint32_t BGZF_BLOCK_SIZE       = 65536;

const int8_t   GZIP_WINDOW_BITS          = -15;
const int8_t   Z_DEFAULT_MEM_LEVEL       = 8;
const uint8_t  BGZF_BLOCK_HEADER_LENGTH  = 18;
const uint8_t  BGZF_BLOCK_FOOTER_LENGTH  = 8;

void BgzfOutputStream::BgzfCompressJob::runJob() {
    
    ogeNameThread("ogeBgzfDeflate");
    compressed_data.resize(BGZF_BLOCK_SIZE);

    // initialize the gzip header
    char* buffer = &compressed_data[0];
    memset(buffer, 0, 18);
    buffer[0]  = 31;    //magic
    buffer[1]  = 139;   //magic
    buffer[2]  = 8;     //deflate method
    buffer[3]  = 4;     //extra fields
    buffer[9]  = 255;   //unknown OS
    buffer[10] = 6;     //XLen
    buffer[12] = 'B';   //xtra field id
    buffer[13] = 'C';   //xtra field id
    buffer[14] = 2;     //xtra field len
    
    // loop to retry for blocks that do not compress enough
    size_t compressedLength = 0;
    
    while ( true ) {
        
        // initialize zstream values
        z_stream zs = {0};
        zs.zalloc    = NULL;
        zs.zfree     = NULL;
        zs.next_in   = (Bytef*)&uncompressed_data[0];
        zs.avail_in  = uncompressed_data.size();
        zs.next_out  = (Bytef*)&compressed_data[BGZF_BLOCK_HEADER_LENGTH];
        zs.avail_out = compressed_data.size() -
        BGZF_BLOCK_HEADER_LENGTH -
        BGZF_BLOCK_FOOTER_LENGTH;
        
        // initialize the zlib compression algorithm
        int status = deflateInit2(&zs,
                                  compression_level,
                                  Z_DEFLATED,
                                  GZIP_WINDOW_BITS,
                                  Z_DEFAULT_MEM_LEVEL,
                                  Z_DEFAULT_STRATEGY);

        if ( status != Z_OK ) {
            cerr << "BGZF writer: zlib deflateInit2 failed" << endl;
            exit(-1);
        }
        
        // compress the data
        status = deflate(&zs, Z_FINISH);
        
        // if not at stream end
        if ( status != Z_STREAM_END ) {
            
            deflateEnd(&zs);
            
            // there was not enough space available in buffer
            // try to reduce the input length & re-start loop
            if ( status == Z_OK ) {
                compression_level++;
                if ( compression_level > Z_BEST_COMPRESSION ){
                    cerr << "BGZF writer: input reduction failed" << endl;
                    exit(-1);
                }

                continue;
            }
            cerr << "BGZF writer: zlib deflate failed. Aborting." << endl;
            exit(-1);
        }
        
        // finalize the compression routine
        status = deflateEnd(&zs);
        if ( status != Z_OK ){
            cerr << "BGZF writer: zlib deflateEnd failed. Aborting." << endl;
            exit(-1);
        }
        
        // update compressedLength
        compressedLength = zs.total_out +
        BGZF_BLOCK_HEADER_LENGTH +
        BGZF_BLOCK_FOOTER_LENGTH;

        compressed_data.resize(compressedLength);
        if ( compressedLength > BGZF_BLOCK_SIZE ) {
            if(compression_level > Z_BEST_COMPRESSION) {
                cerr << "BGZF writer: compression overflow. Aborting." << endl;
                exit(-1);
            }
            else {
                compression_level++;
                continue;
            }
        }

        // quit while loop
        break;
    }

    // store the compressed length
    *((uint16_t *)&buffer[16]) = (uint16_t)(compressedLength - 1);

    // store the CRC32 checksum
    uint32_t crc = crc32(0, NULL, 0);
    crc = crc32(crc, (Bytef*)&uncompressed_data[0], uncompressed_data.size());
    *((uint32_t *) &buffer[compressedLength - 8]) = crc;
    *((uint32_t *) &buffer[compressedLength - 4]) = uncompressed_data.size();

    compressed = true;
    
    stream->flushQueue();
}

void BgzfOutputStream::flushQueue() {
    int error = pthread_mutex_lock(&write_mutex);
    if(0 != error)
    {
        cerr << "Error locking BGZF write mutex. Aborting. (" << error << ")" << endl;
        exit(-1);
    }
    
    while( !job_queue.empty() && job_queue.front()->compressed ) {
        BgzfCompressJob * job = job_queue.front();

        output_stream.write(&job->compressed_data[0], job->compressed_data.size());
        
        job_queue.pop();
        delete job;
    }
    
    error = pthread_mutex_unlock(&write_mutex);
    if(0 != error)
    {
        cerr << "Error unlocking BGZF write mutex. Aborting. (" << error << ")" << endl;
        exit(-1);
    }
}

bool BgzfOutputStream::open(std::string filename) {
    output_stream.open(filename.c_str());
    
    if(output_stream.fail())
        return false;
    
    int ret = pthread_mutex_init(&write_mutex, NULL);
    if(0 != ret) {
        cerr << "Error opening BGZF write mutex. Aborting. (error " << ret << ")." << endl;
        exit(-1);
    }
    
    return true;
}

void BgzfOutputStream::write(const char * data, size_t len) {
    write_buffer.insert(write_buffer.end(), data, data + len);
    
    size_t block_write_size = BGZF_BLOCK_SIZE;
    
    // if we aren't compressing, we need to 
    if(compression_level == Z_NO_COMPRESSION)
        block_write_size -= 1024;
    
    while(write_buffer.size() > block_write_size) {
        BgzfCompressJob * job = new BgzfCompressJob();
        job->compression_level = compression_level;
        job->stream = this;
        job->uncompressed_data.insert(job->uncompressed_data.begin(), write_buffer.begin(), write_buffer.begin() + block_write_size);
        write_buffer.erase(write_buffer.begin(), write_buffer.begin() + block_write_size);

        job_queue.push(job);
        ThreadPool::sharedPool()->addJob(job);
    }
}

void BgzfOutputStream::close() {
    //write last bit of data out
    if(!write_buffer.empty()) {
        BgzfCompressJob * job = new BgzfCompressJob();
        job->compression_level = compression_level;
        job->stream = this;
        job->uncompressed_data.insert(job->uncompressed_data.begin(), write_buffer.begin(), write_buffer.end());
        write_buffer.clear();

        job_queue.push(job);
        ThreadPool::sharedPool()->addJob(job);
    }
    
    //wait for last write to be flushed before closing read.
    while(!job_queue.empty())
        usleep(20000);
    
    int ret = pthread_mutex_lock(&write_mutex);
    if(0 != ret)
    {
        cerr << "Error locking BGZF write mutex at close. Aborting. (error " << ret << ")." << endl;
        exit(-1);
    }

    //write an empty block ("EOF marker")
	static const uint8_t empty_block[29] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0";
    output_stream.write((const char *)empty_block, 28);

    output_stream.close();

    ret = pthread_mutex_unlock(&write_mutex);
    if(0 != ret)
    {
        cerr << "Error unlocking BGZF write mutex at close. Aborting. (error " << ret << ")." << endl;
        exit(-1);
    }

    ret = pthread_mutex_destroy(&write_mutex);
    if(0 != ret)
        cerr << "Error closing BGZF write mutex (error " << ret << "." << endl;
}