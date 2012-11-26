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
    
    data_access_lock.lock();
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
    
    if(stream->use_thread_pool) {
        data_access_lock.unlock();
        compressed_flag_lock.lock();
        compressed = true;
        compressed_flag_lock.unlock();
        
        stream->flushQueue();
    } else {
        stream->output_stream.write(&compressed_data[0], compressed_data.size());
        data_access_lock.unlock();
    }
}

void BgzfOutputStream::flushBlocks() {
    while( !job_queue.empty()) {
        BgzfCompressJob * job = job_queue.front();
        
        job->compressed_flag_lock.lock();
        if(!job_queue.front()->compressed) {
            job->compressed_flag_lock.unlock();
            break;
        }
        job->compressed_flag_lock.unlock();
        
        job->data_access_lock.lock();
        output_stream.write(&job->compressed_data[0], job->compressed_data.size());
        job->data_access_lock.unlock();
        
        job_queue.pop();
        delete job;
    }
}

void BgzfOutputStream::writeEof() {
    //BGZF uses an empty ZIP block to indicate the end of file.
    
    //write an empty block ("EOF marker")
	static const uint8_t empty_block[29] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0";
    output_stream.write((const char *)empty_block, 28);
    
    output_stream.close();
}

void BgzfOutputStream::flushQueue() {
    //if multithreaded, send a signal to the flush queue
    if(use_thread_pool) {
        int error = pthread_mutex_lock(&write_wait_mutex);
        if(0 != error)
        {
            cerr << "Error locking BGZF flush mutex for signalling. Aborting. (" << error << ")" << endl;
            exit(-1);
        }

        error = pthread_cond_signal(&flush_signal);
        if(0 != error)
        {
            cerr << "Error waiting for BGZF flush CV. Aborting. (" << error << ")" << endl;
            exit(-1);
        }

        error = pthread_mutex_unlock(&write_wait_mutex);
        if(0 != error)
        {
            cerr << "Error unlocking BGZF flush mutex for signalling. Aborting. (" << error << ")" << endl;
            exit(-1);
        }
    //if not multithreaded, then do the flush here
    } else {
        flushBlocks();
    }
}

void * BgzfOutputStream::file_write_threadproc(void * data) {
    BgzfOutputStream * stream = (BgzfOutputStream * ) data;
    int error = pthread_mutex_lock(&stream->write_mutex);
    if(0 != error)
    {
        cerr << "Error locking BGZF write mutex. Aborting. (" << error << ")" << endl;
        exit(-1);
    }
    
    //wait for flush events
    while(!stream->closing) {
        error = pthread_mutex_lock(&stream->write_wait_mutex);
        if(0 != error)
        {
            cerr << "Error locking BGZF flush mutex. Aborting. (" << error << ")" << endl;
            exit(-1);
        }
        
        while(stream->job_queue.empty() && !stream->closing) {
            int error = pthread_cond_wait(&stream->flush_signal, &stream->write_wait_mutex);
            if(0 != error)
            {
                cerr << "Error waiting for BGZF flush CV. Aborting. (" << error << ")" << endl;
                exit(-1);
            }
        }

        error = pthread_mutex_unlock(&stream->write_wait_mutex);
        if(0 != error)
        {
            cerr << "Error unlocking BGZF flush mutex. Aborting. (" << error << ")" << endl;
            exit(-1);
        }

        stream->flushBlocks();
    }
    
    stream->writeEof();
    
    error = pthread_mutex_unlock(&stream->write_mutex);
    if(0 != error)
    {
        cerr << "Error unlocking BGZF write mutex. Aborting. (" << error << ")" << endl;
        exit(-1);
    }
    
    return NULL;
}

bool BgzfOutputStream::open(std::string filename) {
    use_thread_pool = OGEParallelismSettings::isMultithreadingEnabled();
    output_stream.open(filename.c_str());
    
    if(output_stream.fail())
        return false;
    
    if(use_thread_pool) {
    int ret = pthread_mutex_init(&write_mutex, NULL);
        if(0 != ret) {
            cerr << "Error opening BGZF write mutex. Aborting. (error " << ret << ")." << endl;
            exit(-1);
        }

        ret = pthread_mutex_init(&write_wait_mutex, NULL);
        if(0 != ret) {
            cerr << "Error opening BGZF write wait mutex. Aborting. (error " << ret << ")." << endl;
            exit(-1);
        }

        ret = pthread_cond_init(&flush_signal, NULL);
        if(0 != ret) {
            cerr << "Error opening BGZF flush signal. Aborting. (error " << ret << ")." << endl;
            exit(-1);
        }
        
        ret = pthread_create(&write_thread, NULL, file_write_threadproc, this);
        if(0 != ret) {
            cerr << "Error creating BGZF write thread. (error " << ret << ")." << endl;
            exit(-1);
        }
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
        job->data_access_lock.lock();
        job->compression_level = compression_level;
        job->stream = this;
        job->uncompressed_data.insert(job->uncompressed_data.begin(), write_buffer.begin(), write_buffer.begin() + block_write_size);
        write_buffer.erase(write_buffer.begin(), write_buffer.begin() + block_write_size);
        job->data_access_lock.unlock();

        job_queue.push(job);
        if(use_thread_pool)
            ThreadPool::sharedPool()->addJob(job);
        else {
            job->runJob();
            delete job;
        }
    }
}

void BgzfOutputStream::close() {
    //write last bit of data out
    if(!write_buffer.empty()) {
        BgzfCompressJob * job = new BgzfCompressJob();
        job->data_access_lock.lock();
        job->compression_level = compression_level;
        job->stream = this;
        job->uncompressed_data.insert(job->uncompressed_data.begin(), write_buffer.begin(), write_buffer.end());
        write_buffer.clear();
        job->data_access_lock.unlock();

        job_queue.push(job);
        if(use_thread_pool)
            ThreadPool::sharedPool()->addJob(job);
        else {
            job->runJob();
            delete job;
        }
    }
    
    //wait for all compression to finish
    if(use_thread_pool) {
        ThreadPool::sharedPool()->waitForJobCompletion();
        
        flushQueue();
        
        //and we are done, wait for write thread to finish.
        closing = true; //signal close to write thread

        //cause another flush so if the thread is waiting on signal, it will see that we changed the closing flag
        flushQueue();
        
        int ret = pthread_join(write_thread, NULL);
        
        if(0 != ret)
            cerr << "Error joining BGZF write thread (error " << ret << "." << endl;
    }
    
    flushBlocks();
    writeEof();
    
    if(use_thread_pool) {
        int ret = pthread_mutex_destroy(&write_mutex);
        if(0 != ret)
            cerr << "Error closing BGZF write mutex (error " << ret << "." << endl;
        
        ret = pthread_mutex_destroy(&write_wait_mutex);
        if(0 != ret)
            cerr << "Error closing BGZF flush wait mutex (error " << ret << "." << endl;
        
        ret = pthread_cond_destroy(&flush_signal);
        if(0 != ret)
            cerr << "Error closing BGZF flush signal (error " << ret << "." << endl;
    }
}