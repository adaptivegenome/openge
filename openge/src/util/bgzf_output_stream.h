#ifndef OGE_BGZF_OUTPUT_STREAM_H
#define OGE_BGZF_OUTPUT_STREAM_H

/*********************************************************************
 *
 * bgzf_output_stream.h: Use BGZF compression on a stream.
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
#include <vector>
#include <stdint.h>
#include "thread_pool.h"

const uint32_t BGZF_BLOCK_SIZE = 65536;

class BgzfOutputStream {
    int compression_level;
    std::ofstream output_stream;
    bool use_threads;

    class BgzfBlock : public ThreadJob {
        BgzfOutputStream * stream;
        char uncompressed_data[BGZF_BLOCK_SIZE];
        char compressed_data[BGZF_BLOCK_SIZE];
        unsigned int uncompressed_size, compressed_size;
        Spinlock data_access_lock;
        SychronizedFlag compress_finished;
    public:
        BgzfBlock(BgzfOutputStream * stream)
        : stream(stream)
        , uncompressed_size(0)
        {
            compress_finished.clear();
        }
        
        unsigned int addData(const char * data, unsigned int length);
        bool isFull();
        bool isCompressed() { return compress_finished.isSet(); }
        void runJob();  //calls compress when run in thread pool
        bool compress();
        bool write(std::ofstream & out);
    };

    BgzfBlock * current_block;
    
    //multithreading:
    pthread_t write_thread;
    condition_variable write_thread_signal;
    mutex write_thread_mutex;
    SychronizedFlag closing;
    SynchronizedQueue<BgzfBlock *> write_queue;
    
    static void * write_threadproc(void * stream_p);
public:
    BgzfOutputStream()
    : compression_level(6)
    , use_threads(true)
    {
        closing.clear();
    }
    bool open(std::string filename);
    void write(const char * data, size_t len);
    void close();
    bool is_open() const { return output_stream.is_open(); }
    bool fail() { return output_stream.fail(); }
    void setCompressionLevel(int level) { compression_level = level; }
};

#endif