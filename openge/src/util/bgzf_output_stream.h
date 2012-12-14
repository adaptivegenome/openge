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
    std::ostream * output_stream;
    std::ofstream output_stream_real;
    bool use_threads;

    class BgzfBlock : public ThreadJob {
        BgzfOutputStream * stream;
        char uncompressed_data[BGZF_BLOCK_SIZE];
        char compressed_data[BGZF_BLOCK_SIZE];
        unsigned int uncompressed_size, compressed_size;
        Spinlock data_access_lock;
    public:
        BgzfBlock(BgzfOutputStream * stream)
        : stream(stream)
        , uncompressed_size(0)
        { }
        
        unsigned int addData(const char * data, unsigned int length);
        bool isFull();
        bool isCompressed() { return isDone(); }
        void runJob();  //calls compress when run in thread pool
        bool compress();
        bool write();
    };

    BgzfBlock * current_block;
    
    //multithreading:
    pthread_t write_thread;
    condition_variable write_thread_signal;
    mutex write_thread_mutex;
    SynchronizedFlag closing;
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
    bool is_open() const { if(output_stream == &output_stream_real) return output_stream_real.is_open(); else return true; }
    bool fail() { return output_stream->fail(); }
    void setCompressionLevel(int level) { compression_level = level; }
};

#endif