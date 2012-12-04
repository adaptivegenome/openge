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

const int BGZF_BLOCK_SIZE = 65536;

class BgzfInputStream
{
    class BgzfBlock : public ThreadJob {
        char compressed_data[BGZF_BLOCK_SIZE];
        char uncompressed_data[BGZF_BLOCK_SIZE];
        unsigned int compressed_size;
        unsigned int uncompressed_size;
        unsigned int read_size;
        
        SynchronizedFlag decompressed, decompression_started;
        Spinlock decompression_start;
        BgzfInputStream * stream;
    public:
        BgzfBlock(BgzfInputStream * stream)
        : read_size(0)
        , stream(stream)
        {
            decompressed.clear();
            decompression_started.clear();
        }
        unsigned int read();
        bool decompress();
        bool isDecompressed() { return decompressed.isSet(); }
        unsigned int readData(void * dest, unsigned int max_size);
        virtual void runJob();
        bool dataRemaining() { return read_size != uncompressed_size; }
    };
public:
    BgzfInputStream()
    {
        eof_seen.clear();
        fail_seen.clear();
    }
    bool open(std::string filename);
    bool read(char * data, size_t len);
    void close();
    bool is_open() { return *input_stream == std::cin || input_stream_real.is_open(); }
    bool eof() { return block_queue.empty() && eof_seen.isSet(); }
    bool fail() { return fail_seen.isSet(); }   //all errors are treated as fatal
protected:
    std::istream * input_stream;
    std::ifstream input_stream_real;
    SynchronizedFlag eof_seen, fail_seen;
    SynchronizedQueue<BgzfBlock *> block_queue;
    
    //multithreading:
    mutex read_signal_lock;
    condition_variable read_signal_cv;
    pthread_t read_thread;
    bool use_threads;
    
    static void * block_readproc(void * stream);
};

#endif