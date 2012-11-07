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
#include "thread_pool.h"

class BgzfOutputStream {
    class BgzfCompressJob : public ThreadJob {
    public:
        Spinlock data_access_lock;  //only needed to get rid of race detection warnings in ThreadSanitizer
        std::vector<char> compressed_data;
        std::vector<char> uncompressed_data;
        bool compressed;
        int compression_level;
        BgzfOutputStream * stream;
        BgzfCompressJob()
        : compressed(false)
        {}
        void runJob();
    };
public:
    BgzfOutputStream()
    : compression_level(6)
    {}

    bool open(std::string filename);
    void write(const char * data, size_t len);
    void close();
    bool is_open() const { return output_stream.is_open(); }
    bool fail() { return output_stream.fail(); }
    void setCompressionLevel(int level) { compression_level = level; }
    void flushQueue();
protected:
    std::ofstream output_stream;
    int compression_level;
    std::vector<char> write_buffer;
    SynchronizedBlockingQueue<BgzfCompressJob *> job_queue;
    pthread_mutex_t write_mutex;
};

#endif