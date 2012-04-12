// ***************************************************************************
// BgzfStream_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 22 March 2012 (Lee Baker) - block cache
// ---------------------------------------------------------------------------
// Based on BGZF routines developed at the Broad Institute.
// Provides the basic functionality for reading & writing BGZF files
// Replaces the old BGZF.* files to avoid clashing with other toolkits
// ***************************************************************************

#ifndef BGZFSTREAM_P_H
#define BGZFSTREAM_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/api_global.h"
#include "api/BamAux.h"
#include "api/IBamIODevice.h"
#include "api/internal/utils/BamThreadPool.h"

#include <semaphore.h>

#include <string>

namespace BamTools {
namespace Internal {

class BgzfStream {
	
	class BgzfBlockCompressJob : public BamThreadJob
	{
	public:
		BgzfBlockCompressJob(BgzfStream * stream, pthread_mutex_t * writing_start_lock, pthread_mutex_t * writing_stop_lock);
    virtual ~BgzfBlockCompressJob();
		void runJob();
		size_t DeflateBlock(int32_t blockLength);
			
		protected:
			RaiiBuffer m_uncompressedBlock;
			RaiiBuffer * m_compressedBlock;
			BgzfStream * stream;
            pthread_mutex_t * writing_start_lock;
            pthread_mutex_t * writing_stop_lock;
			int64_t m_blockOffset;
    public:
      char writing_start_lock_name[32]; //we store the name so we can sem_unlink it when done.
		};
  class BgzfBlockCache
  {
  public:
    BgzfBlockCache(BgzfStream * stream);
    ~BgzfBlockCache();
    
    bool LoadBlock();
    void ReadBlock(void);
    int64_t Tell();
    size_t InflateBlock(const size_t& blockLength);
  protected:
    sem_t * wait_for_read;
    sem_t * wait_for_copy;
    
    pthread_mutex_t block_cache_modify;
    
    pthread_t worker_thread;
    bool thread_exit;
    static void * thread_start(void * block_cache);
    
    struct BgzfBlockCacheBlock
    {
    public:
      RaiiBuffer buffer;
      int32_t m_blockOffset, starting_blockOffset;
      int64_t m_blockAddress, starting_blockAddress;
      int32_t m_blockLength;
      size_t m_tell;
      BgzfBlockCacheBlock() : buffer(Constants::BGZF_MAX_BLOCK_SIZE) {}
    };
    
    // Keep our own copy of these variables so we can use the ReadBlock method unmodified.
    // These variables are identical to BgzfStream internal variables.
    int32_t m_blockLength;
    int32_t m_blockOffset;
    int64_t m_blockAddress;
    RaiiBuffer m_uncompressedBlock;
    RaiiBuffer m_compressedBlock;
    
    bool cache_error; //< When an error is encountered reading from the block cache, signal a fall back to normal reads
    int64_t m_tell; //< Cache a tell response
    
    IBamIODevice* m_device;
    BgzfStream * m_stream;
    SynchronizedQueue<BgzfBlockCacheBlock *> block_queue;
    sem_t block_available;
  };

    // constructor & destructor
    public:
        BgzfStream(void);
        ~BgzfStream(void);

    // main interface methods
    public:
        // closes BGZF file
        void Close(void);
        // returns true if BgzfStream open for IO
        bool IsOpen(void) const;
        // opens the BGZF file
        void Open(const std::string& filename, const IBamIODevice::OpenMode mode);
        // reads BGZF data into a byte buffer
        size_t Read(char* data, const size_t dataLength);
        // seek to position in BGZF file
        void Seek(const int64_t& position);
        // sets IO device (closes previous, if any, but does not attempt to open)
        void SetIODevice(IBamIODevice* device);
        // enable/disable compressed output
        void SetWriteCompressed(bool ok);
        // set how compressed the output should be
        void SetCompressionLevel(int compressionLevel);
        // get file position in BGZF file
        int64_t Tell(void) const;
        // writes the supplied data into the BGZF buffer
        size_t Write(const char* data, const size_t dataLength);

    // internal methods
    private:
        // compresses the current block
        size_t DeflateBlock(int32_t blockLength);
        // flushes the data in the BGZF block
        void FlushBlock(void);
        // de-compresses the current block
        size_t InflateBlock(const size_t& blockLength);
        // reads a BGZF block
        void ReadBlock(void);

    // static 'utility' methods
    public:
        // checks BGZF block header
        static bool CheckBlockHeader(char* header);
        
        // To use one thread pool for all compression in the app, 
        // use the following functions before/after all BamWriter open/closes.
        static void OpenSharedThreadPool();
        static void CloseSharedThreadPool();

    // data members
    public:
        int32_t m_blockLength;
        int32_t m_blockOffset;
        int64_t m_blockAddress;
        int32_t m_writeSemaphoreCount;
        int m_compressionLevel;

        bool m_isWriteCompressed;
        IBamIODevice* m_device;

        RaiiBuffer m_uncompressedBlock;
        RaiiBuffer m_compressedBlock;
        BamThreadPool * m_thread_pool;
        pthread_mutex_t * m_lastWriteLock;
        BgzfBlockCache * m_cache;
};

} // namespace Internal
} // namespace BamTools

#endif // BGZFSTREAM_P_H
