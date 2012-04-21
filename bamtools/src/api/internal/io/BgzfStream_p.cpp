// ***************************************************************************
// BgzfStream_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 17 January 2012(DB)
// ---------------------------------------------------------------------------
// Based on BGZF routines developed at the Broad Institute.
// Provides the basic functionality for reading & writing BGZF files
// Replaces the old BGZF.* files to avoid clashing with other toolkits
// ***************************************************************************

#include "api/BamAux.h"
#include "api/BamConstants.h"
#include "api/internal/io/BamDeviceFactory_p.h"
#include "api/internal/io/BgzfStream_p.h"
#include "api/internal/utils/BamException_p.h"
#include "api/BamParallelismSettings.h"
using namespace BamTools;
using namespace BamTools::Internal;

#ifdef __linux__
#include <sys/prctl.h>
#endif

#include "zlib.h"

#include <cstring>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <errno.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <semaphore.h>
using namespace std;

#define BLOCKCACHE_NUM_BLOCKS 8

#define BLOCKCACHE_WRITE_SEMAPHORE_NAME "write_sem_%x_%d"

// We can run deflate for separate blocks in parallel. This function
// compresses and writes a single 64k block. A pair of semaphores protects
// this job from writing its results out of order- there is essentially a
// chain of semaphores protecting the actual file so that writes are committed
// on at a time, and in sequence.

// ---------------------------
// BgzfBlockCompressJob implementation
// ---------------------------

BgzfStream::BgzfBlockCompressJob::BgzfBlockCompressJob(BgzfStream * stream, pthread_mutex_t *writing_start_lock, pthread_mutex_t * writing_stop_lock)
: m_uncompressedBlock(Constants::BGZF_DEFAULT_BLOCK_SIZE)
, m_compressedBlock(NULL)
, stream(stream)
, writing_start_lock(writing_start_lock)
, writing_stop_lock(writing_stop_lock)
, m_blockOffset(stream->m_blockOffset)
{
	//copy the current buffer to be written into this job, so the main thread can move on
	memcpy(m_uncompressedBlock.Buffer, stream->m_uncompressedBlock.Buffer, Constants::BGZF_DEFAULT_BLOCK_SIZE);
}

BgzfStream::BgzfBlockCompressJob::~BgzfBlockCompressJob()
{
  if(m_compressedBlock)
    delete m_compressedBlock;
}

void BgzfStream::BgzfBlockCompressJob::runJob()
{
#ifdef __linux__
    prctl(PR_SET_NAME,"bt_bgzf_compress",0,0,0);
#endif

    vector<size_t> write_lengths;
    vector<RaiiBuffer *> write_data;
    vector<int32_t> write_blockOffset;

    while ( m_blockOffset > 0 ) {
        // compress the data block
        RaiiBuffer * buffer = new RaiiBuffer(Constants::BGZF_DEFAULT_BLOCK_SIZE);
        m_compressedBlock = buffer;
        const size_t blockLength = DeflateBlock(m_blockOffset);

        write_lengths.push_back(blockLength);
        write_data.push_back(buffer);
        write_blockOffset.push_back(m_blockOffset);
    }
    
#ifdef __linux__
    prctl(PR_SET_NAME,"bt_bgzf_wt_write",0,0,0);
#endif
	//acquire semaphore, now we are the thread that should write next
    if(0 != pthread_mutex_lock(writing_start_lock))
        perror("Error waiting for writing_start_lock");
    if(0 != pthread_mutex_unlock(writing_start_lock))
        perror("Error unlocking for writing_start_lock before closing");
	if(0 != pthread_mutex_destroy(writing_start_lock))
        perror("Error closing writing_start_lock");
    
    delete writing_start_lock;
    
#ifdef __linux__
    prctl(PR_SET_NAME,"bt_bgzf_write",0,0,0);
#endif
    for(size_t i = 0; i < write_lengths.size(); i++)
    {
        m_compressedBlock = write_data[i];
        const size_t blockLength = write_lengths[i];
        m_blockOffset = write_blockOffset[i];

        // flush the data to our output device
        const int64_t numBytesWritten = stream->m_device->Write(m_compressedBlock->Buffer, blockLength);

        // check for device error
        if ( numBytesWritten < 0 ) {
            const string message = string("device error: ") + stream->m_device->GetErrorString();
            throw BamException("BgzfStream::BgzfBlockCompressJob::FlushBlock", message);
        }

        // check that we wrote expected numBytes
        if ( numBytesWritten != static_cast<int64_t>(blockLength) ) {
            stringstream s("");
            s << "expected to write " << blockLength
            << " bytes during flushing, but wrote " << numBytesWritten << " errno " << errno;
            perror("Write failed:");
            throw BamException("BgzfStream::BgzfBlockCompressJob::FlushBlock", s.str());
        }

        delete(m_compressedBlock);
        m_compressedBlock = NULL;
    }

	//allow the next thread to write, we're done
	if(0 != pthread_mutex_unlock(writing_stop_lock))
        perror("Error unlocking writing_stop_lock");
}
	
// compresses the current block
// this is almost identical to the version in BgzfStream, but 
// references member variables of BgzfBlockCompressJob instead
size_t BgzfStream::BgzfBlockCompressJob::DeflateBlock(int32_t blockLength)  {

  // initialize the gzip header
  char* buffer = m_compressedBlock->Buffer;
  memset(buffer, 0, 18);
  buffer[0]  = Constants::GZIP_ID1;
  buffer[1]  = Constants::GZIP_ID2;
  buffer[2]  = Constants::CM_DEFLATE;
  buffer[3]  = Constants::FLG_FEXTRA;
  buffer[9]  = Constants::OS_UNKNOWN;
  buffer[10] = Constants::BGZF_XLEN;
  buffer[12] = Constants::BGZF_ID1;
  buffer[13] = Constants::BGZF_ID2;
  buffer[14] = Constants::BGZF_LEN;

  // loop to retry for blocks that do not compress enough
  int inputLength = blockLength;
  size_t compressedLength = 0;
  const unsigned int bufferSize = Constants::BGZF_MAX_BLOCK_SIZE;

  // since a full block at no compression grows and has to go through the decompression
  // loop twice, we should just anticipate this and go ahead and avoid the failed pass
  //if(stream->m_compressionLevel == 0 && inputLength == Constants::BGZF_MAX_BLOCK_SIZE)
  //  inputLength -= 1024;

  while ( true ) {

    // initialize zstream values
    z_stream zs;
    zs.zalloc    = NULL;
    zs.zfree     = NULL;
    zs.next_in   = (Bytef*)m_uncompressedBlock.Buffer;
    zs.avail_in  = inputLength;
    zs.next_out  = (Bytef*)&buffer[Constants::BGZF_BLOCK_HEADER_LENGTH];
    zs.avail_out = bufferSize -
    Constants::BGZF_BLOCK_HEADER_LENGTH -
    Constants::BGZF_BLOCK_FOOTER_LENGTH;

    // initialize the zlib compression algorithm
    int status = deflateInit2(&zs,
                              stream->m_compressionLevel,
                              Z_DEFLATED,
                              Constants::GZIP_WINDOW_BITS,
                              Constants::Z_DEFAULT_MEM_LEVEL,
                              Z_DEFAULT_STRATEGY);
    if ( status != Z_OK )
      throw BamException("BgzfStream::DeflateBlock", "zlib deflateInit2 failed");

    // compress the data
    status = deflate(&zs, Z_FINISH);

    // if not at stream end
    if ( status != Z_STREAM_END ) {
      deflateEnd(&zs);

      // there was not enough space available in buffer
      // try to reduce the input length & re-start loop
      if ( status == Z_OK ) {
        inputLength -= 1024;
        if ( inputLength < 0 )
          throw BamException("BgzfStream::DeflateBlock", "input reduction failed");
        continue;
      }

      throw BamException("BgzfStream::DeflateBlock", "zlib deflate failed");
    }

    // finalize the compression routine
    status = deflateEnd(&zs);
    if ( status != Z_OK )
      throw BamException("BgzfStream::DeflateBlock", "zlib deflateEnd failed");

    // update compressedLength
    compressedLength = zs.total_out +
    Constants::BGZF_BLOCK_HEADER_LENGTH +
    Constants::BGZF_BLOCK_FOOTER_LENGTH;
    if ( compressedLength > Constants::BGZF_MAX_BLOCK_SIZE )
      throw BamException("BgzfStream::DeflateBlock", "deflate overflow");
    
    // quit while loop
    break;
  }
  
  // store the compressed length
  BamTools::PackUnsignedShort(&buffer[16], static_cast<uint16_t>(compressedLength - 1));
  
  // store the CRC32 checksum
  uint32_t crc = crc32(0, NULL, 0);
  crc = crc32(crc, (Bytef*)m_uncompressedBlock.Buffer, inputLength);
  BamTools::PackUnsignedInt(&buffer[compressedLength - 8], crc);
  BamTools::PackUnsignedInt(&buffer[compressedLength - 4], inputLength);
  
  // ensure that we have less than a block of data left
  int remaining = blockLength - inputLength;
  if ( remaining > 0 ) {
    if ( remaining > inputLength )
      throw BamException("BgzfStream::DeflateBlock", "after deflate, remainder too large");
    memcpy(m_uncompressedBlock.Buffer, m_uncompressedBlock.Buffer + inputLength, remaining);
  }
  
  // update block data
  m_blockOffset = remaining;

  // return result
  return compressedLength;
}

// ---------------------------
// BgzfBlockCache implementation
// ---------------------------

BgzfStream::BgzfBlockCache::BgzfBlockCache(BgzfStream * stream)
: thread_exit(false)
, m_blockLength(0)
, m_blockOffset(0)
, m_blockAddress(0)
, m_uncompressedBlock(Constants::BGZF_DEFAULT_BLOCK_SIZE)
, m_compressedBlock(Constants::BGZF_MAX_BLOCK_SIZE)
, cache_error(false)
, m_tell(0)
, m_stream(stream)
{
  char sem_name[64];
  sprintf(sem_name, "sem_bc_read_%08x", (int32_t) (int64_t) this);
  sem_unlink(sem_name);
  wait_for_read = sem_open(sem_name, O_CREAT | O_EXCL, 0, BLOCKCACHE_NUM_BLOCKS);
  if(wait_for_read == SEM_FAILED)
    perror("Error creating BgzfBlockCache semaphore");
  
  sprintf(sem_name, "sem_bc_q_%08x", (int32_t) (int64_t) this);
  sem_unlink(sem_name);
  wait_for_copy = sem_open(sem_name, O_CREAT | O_EXCL, 0, 0);
  if(wait_for_copy == SEM_FAILED)
    perror("Error creating BgzfBlockCache queue semaphore");

  //pthread_mutex_init(&block_cache_modify, NULL);
  
  int ret = pthread_create(&worker_thread, NULL, thread_start, this);
  if(ret != 0)
    perror("Error creating block cache thread");
}

BgzfStream::BgzfBlockCache::~BgzfBlockCache()
{
  thread_exit = true;
  if(0 != sem_post(wait_for_read))
    perror("Error posting wait_for_read");
  
  if(!cache_error)
    pthread_join(worker_thread, NULL);
  
  //flush queues so we don't leak memory
  while(0 < block_queue.size())
    delete block_queue.pop();
  
  if(0 != sem_close(wait_for_read))
    perror("Error closing wait_for_read");
  if(0 != sem_close(wait_for_copy))
    perror("Error closing wait_for_copy");
  
  char sem_name[64];
  sprintf(sem_name, "sem_bc_read_%08x", (int32_t) (int64_t) this);
  sem_unlink(sem_name);
  sprintf(sem_name, "sem_bc_q_%08x", (int32_t) (int64_t) this);
  sem_unlink(sem_name);
}


// reads a BGZF block
void BgzfStream::BgzfBlockCache::ReadBlock(void) {
  
  BT_ASSERT_X( m_device, "BgzfStream::ReadBlock() - trying to read from null IO device");
  
  // store block's starting address
  const int64_t blockAddress = m_device->Tell();
  
  // read block header from file
  char header[Constants::BGZF_BLOCK_HEADER_LENGTH];
  int64_t numBytesRead = m_device->Read(header, Constants::BGZF_BLOCK_HEADER_LENGTH);
  
  // check for device error
  if ( numBytesRead < 0 ) {
    const string message = string("device error: ") + m_device->GetErrorString();
    throw BamException("BgzfStream::ReadBlock", message);
  }
  
  // if block header empty
  if ( numBytesRead == 0 ) {
    m_blockLength = 0;
    return;
  }
  
  // if block header invalid size
  if ( numBytesRead != static_cast<int8_t>(Constants::BGZF_BLOCK_HEADER_LENGTH) )
    throw BamException("BgzfStream::ReadBlock", "invalid block header size");
  
  // validate block header contents
  if ( !BgzfStream::CheckBlockHeader(header) )
    throw BamException("BgzfStream::ReadBlock", "invalid block header contents");
  
  // copy header contents to compressed buffer
  const size_t blockLength = BamTools::UnpackUnsignedShort(&header[16]) + 1;
  memcpy(m_compressedBlock.Buffer, header, Constants::BGZF_BLOCK_HEADER_LENGTH);
  
  // read remainder of block
  const size_t remaining = blockLength - Constants::BGZF_BLOCK_HEADER_LENGTH;
  numBytesRead = m_device->Read(&m_compressedBlock.Buffer[Constants::BGZF_BLOCK_HEADER_LENGTH], remaining);
  
  // check for device error
  if ( numBytesRead < 0 ) {
    const string message = string("device error: ") + m_device->GetErrorString();
    throw BamException("BgzfStream::ReadBlock", message);
  }
  
  // check that we read in expected numBytes
  if ( numBytesRead != static_cast<int64_t>(remaining) )
    throw BamException("BgzfStream::ReadBlock", "could not read data from block");
  
  // decompress block data
  const size_t newBlockLength = InflateBlock(blockLength);
  
  // update block data
  if ( m_blockLength != 0 )
    m_blockOffset = 0;
  m_blockAddress = blockAddress;
  m_blockLength  = newBlockLength;
}

// decompresses the current block
size_t BgzfStream::BgzfBlockCache::InflateBlock(const size_t& blockLength) {
  
  // setup zlib stream object
  z_stream zs;
  zs.zalloc    = NULL;
  zs.zfree     = NULL;
  zs.next_in   = (Bytef*)m_compressedBlock.Buffer + 18;
  zs.avail_in  = blockLength - 16;
  zs.next_out  = (Bytef*)m_uncompressedBlock.Buffer;
  zs.avail_out = Constants::BGZF_DEFAULT_BLOCK_SIZE;
  
  // initialize
  int status = inflateInit2(&zs, Constants::GZIP_WINDOW_BITS);
  if ( status != Z_OK )
    throw BamException("BgzfStream::InflateBlock", "zlib inflateInit failed");
  
  // decompress
  status = inflate(&zs, Z_FINISH);
  if ( status != Z_STREAM_END ) {
    inflateEnd(&zs);
    perror("Inflate failed");
    throw BamException("BgzfStream::InflateBlock", "zlib inflate failed");
  }
  
  // finalize
  status = inflateEnd(&zs);
  if ( status != Z_OK ) {
    inflateEnd(&zs);
    throw BamException("BgzfStream::InflateBlock", "zlib inflateEnd failed");
  }
  
  // return result
  return zs.total_out;
}

int64_t BgzfStream::BgzfBlockCache::Tell()
{
  return m_tell;
}
  
bool BgzfStream::BgzfBlockCache::LoadBlock()
{
  if(!cache_error)
  {
    if( 0 != sem_wait(wait_for_copy))
      perror("Copy semaphore wait failed");
    
    BgzfBlockCacheBlock * block = block_queue.pop();
    if(0 != sem_post(wait_for_read))
      perror("Read semaphore post failed");
    
    //fprintf(stderr, "%llx %llx:%llx %d:%d", (uint64_t)m_device, m_stream->m_blockAddress, block->starting_blockAddress,  m_stream->m_blockOffset, block->starting_blockOffset);
    assert(block);
    
    //fprintf(stderr, "Get  %2x: %8lld:%5d  %8lld:%5d Len: %5d\n", (int)((int64_t) m_stream & 0xff), block->starting_blockAddress, block->starting_blockOffset, block->m_blockAddress, block->m_blockOffset, m_blockLength);
    
    //ensure that we are looking at the correct block
    if(block->starting_blockAddress != m_stream->m_blockAddress)
    {
      cerr << "Error! wrong address detected for read. Falling back to normal reads." << endl;
      cache_error = true;
      thread_exit = true;
      pthread_join(worker_thread, NULL);
      
      m_stream->m_device->Seek(m_tell);
    }
    else
    {
      //now copy out block data
      m_stream->m_blockOffset = block->m_blockOffset;
      m_stream->m_blockLength = block->m_blockLength;
      m_stream->m_blockAddress = block->m_blockAddress;
      m_tell = block->m_tell;
      memcpy(m_stream->m_uncompressedBlock.Buffer, block->buffer.Buffer, block->m_blockLength);
    }

    delete block;
  }
  
  //when an error occurs with the caching system, we need to fall
  // back to a normal, on-demand reading mechanism. This could occur
  // when a seek happens, or something else unexpected.
  if(cache_error)
  {
    //temporary- load block on demand
    m_device = m_stream->m_device;
    m_blockOffset = m_stream->m_blockOffset;
    m_blockLength = m_stream->m_blockLength;
    m_blockAddress = m_stream->m_blockAddress;
    
    //cerr << "Read " << (int64_t) m_blockAddress << " : " << (int64_t) m_blockOffset << endl;
    
    ReadBlock();
    
    m_stream->m_blockOffset = m_blockOffset;
    m_stream->m_blockLength = m_blockLength;
    m_stream->m_blockAddress = m_blockAddress;
    m_tell = m_blockAddress;
    memcpy(m_stream->m_uncompressedBlock.Buffer, m_uncompressedBlock.Buffer, m_blockLength);
  }
  
  if ( m_stream->m_blockLength != 0 )
    m_stream->m_blockOffset = 0;
  
  return true;
}

void * BgzfStream::BgzfBlockCache::thread_start(void * block_cache)
{
#ifdef __linux__
  prctl(PR_SET_NAME,"bt_bgzf_decompress",0,0,0);
#endif

  BgzfBlockCache * cache = (BgzfBlockCache*) block_cache;

  while(true)
  {
    if(0 != sem_wait(cache->wait_for_read))
      perror("Read semaphore wait failed");

    if(cache->thread_exit)
      break;

    BgzfBlockCacheBlock * block = new BgzfBlockCacheBlock();

    cache->m_device = cache->m_stream->m_device;

    block->starting_blockOffset = cache->m_blockOffset;
    block->starting_blockAddress = cache->m_blockAddress;

    cache->ReadBlock();

    block->m_blockOffset = cache->m_blockOffset;
    block->m_blockLength = cache->m_blockLength;
    block->m_blockAddress = cache->m_blockAddress;
    block->m_tell = cache->m_blockAddress;
    memcpy(block->buffer.Buffer, cache->m_uncompressedBlock.Buffer, cache->m_uncompressedBlock.NumBytes);

    //fprintf(stderr, "Read %2x: %8lld:%5d  %8lld:%5d Len: %5d\n", (int)((int64_t) cache->m_stream & 0xff), block->starting_blockAddress, block->starting_blockOffset, block->m_blockAddress, block->m_blockOffset, cache->m_blockLength);

    cache->block_queue.push(block);

    if(0 != sem_post(cache->wait_for_copy))
      perror("Read semaphore post failed");
  }
  
  return NULL;
}

// ---------------------------
// BgzfStream implementation
// ---------------------------

// constructor
BgzfStream::BgzfStream(void)
  : m_blockLength(0)
  , m_blockOffset(0)
  , m_blockAddress(0)
  , m_writeSemaphoreCount(0)
  , m_compressionLevel(Z_DEFAULT_COMPRESSION)
  , m_isWriteCompressed(true)
  , m_device(0)
  , m_uncompressedBlock(Constants::BGZF_DEFAULT_BLOCK_SIZE)
  , m_compressedBlock(Constants::BGZF_MAX_BLOCK_SIZE)
  , m_thread_pool(NULL)
  , m_lastWriteLock(NULL)
  , m_cache(NULL)
{ }

// destructor
BgzfStream::~BgzfStream(void) {
    Close();

    if(m_lastWriteLock)
        pthread_mutex_destroy(m_lastWriteLock);
    delete m_lastWriteLock;
    m_lastWriteLock = NULL;

}

// checks BGZF block header
bool BgzfStream::CheckBlockHeader(char* header) {
    return (header[0] == Constants::GZIP_ID1 &&
            header[1] == Constants::GZIP_ID2 &&
            header[2] == Z_DEFLATED &&
            (header[3] & Constants::FLG_FEXTRA) != 0 &&
            BamTools::UnpackUnsignedShort(&header[10]) == Constants::BGZF_XLEN &&
            header[12] == Constants::BGZF_ID1 &&
            header[13] == Constants::BGZF_ID2 &&
            BamTools::UnpackUnsignedShort(&header[14]) == Constants::BGZF_LEN );
}

// closes BGZF file
void BgzfStream::Close(void) {

    // skip if no device open
    if ( m_device == 0 ) return;

    // if writing to file, flush the current BGZF block,
    // then write an empty block (as EOF marker)
    if ( m_device->IsOpen() && (m_device->Mode() == IBamIODevice::WriteOnly) ) {
        FlushBlock();
        if(0 != pthread_mutex_lock(m_lastWriteLock))
            perror("Error waiting for m_lastWriteLock to close file");
        const size_t blockLength = DeflateBlock(0);
        m_device->Write(m_compressedBlock.Buffer, blockLength);
    } else {
    if(m_lastWriteLock && 0 != pthread_mutex_lock(m_lastWriteLock))
        perror("Error waiting for m_lastWriteLock to close file");
    }
    
    if(m_lastWriteLock && 0 != pthread_mutex_unlock(m_lastWriteLock))
        perror("Error waiting for m_lastWriteLock unlock to close file");

    
    if(m_lastWriteLock)
        pthread_mutex_destroy(m_lastWriteLock);
    delete m_lastWriteLock;
    m_lastWriteLock = NULL;
    
    if(m_cache) delete m_cache;

    // close device
    m_device->Close();
    delete m_device;
    m_device = 0;

    // ensure our buffers are cleared out
    m_uncompressedBlock.Clear();
    m_compressedBlock.Clear();

    // reset state
    m_blockLength = 0;
    m_blockOffset = 0;
    m_blockAddress = 0;
    m_isWriteCompressed = true;
}

// compresses the current block
size_t BgzfStream::DeflateBlock(int32_t blockLength) {

    // initialize the gzip header
    char* buffer = m_compressedBlock.Buffer;
    memset(buffer, 0, 18);
    buffer[0]  = Constants::GZIP_ID1;
    buffer[1]  = Constants::GZIP_ID2;
    buffer[2]  = Constants::CM_DEFLATE;
    buffer[3]  = Constants::FLG_FEXTRA;
    buffer[9]  = Constants::OS_UNKNOWN;
    buffer[10] = Constants::BGZF_XLEN;
    buffer[12] = Constants::BGZF_ID1;
    buffer[13] = Constants::BGZF_ID2;
    buffer[14] = Constants::BGZF_LEN;

    // loop to retry for blocks that do not compress enough
    int inputLength = blockLength;
    size_t compressedLength = 0;
    const unsigned int bufferSize = Constants::BGZF_MAX_BLOCK_SIZE;

    while ( true ) {

        // initialize zstream values
        z_stream zs;
        zs.zalloc    = NULL;
        zs.zfree     = NULL;
        zs.next_in   = (Bytef*)m_uncompressedBlock.Buffer;
        zs.avail_in  = inputLength;
        zs.next_out  = (Bytef*)&buffer[Constants::BGZF_BLOCK_HEADER_LENGTH];
        zs.avail_out = bufferSize -
                       Constants::BGZF_BLOCK_HEADER_LENGTH -
                       Constants::BGZF_BLOCK_FOOTER_LENGTH;

        // initialize the zlib compression algorithm
        int status = deflateInit2(&zs,
                                  m_compressionLevel,
                                  Z_DEFLATED,
                                  Constants::GZIP_WINDOW_BITS,
                                  Constants::Z_DEFAULT_MEM_LEVEL,
                                  Z_DEFAULT_STRATEGY);
        if ( status != Z_OK )
            throw BamException("BgzfStream::DeflateBlock", "zlib deflateInit2 failed");

        // compress the data
        status = deflate(&zs, Z_FINISH);

        // if not at stream end
        if ( status != Z_STREAM_END ) {

            deflateEnd(&zs);

            // there was not enough space available in buffer
            // try to reduce the input length & re-start loop
            if ( status == Z_OK ) {
                inputLength -= 1024;
                if ( inputLength < 0 )
                    throw BamException("BgzfStream::DeflateBlock", "input reduction failed");
                continue;
            }

            throw BamException("BgzfStream::DeflateBlock", "zlib deflate failed");
        }

        // finalize the compression routine
        status = deflateEnd(&zs);
        if ( status != Z_OK )
            throw BamException("BgzfStream::DeflateBlock", "zlib deflateEnd failed");

        // update compressedLength
        compressedLength = zs.total_out +
                           Constants::BGZF_BLOCK_HEADER_LENGTH +
                           Constants::BGZF_BLOCK_FOOTER_LENGTH;
        if ( compressedLength > Constants::BGZF_MAX_BLOCK_SIZE )
            throw BamException("BgzfStream::DeflateBlock", "deflate overflow");

        // quit while loop
        break;
    }

    // store the compressed length
    BamTools::PackUnsignedShort(&buffer[16], static_cast<uint16_t>(compressedLength - 1));

    // store the CRC32 checksum
    uint32_t crc = crc32(0, NULL, 0);
    crc = crc32(crc, (Bytef*)m_uncompressedBlock.Buffer, inputLength);
    BamTools::PackUnsignedInt(&buffer[compressedLength - 8], crc);
    BamTools::PackUnsignedInt(&buffer[compressedLength - 4], inputLength);

    // ensure that we have less than a block of data left
    int remaining = blockLength - inputLength;
    if ( remaining > 0 ) {
        if ( remaining > inputLength )
            throw BamException("BgzfStream::DeflateBlock", "after deflate, remainder too large");
        memcpy(m_uncompressedBlock.Buffer, m_uncompressedBlock.Buffer + inputLength, remaining);
    }

    // update block data
    m_blockOffset = remaining;

    // return result
    return compressedLength;
}

// flushes the data in the BGZF block
void BgzfStream::FlushBlock(void) {

    BT_ASSERT_X( m_device, "BgzfStream::FlushBlock() - attempting to flush to null device" );

    if(BamParallelismSettings::isMultithreadingEnabled()) {
      pthread_mutex_t * write_lock = new pthread_mutex_t;
      pthread_mutex_init(write_lock, NULL);
        
      if(0 != pthread_mutex_lock(write_lock))
        perror("Error locking for writing_stop_lock at beginning of compress job");

      BgzfBlockCompressJob * job = new BgzfBlockCompressJob(this, m_lastWriteLock, write_lock);

      m_lastWriteLock = write_lock;
        
      m_thread_pool->addJob(job);
      
      // update block data
      m_blockAddress += m_blockOffset ;
      m_blockOffset = 0;
    } else {
      // flush all of the remaining blocks
      while ( m_blockOffset > 0 ) {
          // compress the data block
          const size_t blockLength = DeflateBlock(m_blockOffset);

          // flush the data to our output device
          const int64_t numBytesWritten = m_device->Write(m_compressedBlock.Buffer, blockLength);

          // check for device error
          if ( numBytesWritten < 0 ) {
              const string message = string("device error: ") + m_device->GetErrorString();
              throw BamException("BgzfStream::FlushBlock", message);
          }

          // check that we wrote expected numBytes
          if ( numBytesWritten != static_cast<int64_t>(blockLength) ) {
              stringstream s("");
              s << "expected to write " << blockLength
                << " bytes during flushing, but wrote " << numBytesWritten;
              throw BamException("BgzfStream::FlushBlock", s.str());
          }

          // update block data
          m_blockAddress += blockLength;
        }
    }
}

// decompresses the current block
size_t BgzfStream::InflateBlock(const size_t& blockLength) {

    // setup zlib stream object
    z_stream zs;
    zs.zalloc    = NULL;
    zs.zfree     = NULL;
    zs.next_in   = (Bytef*)m_compressedBlock.Buffer + 18;
    zs.avail_in  = blockLength - 16;
    zs.next_out  = (Bytef*)m_uncompressedBlock.Buffer;
    zs.avail_out = Constants::BGZF_DEFAULT_BLOCK_SIZE;

    // initialize
    int status = inflateInit2(&zs, Constants::GZIP_WINDOW_BITS);
    if ( status != Z_OK )
        throw BamException("BgzfStream::InflateBlock", "zlib inflateInit failed");

    // decompress
    status = inflate(&zs, Z_FINISH);
    if ( status != Z_STREAM_END ) {
        inflateEnd(&zs);
        throw BamException("BgzfStream::InflateBlock", "zlib inflate failed");
    }

    // finalize
    status = inflateEnd(&zs);
    if ( status != Z_OK ) {
        inflateEnd(&zs);
        throw BamException("BgzfStream::InflateBlock", "zlib inflateEnd failed");
    }

    // return result
    return zs.total_out;
}

bool BgzfStream::IsOpen(void) const {
    if ( m_device == 0 )
        return false;
    return m_device->IsOpen();
}

void BgzfStream::Open(const string& filename, const IBamIODevice::OpenMode mode) {

    // close current device if necessary
    Close();
    BT_ASSERT_X( (m_device == 0), "BgzfStream::Open() - unable to properly close previous IO device" );

    if(BamParallelismSettings::isMultithreadingEnabled())
        m_thread_pool = new BamThreadPool();
	
    // retrieve new IO device depending on filename
    m_device = BamDeviceFactory::CreateDevice(filename);
    BT_ASSERT_X( m_device, "BgzfStream::Open() - unable to create IO device from filename" );

    // if device fails to open
    if ( !m_device->Open(mode) ) {
        const string deviceError = m_device->GetErrorString();
        const string message = string("could not open BGZF stream: \n\t") + deviceError;
        throw BamException("BgzfStream::Open", message);
    }
  
    if(mode == IBamIODevice::ReadOnly && BamParallelismSettings::isMultithreadingEnabled())
      m_cache = new BgzfBlockCache(this);
    
    m_lastWriteLock = new pthread_mutex_t;
    if(0 != pthread_mutex_init(m_lastWriteLock, NULL))
        perror("Error creating first bgzf write lock");
}

// reads BGZF data into a byte buffer
size_t BgzfStream::Read(char* data, const size_t dataLength) {
    if ( dataLength == 0 )
        return 0;

    // if stream not open for reading
    BT_ASSERT_X( m_device, "BgzfStream::Read() - trying to read from null device");
    if ( !m_device->IsOpen() || (m_device->Mode() != IBamIODevice::ReadOnly) )
        return 0;

    // read blocks as needed until desired data length is retrieved
    char* output = data;
    size_t numBytesRead = 0;
    while ( numBytesRead < dataLength ) {

        // determine bytes available in current block
        int bytesAvailable = m_blockLength - m_blockOffset;

        // read (and decompress) next block if needed
        if ( bytesAvailable <= 0 ) {
            if(m_cache)
                m_cache->LoadBlock();
            else
                ReadBlock();
            bytesAvailable = m_blockLength - m_blockOffset;
            if ( bytesAvailable <= 0 )
                break;
        }

        // copy data from uncompressed source buffer into data destination buffer
        const size_t copyLength = min( (dataLength-numBytesRead), (size_t)bytesAvailable );
        memcpy(output, m_uncompressedBlock.Buffer + m_blockOffset, copyLength);

        // update counters
        m_blockOffset += copyLength;
        output        += copyLength;
        numBytesRead  += copyLength;
    }

    // update block data
    if ( m_blockOffset == m_blockLength ) {
        m_blockAddress = m_cache ? m_cache->Tell() : m_device->Tell();
        m_blockOffset  = 0;
        m_blockLength  = 0;
    }

    // return actual number of bytes read
    return numBytesRead;
}

// reads a BGZF block
void BgzfStream::ReadBlock(void) {
  
    BT_ASSERT_X( !m_cache, "BgzfStream::Read() - trying to read when cache exists");
    BT_ASSERT_X( m_device, "BgzfStream::ReadBlock() - trying to read from null IO device");

    // store block's starting address
    const int64_t blockAddress = m_device->Tell();

    // read block header from file
    char header[Constants::BGZF_BLOCK_HEADER_LENGTH];
    int64_t numBytesRead = m_device->Read(header, Constants::BGZF_BLOCK_HEADER_LENGTH);

    // check for device error
    if ( numBytesRead < 0 ) {
        const string message = string("device error: ") + m_device->GetErrorString();
        throw BamException("BgzfStream::ReadBlock", message);
    }

    // if block header empty
    if ( numBytesRead == 0 ) {
        m_blockLength = 0;
        return;
    }

    // if block header invalid size
    if ( numBytesRead != static_cast<int8_t>(Constants::BGZF_BLOCK_HEADER_LENGTH) )
        throw BamException("BgzfStream::ReadBlock", "invalid block header size");

    // validate block header contents
    if ( !BgzfStream::CheckBlockHeader(header) )
        throw BamException("BgzfStream::ReadBlock", "invalid block header contents");

    // copy header contents to compressed buffer
    const size_t blockLength = BamTools::UnpackUnsignedShort(&header[16]) + 1;
    memcpy(m_compressedBlock.Buffer, header, Constants::BGZF_BLOCK_HEADER_LENGTH);

    // read remainder of block
    const size_t remaining = blockLength - Constants::BGZF_BLOCK_HEADER_LENGTH;
    numBytesRead = m_device->Read(&m_compressedBlock.Buffer[Constants::BGZF_BLOCK_HEADER_LENGTH], remaining);

    // check for device error
    if ( numBytesRead < 0 ) {
        const string message = string("device error: ") + m_device->GetErrorString();
        throw BamException("BgzfStream::ReadBlock", message);
    }

    // check that we read in expected numBytes
    if ( numBytesRead != static_cast<int64_t>(remaining) )
        throw BamException("BgzfStream::ReadBlock", "could not read data from block");

    // decompress block data
    const size_t newBlockLength = InflateBlock(blockLength);

    // update block data
    if ( m_blockLength != 0 )
        m_blockOffset = 0;
    m_blockAddress = blockAddress;
    m_blockLength  = newBlockLength;
}

// seek to position in BGZF file
void BgzfStream::Seek(const int64_t& position) {

    BT_ASSERT_X( m_device, "BgzfStream::Seek() - trying to seek on null IO device");

    // skip if device is not open
    if ( !IsOpen() ) return;

    // determine adjusted offset & address
    int     blockOffset  = (position & 0xFFFF);
    int64_t blockAddress = (position >> 16) & 0xFFFFFFFFFFFFLL;

    // attempt seek in file
    if ( m_device->IsRandomAccess() && m_device->Seek(blockAddress) ) {

        // update block data & return success
        m_blockLength  = 0;
        m_blockAddress = blockAddress;
        m_blockOffset  = blockOffset;
    }
    else {
        stringstream s("");
        s << "unable to seek to position: " << position;
        throw BamException("BgzfStream::Seek", s.str());
    }
}

void BgzfStream::SetWriteCompressed(bool ok) {
  m_isWriteCompressed = ok;
  if(ok)
    m_compressionLevel = Z_DEFAULT_COMPRESSION;
  else
    m_compressionLevel = Z_NO_COMPRESSION;
}

void BgzfStream::SetCompressionLevel(int compressionLevel) {
  m_compressionLevel = compressionLevel;
  m_isWriteCompressed = (compressionLevel != Z_NO_COMPRESSION);
}

// get file position in BGZF file
int64_t BgzfStream::Tell(void) const {
    if ( !IsOpen() )
        return 0;
    return ( (m_blockAddress << 16) | (m_blockOffset & 0xFFFF) );
}

// writes the supplied data into the BGZF buffer
size_t BgzfStream::Write(const char* data, const size_t dataLength) {

    BT_ASSERT_X( m_device, "BgzfStream::Write() - trying to write to null IO device");
    BT_ASSERT_X( (m_device->Mode() == IBamIODevice::WriteOnly),
                 "BgzfStream::Write() - trying to write to non-writable IO device");

    // skip if file not open for writing
    if ( !IsOpen() )
        return 0;

    // write blocks as needed til all data is written
    size_t numBytesWritten = 0;
    const char* input = data;
    const size_t blockLength = Constants::BGZF_DEFAULT_BLOCK_SIZE;
    while ( numBytesWritten < dataLength ) {

        // copy data contents to uncompressed output buffer
        unsigned int copyLength = min(blockLength - m_blockOffset, dataLength - numBytesWritten);
        char* buffer = m_uncompressedBlock.Buffer;
        memcpy(buffer + m_blockOffset, input, copyLength);

        // update counter
        m_blockOffset   += copyLength;
        input           += copyLength;
        numBytesWritten += copyLength;

        // flush (& compress) output buffer when full
        if ( m_blockOffset == static_cast<int32_t>(blockLength) )
            FlushBlock();
    }

    // return actual number of bytes written
    return numBytesWritten;
}
