// ***************************************************************************
// bamtools_mergesort.cpp (c) 2010 Derek Barnett, Erik Garrison, Lee C. Baker
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 6 April 2012
// ---------------------------------------------------------------------------
// Merges multiple BAM files into one sorted BAM file
//
// This file contains the bulk of SortTool and MergeTool, running in separate
// threads and paired together with a FIFO. To merge and sort files, the following steps
// are taken:
// # Files are read in sequentially, concatenated, and fed to a FIFO (Thread #1)
// # The merged file is read from the FIFO, split into chunks, and passed to worker threads. (Thread #2)
// # Worker threads compress (optional) and then write these sorted chunks to disk. (Thread #3,4...N)
// # The sorted chunks on disk are combined into one final sorted file. (Thread #2)
//
// TODO:
// * SORT_DEFAULT_MAX_BUFFER_MEMORY is currently unused
// * Should read in the files to be merged separate threads
// * Document functions better
// ***************************************************************************

#include "bamtools_mergesort.h"

#include <api/SamConstants.h>
#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
#include <api/algorithms/Sort.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_utilities.h>
#include "api/internal/utils/BamThreadPool.h"
#include "api/BamParallelismSettings.h"
#include "api/internal/io/BgzfStream_p.h"

#ifdef __linux__
#include <sys/prctl.h>
#endif

//for mkfifo
#include <sys/types.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <semaphore.h>

#include <sys/time.h>
#include <sys/resource.h>

#include <pthread.h>

using namespace BamTools;
using namespace BamTools::Algorithms;

#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <iostream>
using namespace std;

const unsigned int SORT_DEFAULT_MAX_BUFFER_COUNT  = 500000;  // max numberOfAlignments for buffer
const unsigned int SORT_DEFAULT_MAX_BUFFER_MEMORY = 1024;    // Mb
const unsigned int MERGESORT_DEFAULT_SHOULD_COMPRESS_TEMP_FILES = 0;
const unsigned int MERGESORT_DEFAULT_NUM_THREADS = 0;
const unsigned int MERGESORT_MIN_SORT_SIZE = 30000;    //don't parallelize sort jobs smaller than this many alignments

// ---------------------------------------------
// MergeSettings implementation

struct MergeSortTool::MergeSortSettings {
  
  // flags from merge
  bool HasInputBamFilename;
  bool HasOutputBamFilename;
  bool HasCompressionLevel;
  bool HasRegion;
  bool HasMaxBufferCount;
  bool HasMaxBufferMemory;
  bool IsSortingByName;
  bool HasCompressTempFiles;
  bool HasNumThreads;
  unsigned int CompressTempFiles;
  bool HasVerbose;
  bool HasNoThreads;
  
  // filenames
  vector<string> InputFiles;
  
  // other parameters
  string OutputFilename;
  string Region;
  unsigned int MaxBufferCount;
  unsigned int MaxBufferMemory;
  unsigned int NumThreads;
  unsigned int CompressionLevel;
  
  // constructor
  MergeSortSettings(void)
  : HasInputBamFilename(false)
  , HasOutputBamFilename(false)
  , HasCompressionLevel(false)
  , HasRegion(false)
  , HasMaxBufferCount(false)
  , HasMaxBufferMemory(false)
  , IsSortingByName(false)
  , HasCompressTempFiles(false)
  , HasNumThreads(false)
  , CompressTempFiles(true)
  , HasVerbose(false)
  , HasNoThreads(false)
  
  //, OutputBamFilename(Options::StandardOut())
  , MaxBufferCount(SORT_DEFAULT_MAX_BUFFER_COUNT)
  , MaxBufferMemory(SORT_DEFAULT_MAX_BUFFER_MEMORY)
  , NumThreads(-1)
  , CompressionLevel(0)
  { }
};  

// ---------------------------------------------
// MergeToolPrivate implementation

struct MergeSortTool::MergeSortToolPrivate {
public:
  MergeSortToolPrivate(MergeSortTool::MergeSortSettings* settings)
  : m_tempFilenameStub("bam_mergesort_temp_")
  , m_numberOfAlignments(0)
  , m_settings(settings)
  , use_spinlocks(false)
  , alignment_available_spinlock(0)
  , alignment_consumed_spinlock(0)
  { 
      use_spinlocks = true; BamParallelismSettings::availableCores() > 6;
  }

  ~MergeSortToolPrivate(void) { }

  // interface
public:
  bool Run(void);
  bool RunSort(void);
  bool RunMerge(void);

  static void * RunMergeThread(void *);
  static void * RunSortThread(void *);

  //from SortTool:
  bool CreateSortedTempFile(vector<BamAlignment *> * buffer);
  bool GenerateSortedRuns(void);
  bool MergeSortedRuns(void);
  bool WriteTempFile(const vector<BamAlignment> & buffer, const string& tempFilename);
  bool WriteTempFile(const vector<BamAlignment *> & buffer, const string& tempFilename);
  template<class T> void SortBuffer(vector<T>& buffer);

  // data members
private:
    BamAlignment * GetMergedAlignment();
    void PostMergedAlignment(BamAlignment * al);
  
  string m_tempFilenameStub;
  int m_numberOfRuns;
  int64_t m_numberOfAlignments;
  SamHeader m_header;
  RefVector m_references;
  vector<string> m_tempFilenames;
  bool sort_retval, merge_retval;
  BamThreadPool * thread_pool;
  BamThreadPool * sort_thread_pool;
  MergeSortTool::MergeSortSettings* m_settings;
  SynchronizedQueue <BamAlignment *> merged_alignment_queue;
  sem_t * alignment_available; 
  bool use_spinlocks;
  
  // for performance reasons, we may use these for lockless coordination of access to merged_alignment_queue
  volatile int64_t alignment_available_spinlock, alignment_consumed_spinlock; 
public:
  class TempFileWriteJob : public BamThreadJob
  {
  public:
    TempFileWriteJob(MergeSortTool::MergeSortToolPrivate * tool, vector<BamAlignment *> * buffer, string filename);
    virtual void runJob();
  protected:
    string filename;
    vector<BamAlignment *> * buffer;
    MergeSortTool::MergeSortToolPrivate * tool;
  };
    
  template<class T>
  class SortJob : public BamThreadJob
  {
  public:
    SortJob(typename vector<T>::iterator begin,typename vector<T>::iterator end, pthread_mutex_t & completion_lock,const MergeSortSettings * settings);
    virtual void runJob();
    protected:
    typename vector<T>::iterator begin;
    typename vector<T>::iterator end;
    pthread_mutex_t & completion_lock;
    const MergeSortSettings * m_settings;
  };
};

void * MergeSortTool::MergeSortToolPrivate::RunMergeThread(void * data)
{
#ifdef __linux__
  prctl(PR_SET_NAME,"bt_merge",0,0,0);
#endif
  MergeSortTool::MergeSortToolPrivate * tool = (MergeSortTool::MergeSortToolPrivate *)data;
  tool->merge_retval = tool->RunMerge();

  return NULL;
}

void * MergeSortTool::MergeSortToolPrivate::RunSortThread(void * data)
{
#ifdef __linux__
  prctl(PR_SET_NAME,"bt_sort",0,0,0);
#endif

  MergeSortTool::MergeSortToolPrivate * tool = (MergeSortTool::MergeSortToolPrivate *)data;
  tool->sort_retval = tool->RunSort();

  return NULL;
}

const char * ALIGNMENT_QUEUE_NAME = "mergesort_align_queue";

//actual run creates threads for other 'run'commands
bool MergeSortTool::MergeSortToolPrivate::Run(void)
{
  timeval start_time;
  gettimeofday(&start_time, NULL);

  if(!m_settings->HasNoThreads) {
    BamParallelismSettings::enableMultithreading();
    if(m_settings->HasNumThreads)
        BamParallelismSettings::setNumberThreads(m_settings->NumThreads);
    thread_pool = new BamThreadPool();
    sort_thread_pool = new BamThreadPool();
    if(m_settings->HasVerbose) {
      if(m_settings->HasNumThreads)
        cerr << "Specified " << m_settings->NumThreads << " cores for use in thread pool." << endl;
      else
        cerr << "Detected " << BamParallelismSettings::availableCores() << " cores for use in thread pool." << endl;
    }
  } else {
    BamParallelismSettings::disableMultithreading();
    if(m_settings->HasVerbose)
      cerr << "Thread pool use disabled." << endl;
  }

  // if compression level isn't set, then use 0 compression for stdout, and 6 when
  // writing to files.
  if(!m_settings->HasCompressionLevel) {
    if(m_settings->HasOutputBamFilename)
      m_settings->CompressionLevel = 6;
    else
      m_settings->CompressionLevel = 0;
  }

  sem_unlink(ALIGNMENT_QUEUE_NAME);
  alignment_available = sem_open(ALIGNMENT_QUEUE_NAME, O_CREAT | O_EXCL, 0, 0);
  if(alignment_available == 0)
    perror("Failed creating alignment queue semaphore.");

  // We use a fifo to communicate data between the merge thread and the sort thread.
  // This enables us to use much of the merge / sort code as is.

  //create threads to do the work
  pthread_t merge_thread, sort_thread;
  pthread_create(&merge_thread, NULL, MergeSortTool::MergeSortToolPrivate::RunMergeThread, this);
  pthread_create(&sort_thread, NULL, MergeSortTool::MergeSortToolPrivate::RunSortThread, this);

  //wait for both threads to finish
  pthread_join(merge_thread, NULL);
  pthread_join(sort_thread, NULL);

  if(0 != sem_close(alignment_available))
    perror("Error closing alignment_available semaphore");

  if(0 != sem_unlink(ALIGNMENT_QUEUE_NAME))
    perror("Error unlinking alignment queue");

  if(!m_settings->HasNoThreads) {
    delete thread_pool;
    delete sort_thread_pool;
  }

  if(m_settings->HasVerbose)
  {
    //print cpu usage and time
    rusage r;
    getrusage(RUSAGE_SELF, &r);
    timeval stop_time;
    gettimeofday(&stop_time, NULL);
    timeval real_time;
    real_time.tv_sec = stop_time.tv_sec - start_time.tv_sec;
    real_time.tv_usec = stop_time.tv_usec - start_time.tv_usec;
    fprintf(stderr, "Elapsed time: %3ldm%06.3fs\n", real_time.tv_sec /60, float(real_time.tv_sec %60) + (1.e-6 * real_time.tv_usec) );
    fprintf(stderr, "CPU time: %3ldm%06.3fs (user) / %3ldm%06.3fs (sys)\n", r.ru_utime.tv_sec /60, float(r.ru_utime.tv_sec %60) + (1.e-6 * r.ru_utime.tv_usec), r.ru_stime.tv_sec /60, float(r.ru_stime.tv_sec %60) + (1.e-6 * r.ru_stime.tv_usec));
#ifndef __linux__
    // on linux, the max ram consumption is in KB. On Mac, it is in bytes.
    r.ru_maxrss /= 1024;
#endif
    fprintf(stderr, "Max mem: %6ld MB\n", r.ru_maxrss /1024);
  }

  return sort_retval && merge_retval;
}

// This function is mostly copied from SortToolPrivate in bamtools_sort.c
bool MergeSortTool::MergeSortToolPrivate::RunMerge(void) {

  // set to default input if none provided
  if ( !m_settings->HasInputBamFilename )
    m_settings->InputFiles.push_back(Options::StandardIn());

  if(m_settings->HasVerbose)
    cerr << "Merging " << m_settings->InputFiles.size() << " input files." << endl;

  // opens the BAM files (by default without checking for indexes)
  BamMultiReader reader;
  if ( !reader.Open(m_settings->InputFiles) ) {
    cerr << "bamtools merge ERROR: could not open input BAM file(s)... Aborting." << endl;
    return false;
  }

  //save header information for temp files and final merge
  m_header = reader.GetHeader();
  m_references = reader.GetReferenceData();

  //the first semaphore post/wait signals header information is valid
  if(0 != sem_post(alignment_available))
    perror("Error posting alignment_available (open)");

  // if no region specified, store entire contents of file(s)
  if ( !m_settings->HasRegion ) {
    BamAlignment * al;
    while ( NULL != (al = reader.GetNextAlignmentCore()) ) {
        PostMergedAlignment(al);
    }
  }

  // otherwise attempt to use region as constraint
  else {
    // if region string parses OK
    BamRegion region;
    if ( Utilities::ParseRegionString(m_settings->Region, reader, region) ) {

      // attempt to find index files
      reader.LocateIndexes();

      // if index data available for all BAM files, we can use SetRegion
      if ( reader.HasIndexes() ) {

        // attempt to use SetRegion(), if failed report error
        if ( !reader.SetRegion(region.LeftRefID,
                               region.LeftPosition,
                               region.RightRefID,
                               region.RightPosition) )
        {
          cerr << "bamtools merge ERROR: set region failed. Check that REGION describes a valid range"
          << endl;
          reader.Close();
          return false;
        }

        // everything checks out, just iterate through specified region, storing alignments
        BamAlignment * al;
        while ( NULL != (al = reader.GetNextAlignmentCore()) ) {
            PostMergedAlignment(al);
        }
      }

      // no index data available, we have to iterate through until we
      // find overlapping alignments
      else {
        BamAlignment * al;
        while ( NULL != (al = reader.GetNextAlignmentCore()) ) {
          if ( (al->RefID >= region.LeftRefID)  && ( (al->Position + al->Length) >= region.LeftPosition ) &&
              (al->RefID <= region.RightRefID) && ( al->Position <= region.RightPosition) ) {
              PostMergedAlignment(al);
          }
        }
      }
    }

    // error parsing REGION string
    else {
      cerr << "bamtools merge ERROR: could not parse REGION - " << m_settings->Region << endl;
      cerr << "Check that REGION is in valid format (see documentation) and that the coordinates are valid"
      << endl;
      reader.Close();
      return false;
    }
  }

  //post with an empty queue to signal EOF to the sorting thread
  if(use_spinlocks)
    alignment_available_spinlock++;
  else
    if(0 != sem_post(alignment_available))
      perror("Error posting alignment_available (close)");

  // clean & exit
  reader.Close();

  return true;
}

void MergeSortTool::MergeSortToolPrivate::PostMergedAlignment( BamAlignment * al)
{
    merged_alignment_queue.push(al);
    if(use_spinlocks)
        alignment_available_spinlock++;
    else
        if(0 != sem_post(alignment_available))
            perror("Error posting alignment_available");
    
    //if the queue gets too big, our memory usage may get out of control.
    if(alignment_available_spinlock % 4000 == 0)    //don't check often
        if(merged_alignment_queue.size() > 40000)
            usleep(20000);
}

BamAlignment * MergeSortTool::MergeSortToolPrivate::GetMergedAlignment()
{
  if(use_spinlocks) {
    //give up our timeslice if the queue is empty
    while(alignment_available_spinlock == alignment_consumed_spinlock) sched_yield();
    alignment_consumed_spinlock++;
  } else
    if(0 != sem_wait(alignment_available))
      perror("Error waiting for alignment_available semaphore");

  bool ret = !merged_alignment_queue.empty();
  if(ret)
    return merged_alignment_queue.pop();
  else
    return NULL;
}

// generates mutiple sorted temp BAM files from single unsorted BAM file
bool MergeSortTool::MergeSortToolPrivate::GenerateSortedRuns(void) {

  if(m_settings->HasVerbose)
    cerr << "Generating sorted temp files...";

  //wait for header information
  if(0 != sem_wait(alignment_available))
    perror("Error waiting for alignment_available semaphore");

  // get basic data that will be shared by all temp/output files
  m_header.SortOrder = ( m_settings->IsSortingByName
                      ? Constants::SAM_HD_SORTORDER_QUERYNAME
                      : Constants::SAM_HD_SORTORDER_COORDINATE );

  // set up alignments buffer
  vector<BamAlignment *> * buffer = new vector<BamAlignment *>();
  buffer->reserve( (size_t)(m_settings->MaxBufferCount*1.1) );
  bool bufferFull = false;
  BamAlignment * al = NULL;

  // if sorting by name, we need to generate full char data
  // so can't use GetNextAlignmentCore()
  if ( m_settings->IsSortingByName ) {

    // iterate through file
    while (true) {
      al = GetMergedAlignment();
      if(!al)
        break;

      // check buffer's usage
      bufferFull = ( buffer->size() >= m_settings->MaxBufferCount );

      // store alignments until buffer is "full"
      if ( !bufferFull )
        buffer->push_back(al);

      // if buffer is "full"
      else {
        // so create a sorted temp file with current buffer contents
        // then push "al" into fresh buffer
        CreateSortedTempFile(buffer);
        buffer = new vector<BamAlignment *>();
        buffer->reserve( (size_t)(m_settings->MaxBufferCount*1.1) );
        buffer->push_back(al);
      }
    }
  }

  // sorting by position, can take advantage of GNACore() speedup
  else {
    // iterate through file
    while (true) {
      al = GetMergedAlignment();
      if(!al)
        break;

      // check buffer's usage
      bufferFull = ( buffer->size() >= m_settings->MaxBufferCount );

      // store alignments until buffer is "full"
      if ( !bufferFull )
        buffer->push_back(al);

      // if buffer is "full"
      else {
        // create a sorted temp file with current buffer contents
        // then push "al" into fresh buffer
        CreateSortedTempFile(buffer);
        buffer = new vector<BamAlignment *>();
        buffer->push_back(al);
          cerr << "." ;
      }
    }
  }

  // handle any leftover buffer contents
  if ( !buffer->empty() )
    CreateSortedTempFile(buffer);

  if(m_settings->HasVerbose)
    cerr << "waiting for files to be compressed / written...";

  //wait for all temp files to be created in other threads
  if(!m_settings->HasNoThreads)
    thread_pool->waitForJobCompletion();

  if(m_settings->HasVerbose)
    cerr << "done." << endl;

  return true;
}

MergeSortTool::MergeSortToolPrivate::TempFileWriteJob::TempFileWriteJob(MergeSortTool::MergeSortToolPrivate * tool, vector<BamAlignment *> * buffer, string filename) :
filename(filename), buffer(buffer), tool(tool)
{
}

void MergeSortTool::MergeSortToolPrivate::TempFileWriteJob::runJob()
{
#ifdef __linux__
    prctl(PR_SET_NAME,"bt_temp_sort",0,0,0);
#endif
    // do sorting
    tool->SortBuffer(*buffer);
    
#ifdef __linux__
    prctl(PR_SET_NAME,"bt_temp_write",0,0,0);
#endif
  // as noted in the comments of the original file, success is never
  // used so we never return it.
  bool success = tool->WriteTempFile( *buffer, filename );
  if(!success)
    cerr << "Problem writing out temporary file " << filename;

#ifdef __linux__
    prctl(PR_SET_NAME,"bt_temp_cleanup",0,0,0);
#endif
  for(size_t i = 0; i < buffer->size(); i++)
    delete (*buffer)[i];
  delete buffer;
}

bool MergeSortTool::MergeSortToolPrivate::CreateSortedTempFile(vector<BamAlignment* > * buffer) {
  //make filename
  stringstream tempStr;
  tempStr << m_tempFilenameStub << m_numberOfRuns;
  m_tempFilenames.push_back(tempStr.str());

  ++m_numberOfRuns;

  bool success = true;

  if(m_settings->HasNoThreads) {
    vector<BamAlignment> sort_buffer;
    sort_buffer.reserve(buffer->size());
    for(size_t i = 0; i < buffer->size(); i++)
    {
      BamAlignment * al = buffer->at(i);
      sort_buffer.push_back(*al);
      delete al;
    }

    // do sorting
    SortBuffer(sort_buffer);

    // write sorted contents to temp file, store success/fail
    success = WriteTempFile( sort_buffer, tempStr.str() );
    delete buffer;

  } else {
    TempFileWriteJob * job = new TempFileWriteJob(this, buffer, tempStr.str());
    thread_pool->addJob(job);
  }

  // return success/fail of writing to temp file
  // TODO: a failure returned here is not actually caught and handled anywhere
  return success;
}

// merges sorted temp BAM files into single sorted output BAM file
bool MergeSortTool::MergeSortToolPrivate::MergeSortedRuns(void) {

  // In case we just have one temp file, no use in reading it in and then
  // writing the same thing back out.
  if(m_tempFilenames.size() == 1 && m_settings->HasCompressTempFiles && m_settings->CompressTempFiles) {
    if(m_settings->HasVerbose)
      cerr << "Since only one temp file was create, we will use it as the final output." << endl;

    if(0 != rename(m_tempFilenames[0].c_str(), m_settings->OutputFilename.c_str())) {
      cerr << "bamtools mergesort ERROR: could not move temporary file for output... Aborting." << endl;
      return false;
    }
  } else {
    cerr << "Combining temp files for final output...";

    // open up multi reader for all of our temp files
    // this might get broken up if we do a multi-pass system later ??
    BamMultiReader multiReader;
    if ( !multiReader.Open(m_tempFilenames) ) {
      cerr << "bamtools mergesort ERROR: could not open BamMultiReader for merging temp files... Aborting."
      << endl;
      return false;
    }

    // open writer for our completely sorted output BAM file
    BamWriter mergedWriter;
    mergedWriter.SetCompressionLevel(m_settings->CompressionLevel);

    if ( !mergedWriter.Open(m_settings->OutputFilename, m_header, m_references) ) {
      cerr << "bamtools mergesort ERROR: could not open " << m_settings->OutputFilename
      << " for writing... Aborting." << endl;
      multiReader.Close();
      return false;
    }

    // while data available in temp files
    BamAlignment al;
    int count = 0;
    while ( multiReader.GetNextAlignmentCore(al) ) {
      mergedWriter.SaveAlignment(al);
      if(count++ % 1000000 == 0)
          cerr << ".";  //progress indicator every 1M alignments written.
    }

    // close files
    multiReader.Close();
    mergedWriter.Close();
    cerr << "done." << endl << "Clearing " << m_tempFilenames.size() << " temp files...";
  }

  // delete all temp files
  vector<string>::const_iterator tempIter = m_tempFilenames.begin();
  vector<string>::const_iterator tempEnd  = m_tempFilenames.end();
  for ( ; tempIter != tempEnd; ++tempIter ) {
    const string& tempFilename = (*tempIter);
    remove(tempFilename.c_str());
  }
  cerr << "done." << endl;
  // return success
  return true;
}

bool MergeSortTool::MergeSortToolPrivate::RunSort(void) {
  // this does a single pass, chunking up the input file into smaller sorted temp files,
  // then write out using BamMultiReader to handle merging

  if ( GenerateSortedRuns() )
    return MergeSortedRuns();
  else
    return false;
}

template<class T> MergeSortTool::MergeSortToolPrivate::SortJob<T>::SortJob(typename vector<T>::iterator begin, typename vector<T>::iterator end, pthread_mutex_t & completion_lock, const MergeSortSettings * settings)
: begin(begin)
, end(end)
, completion_lock(completion_lock)
, m_settings(settings)
{
}

template<class T> void MergeSortTool::MergeSortToolPrivate::SortJob<T>::runJob()
{
#ifdef __linux__
  prctl(PR_SET_NAME,"bt_tempfile_sort",0,0,0);
#endif
  if(m_settings->IsSortingByName)
    std::stable_sort( begin, end, Sort::ByName() );
  else
    std::stable_sort( begin, end, Sort::ByPosition() );

  pthread_mutex_unlock(&completion_lock);
}

//this function is designed to accept either BamAlignment or BamAlignment* as T
template<class T>
void MergeSortTool::MergeSortToolPrivate::SortBuffer(vector<T>& buffer) {

  if(m_settings->HasNoThreads)
  {
    if(m_settings->IsSortingByName)
      std::stable_sort( buffer.begin(), buffer.end(), Sort::ByName() );
    else
      std::stable_sort( buffer.begin(), buffer.end(), Sort::ByPosition() );
  } else {
      int divisions = buffer.size() / MERGESORT_MIN_SORT_SIZE;
      if(divisions > BamParallelismSettings::availableCores())
          divisions = BamParallelismSettings::availableCores();
      if(divisions < 1)
          divisions = 1;
 
    pthread_mutex_t * locks = new pthread_mutex_t[divisions];
    SortJob<T> ** jobs = new SortJob<T> *[divisions];
    size_t section_length = buffer.size() / divisions;

    //start jobs
    for(int ctr = 0; ctr < divisions; ctr++) {
      if(0 != pthread_mutex_init(&locks[ctr], NULL)) {
        perror("Error initializing a sort mutex");
        assert(0);
      }

      if(0 != pthread_mutex_lock(&locks[ctr])) {
        perror("Error locking(1) a sort mutex");
        assert(0);
      }

      typename vector<T>::iterator begin = buffer.begin() + ctr * section_length;
      typename vector<T>::iterator end = (ctr == divisions - 1) ? buffer.end() : (buffer.begin() + (1+ctr) * section_length);
#if 1
      jobs[ctr] = new SortJob<T>(begin, end, locks[ctr], m_settings);
      sort_thread_pool->addJob(jobs[ctr]);
#else
      if(m_settings->IsSortingByName)
        std::stable_sort( begin, end, Sort::ByName() );
      else
        std::stable_sort( begin, end, Sort::ByPosition() );
#endif
    }

    if(0 != pthread_mutex_lock(&locks[0])) {
      perror("Error locking(2) a sort mutex");
      assert(0);
    }

    //now, rejoin
    for(int ctr = 1; ctr < divisions; ctr++) {
        
        if(0 != pthread_mutex_lock(&locks[ctr])) {
            perror("Error locking(3) a sort mutex");
            assert(0);
        }
        if(0 != pthread_mutex_unlock(&locks[ctr])) {
            perror("Error unlocking(3) a sort mutex");
            assert(0);
        }

      if(0 != pthread_mutex_destroy(&locks[ctr])) {
        perror("Error destroying a sort mutex");
        assert(0);
      }

      typename vector<T>::iterator midpoint = buffer.begin() + ctr * section_length;
      typename vector<T>::iterator end = (ctr == divisions - 1) ? buffer.end() : (buffer.begin() + (1+ctr) * section_length);

      if(m_settings->IsSortingByName)
        std::inplace_merge(buffer.begin(), midpoint, end, Sort::ByName() );
      else
        std::inplace_merge(buffer.begin(), midpoint, end, Sort::ByPosition() );
    }

    delete [] jobs;
    delete [] locks;
  }
}

bool MergeSortTool::MergeSortToolPrivate::WriteTempFile(const vector<BamAlignment>& buffer,
                                              const string& tempFilename)
{
  // open temp file for writing
  BamWriter tempWriter;

  if(!m_settings->HasCompressTempFiles || !m_settings->CompressTempFiles)
    tempWriter.SetCompressionMode(BamWriter::Uncompressed);
  else
    tempWriter.SetCompressionMode(BamWriter::Compressed);

  if ( !tempWriter.Open(tempFilename, m_header, m_references) ) {
    cerr << "bamtools sort ERROR: could not open " << tempFilename
    << " for writing." << endl;
    return false;
  }

  // write data
  vector<BamAlignment>::const_iterator buffIter = buffer.begin();
  vector<BamAlignment>::const_iterator buffEnd  = buffer.end();
  for ( ; buffIter != buffEnd; ++buffIter )  {
    const BamAlignment& al = (*buffIter);
    tempWriter.SaveAlignment(al);
  }

  // close temp file & return success
  tempWriter.Close();
  return true;
}

bool MergeSortTool::MergeSortToolPrivate::WriteTempFile(const vector<BamAlignment *>& buffer,
                                                        const string& tempFilename)
{
    // open temp file for writing
    BamWriter tempWriter;
    
    if(!m_settings->HasCompressTempFiles || !m_settings->CompressTempFiles)
        tempWriter.SetCompressionMode(BamWriter::Uncompressed);
    else
        tempWriter.SetCompressionMode(BamWriter::Compressed);
    
    if ( !tempWriter.Open(tempFilename, m_header, m_references) ) {
        cerr << "bamtools sort ERROR: could not open " << tempFilename
        << " for writing." << endl;
        return false;
    }
    
    // write data
    vector<BamAlignment *>::const_iterator buffIter = buffer.begin();
    vector<BamAlignment *>::const_iterator buffEnd  = buffer.end();
    for ( ; buffIter != buffEnd; ++buffIter )  {
        const BamAlignment * al = (*buffIter);
        tempWriter.SaveAlignment(*al);
    }
    
    // close temp file & return success
    tempWriter.Close();
    return true;
}


// ---------------------------------------------
// MergeTool implementation

MergeSortTool::MergeSortTool(void)
: AbstractTool()
, m_settings(new MergeSortSettings)
, m_impl(0)
{
  // set program details
  Options::SetProgramInfo("bamtools mergesort", "merges multiple BAM files into one sorted file", "[-in <filename> -in <filename> ...] [-out <filename> | [-forceCompression]] [-region <REGION>] [sortOptions]");

  // set up options
  OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");
  Options::AddValueOption("-in",  "BAM filename", "the input BAM file(s)", "", m_settings->HasInputBamFilename,  m_settings->InputFiles,     IO_Opts);
  Options::AddValueOption("-out", "BAM filename", "the output BAM file",   "", m_settings->HasOutputBamFilename, m_settings->OutputFilename, IO_Opts);
  Options::AddValueOption("-compressionLevel","0-9", "if results are sent to stdout (like when piping to another tool), default behavior is to leave output uncompressed(0); otherwise, compression is set to the default level (6)", "", m_settings->HasCompressionLevel, m_settings->CompressionLevel, IO_Opts);
  Options::AddValueOption("-region", "REGION", "genomic region. See README for more details", "", m_settings->HasRegion, m_settings->Region, IO_Opts);
  Options::AddOption("-v", "(verbose) print additional information while processing files", m_settings->HasVerbose, IO_Opts);

  OptionGroup* SortOpts = Options::CreateOptionGroup("Sorting Methods");
  Options::AddOption("-byname", "sort by alignment name", m_settings->IsSortingByName, SortOpts);

  OptionGroup* MemOpts = Options::CreateOptionGroup("Memory/CPU Settings");
  Options::AddValueOption("-n",   "count", "max number of alignments per tempfile", "",
                          m_settings->HasMaxBufferCount,  m_settings->MaxBufferCount,
                          MemOpts, SORT_DEFAULT_MAX_BUFFER_COUNT);
  Options::AddValueOption("-mem", "Mb", "max memory to use", "",
                          m_settings->HasMaxBufferMemory, m_settings->MaxBufferMemory,
                          MemOpts, SORT_DEFAULT_MAX_BUFFER_MEMORY);
  Options::AddValueOption("-compressTempFiles", "1=yes, 0=no", "temp files should be compressed (default yes)", "",
                          m_settings->HasCompressTempFiles, m_settings->CompressTempFiles,
                          MemOpts, MERGESORT_DEFAULT_SHOULD_COMPRESS_TEMP_FILES);
  Options::AddValueOption("-threads", "count", "number of threads per thread pool", "",
                          m_settings->HasNumThreads,  m_settings->NumThreads,
                          MemOpts, MERGESORT_DEFAULT_NUM_THREADS );
  Options::AddOption("-noThreads", "disable use of the thread pool for parallel processing", m_settings->HasNoThreads, MemOpts);
}

MergeSortTool::~MergeSortTool(void) {
  delete m_settings;
  m_settings = 0;

  delete m_impl;
  m_impl = 0;
}

int MergeSortTool::Help(void) {
  Options::DisplayHelp();
  return 0;
}

int MergeSortTool::Run(int argc, char* argv[]) {
  // parse command line arguments
  Options::Parse(argc, argv, 1);

  // initialize MergeTool with settings
  m_impl = new MergeSortToolPrivate(m_settings);

  // run MergeTool, return success/fail
  if ( m_impl->Run() )
    return 0;
  else
    return 1;
}

