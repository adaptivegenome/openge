#include "commands.h"
#include "../util/thread_pool.h"
#include <cassert>
using namespace std;

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

#include <api/SamConstants.h>
#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
#include <api/algorithms/Sort.h>
#include <api/BamParallelismSettings.h>

#include "../algorithms/algorithm_module_adaptor.h"
#include "../algorithms/mark_duplicates.h"
#include "../algorithms/file_writer.h"

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

namespace po = boost::program_options;

#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <iostream>
using namespace std;

const unsigned int SORT_DEFAULT_MAX_BUFFER_COUNT  = 500000;  // max numberOfAlignments for buffer
const unsigned int SORT_DEFAULT_MAX_BUFFER_MEMORY = 1024;    // Mb
const unsigned int MERGESORT_MIN_SORT_SIZE = 30000;    //don't parallelize sort jobs smaller than this many alignments


// ---------------------------------------------
// MergeToolPrivate implementation

class MergeSortCommand::MergeSortCommandImplementation {
public:
    MergeSortCommandImplementation(MergeSortCommand & command )
    : m_tempFilenameStub("bam_mergesort_temp_")
    , m_numberOfAlignments(0)
    , command(command)
    , use_spinlocks(false)
    , alignment_available_spinlock(0)
    , alignment_consumed_spinlock(0)
    { 
        use_spinlocks = true; //ThreadPool::availableCores() > 6;
    }
    
    ~MergeSortCommandImplementation(void) { }
    
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
    
    //this has been copied from bamtools utilities, since it isn't in the API
    bool ParseRegionString(const string& regionString,
                           const BamMultiReader& reader,
                           BamRegion& region);
    
    // data members
private:
    BamAlignment * GetMergedAlignment();
    void PostMergedAlignment(const BamAlignment & al);
    void PostMergedAlignment(BamAlignment * al);
    
    string m_tempFilenameStub;
    int m_numberOfRuns;
    int64_t m_numberOfAlignments;
    SamHeader m_header;
    RefVector m_references;
    vector<string> m_tempFilenames;
    bool sort_retval, merge_retval;
    MergeSortCommand & command;
    ThreadPool * thread_pool;
    ThreadPool * sort_thread_pool;
    SynchronizedQueue <BamAlignment *> merged_alignment_queue;
    sem_t * alignment_available; 
    bool use_spinlocks;
    
    //options:
    bool sort_by_names;
    bool mark_duplicates;
    bool remove_duplicates;
    bool compresstempfiles;
    int compression_level;
    size_t alignments_per_tempfile;
    
    // for performance reasons, we may use these for lockless coordination of access to merged_alignment_queue
    volatile int64_t alignment_available_spinlock, alignment_consumed_spinlock; 
public:
    class TempFileWriteJob : public ThreadJob
    {
    public:
        TempFileWriteJob(MergeSortCommand::MergeSortCommandImplementation * tool, vector<BamAlignment *> * buffer, string filename);
        virtual void runJob();
    protected:
        string filename;
        vector<BamAlignment *> * buffer;
        MergeSortCommand::MergeSortCommandImplementation * tool;
    };
    
    template<class T>
    class SortJob : public ThreadJob
    {
    public:
        SortJob(typename vector<T>::iterator begin,typename vector<T>::iterator end, pthread_mutex_t & completion_lock, const MergeSortCommandImplementation & implementation);
        virtual void runJob();
    protected:
        typename vector<T>::iterator begin;
        typename vector<T>::iterator end;
        pthread_mutex_t & completion_lock;
        const MergeSortCommandImplementation & m_implementation;
    };
};

void * MergeSortCommand::MergeSortCommandImplementation::RunMergeThread(void * data)
{
#ifdef __linux__
    prctl(PR_SET_NAME,"bt_merge",0,0,0);
#endif
    MergeSortCommand::MergeSortCommandImplementation * tool = (MergeSortCommand::MergeSortCommandImplementation *)data;
    tool->merge_retval = tool->RunMerge();
    
    return NULL;
}

void * MergeSortCommand::MergeSortCommandImplementation::RunSortThread(void * data)
{
#ifdef __linux__
    prctl(PR_SET_NAME,"bt_sort",0,0,0);
#endif
    
    MergeSortCommand::MergeSortCommandImplementation * tool = (MergeSortCommand::MergeSortCommandImplementation *)data;
    tool->sort_retval = tool->RunSort();
    
    return NULL;
}

const char * ALIGNMENT_QUEUE_NAME = "mergesort_align_queue";

//actual run creates threads for other 'run'commands
bool MergeSortCommand::MergeSortCommandImplementation::Run(void)
{
    timeval start_time;
    gettimeofday(&start_time, NULL);
    
    // options
    sort_by_names =command.vm.count("byname") != 0;
    compression_level = command.vm["compression"].as<int>();
    compresstempfiles = command.vm.count("compresstempfiles") != 0;
    alignments_per_tempfile = command.vm["n"].as<int>();

    mark_duplicates = command.vm.count("markduplicates") != 0;
    remove_duplicates = command.vm.count("removeduplicates") != 0;
    if(remove_duplicates)
        mark_duplicates = true;

    if(!command.nothreads) {
        thread_pool = new ThreadPool();
        sort_thread_pool = new ThreadPool();
    } else if(command.verbose)
        cerr << "Thread pool use disabled." << endl;

    sem_unlink(ALIGNMENT_QUEUE_NAME);
    alignment_available = sem_open(ALIGNMENT_QUEUE_NAME, O_CREAT | O_EXCL, 0, 0);
    if(alignment_available == 0)
        perror("Failed creating alignment queue semaphore.");

    //create threads to do the work
    pthread_t merge_thread, sort_thread;
    pthread_create(&merge_thread, NULL, MergeSortCommand::MergeSortCommandImplementation::RunMergeThread, this);
    pthread_create(&sort_thread, NULL, MergeSortCommand::MergeSortCommandImplementation::RunSortThread, this);
    
    //wait for both threads to finish
    pthread_join(merge_thread, NULL);
    pthread_join(sort_thread, NULL);
    
    if(0 != sem_close(alignment_available))
        perror("Error closing alignment_available semaphore");
    
    if(0 != sem_unlink(ALIGNMENT_QUEUE_NAME))
        perror("Error unlinking alignment queue");
    
    if(!command.nothreads) {
        delete thread_pool;
        delete sort_thread_pool;
    }
    
    if(command.verbose)
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
bool MergeSortCommand::MergeSortCommandImplementation::RunMerge(void) {
    if(command.verbose)
        cerr << "Merging " << command.input_filenames.size() << " input files." << endl;
    
    // opens the BAM files (by default without checking for indexes)
    BamMultiReader reader;
    if ( !reader.Open(command.input_filenames) ) {
        cerr << "bamtools merge ERROR: could not open input BAM file(s)... Aborting." << endl;
        return false;
    }
    
    //save header information for temp files and final merge
    m_header = reader.GetHeader();
    m_references = reader.GetReferenceData();
    
    //the first semaphore post/wait signals header information is valid
    if(0 != sem_post(alignment_available))
        perror("Error posting alignment_available (open)");
    
    bool hasregion = command.vm.count("region") != 0;
    
    // if no region specified, store entire contents of file(s)
    if ( !hasregion ) {
        BamAlignment * al;
        while ( NULL != (al = reader.GetNextAlignmentCore()) ) {
            PostMergedAlignment(al);
        }
    }
    
    // otherwise attempt to use region as constraint
    else {
        string region_string = command.vm["region"].as<string>();
        // if region string parses OK
        BamRegion region;
        if ( ParseRegionString(region_string, reader, region) ) {
            
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
            cerr << "bamtools merge ERROR: could not parse REGION - " << region_string << endl;
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

void MergeSortCommand::MergeSortCommandImplementation::PostMergedAlignment( BamAlignment * al)
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

void MergeSortCommand::MergeSortCommandImplementation::PostMergedAlignment(const BamAlignment & al)
{
    merged_alignment_queue.push(new BamAlignment(al));
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

BamAlignment * MergeSortCommand::MergeSortCommandImplementation::GetMergedAlignment()
{
    if(use_spinlocks) {
        //give up our timeslice if the queue is empty
        while(alignment_available_spinlock == alignment_consumed_spinlock) usleep(2000);
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
bool MergeSortCommand::MergeSortCommandImplementation::GenerateSortedRuns(void) {
    
    if(command.verbose)
        cerr << "Generating sorted temp files...";
    
    //wait for header information
    if(0 != sem_wait(alignment_available))
        perror("Error waiting for alignment_available semaphore");
    
    // get basic data that will be shared by all temp/output files
    m_header.SortOrder = ( sort_by_names
                          ? Constants::SAM_HD_SORTORDER_QUERYNAME
                          : Constants::SAM_HD_SORTORDER_COORDINATE );
    
    // set up alignments buffer
    vector<BamAlignment *> * buffer = new vector<BamAlignment *>();
    buffer->reserve( (size_t)(alignments_per_tempfile*1.1) );
    bool bufferFull = false;
    BamAlignment * al = NULL;
    
    // if sorting by name, we need to generate full char data
    // so can't use GetNextAlignmentCore()
    if ( sort_by_names ) {
        
        // iterate through file
        while (true) {
            al = GetMergedAlignment();
            if(!al)
                break;
            
            // check buffer's usage
            bufferFull = ( buffer->size() >= alignments_per_tempfile );
            
            // store alignments until buffer is "full"
            if ( !bufferFull )
                buffer->push_back(al);
            
            // if buffer is "full"
            else {
                // so create a sorted temp file with current buffer contents
                // then push "al" into fresh buffer
                CreateSortedTempFile(buffer);
                buffer = new vector<BamAlignment *>();
                buffer->reserve( (size_t)(alignments_per_tempfile*1.1) );
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
            bufferFull = ( buffer->size() >= alignments_per_tempfile );
            
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
    
    if(command.verbose)
        cerr << "waiting for files to be compressed / written...";
    
    //wait for all temp files to be created in other threads
    if(!command.nothreads)
        thread_pool->waitForJobCompletion();
    
    if(command.verbose)
        cerr << "done." << endl;
    
    return true;
}


// this has been copied from bamtools utilities, since it isn't in the API. Original file is bamtools_utilities.cpp
bool MergeSortCommand::MergeSortCommandImplementation::ParseRegionString(const string& regionString,
                                  const BamMultiReader& reader,
                                  BamRegion& region)
{
    // -------------------------------
    // parse region string
    
    // check first for empty string
    if ( regionString.empty() ) 
        return false;   
    
    // non-empty string, look for a colom
    size_t foundFirstColon = regionString.find(':');
    
    // store chrom strings, and numeric positions
    string startChrom;
    string stopChrom;
    int startPos;
    int stopPos;
    
    // no colon found
    // going to use entire contents of requested chromosome 
    // just store entire region string as startChrom name
    // use BamReader methods to check if its valid for current BAM file
    if ( foundFirstColon == string::npos ) {
        startChrom = regionString;
        startPos   = 0;
        stopChrom  = regionString;
        stopPos    = -1;
    }
    
    // colon found, so we at least have some sort of startPos requested
    else {
        
        // store start chrom from beginning to first colon
        startChrom = regionString.substr(0,foundFirstColon);
        
        // look for ".." after the colon
        size_t foundRangeDots = regionString.find("..", foundFirstColon+1);
        
        // no dots found
        // so we have a startPos but no range
        // store contents before colon as startChrom, after as startPos
        if ( foundRangeDots == string::npos ) {
            startPos   = atoi( regionString.substr(foundFirstColon+1).c_str() ); 
            stopChrom  = startChrom;
            stopPos    = -1;
        } 
        
        // ".." found, so we have some sort of range selected
        else {
            
            // store startPos between first colon and range dots ".."
            startPos = atoi( regionString.substr(foundFirstColon+1, foundRangeDots-foundFirstColon-1).c_str() );
            
            // look for second colon
            size_t foundSecondColon = regionString.find(':', foundRangeDots+1);
            
            // no second colon found
            // so we have a "standard" chrom:start..stop input format (on single chrom)
            if ( foundSecondColon == string::npos ) {
                stopChrom  = startChrom;
                stopPos    = atoi( regionString.substr(foundRangeDots+2).c_str() );
            }
            
            // second colon found
            // so we have a range requested across 2 chrom's
            else {
                stopChrom  = regionString.substr(foundRangeDots+2, foundSecondColon-(foundRangeDots+2));
                stopPos    = atoi( regionString.substr(foundSecondColon+1).c_str() );
            }
        }
    }
    
    // -------------------------------
    // validate reference IDs & genomic positions
    
    const RefVector references = reader.GetReferenceData();
    
    // if startRefID not found, return false
    int startRefID = reader.GetReferenceID(startChrom);
    if ( startRefID == -1 ) return false;
    
    // startPos cannot be greater than or equal to reference length
    const RefData& startReference = references.at(startRefID);
    if ( startPos >= startReference.RefLength ) return false;
    
    // if stopRefID not found, return false
    int stopRefID = reader.GetReferenceID(stopChrom);
    if ( stopRefID == -1 ) return false;
    
    // stopPosition cannot be larger than reference length
    const RefData& stopReference = references.at(stopRefID);
    if ( stopPos > stopReference.RefLength ) return false;
    
    // if no stopPosition specified, set to reference end
    if ( stopPos == -1 ) stopPos = stopReference.RefLength;
    
    // -------------------------------
    // set up Region struct & return
    
    region.LeftRefID     = startRefID;
    region.LeftPosition  = startPos;
    region.RightRefID    = stopRefID;;
    region.RightPosition = stopPos;
    return true;
}

MergeSortCommand::MergeSortCommandImplementation::TempFileWriteJob::TempFileWriteJob(MergeSortCommand::MergeSortCommandImplementation * tool, vector<BamAlignment *> * buffer, string filename) :
filename(filename), buffer(buffer), tool(tool)
{
}

void MergeSortCommand::MergeSortCommandImplementation::TempFileWriteJob::runJob()
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

bool MergeSortCommand::MergeSortCommandImplementation::CreateSortedTempFile(vector<BamAlignment* > * buffer) {
    //make filename
    stringstream tempStr;
    tempStr << m_tempFilenameStub << m_numberOfRuns;
    m_tempFilenames.push_back(tempStr.str());
    
    ++m_numberOfRuns;
    
    bool success = true;
    
    if(command.nothreads) {
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
bool MergeSortCommand::MergeSortCommandImplementation::MergeSortedRuns(void) {
    string filename_out = command.vm["out"].as<string>();

    cerr << "Combining temp files for final output...";
    
    // open up multi reader for all of our temp files
    BamMultiReader multiReader;
    if ( !multiReader.Open(m_tempFilenames) ) {
        cerr << "mergesort ERROR: could not open BamMultiReader for merging temp files... Aborting."
        << endl;
        return false;
    }
    
    AlgorithmModuleAdaptor adaptor(multiReader.GetHeader(), multiReader.GetReferenceData());
    MarkDuplicates md;
    FileWriter writer;
    
    md.removeDuplicates = remove_duplicates;
    md.verbose = command.verbose;
    writer.setFilename(filename_out);
    writer.setCompressionLevel(compression_level);
    
    if(mark_duplicates)
    {
        adaptor.addSink(&md);
        md.addSink(&writer);

        adaptor.startAsync();
        md.startAsync();
    } else {
        adaptor.addSink(&writer);

        adaptor.startAsync();
    }
    writer.startAsync();

    // while data available in temp files
    BamAlignment * al;
    int count = 0;
    while (NULL != (al = multiReader.GetNextAlignment()) ) {
        adaptor.putAlignment(al);
        if(count++ % 1000000 == 0)
            cerr << ".";  //progress indicator every 1M alignments written.
    }

    // close files
    multiReader.Close();
    
    adaptor.done();
    
    writer.finishAsync();
    if(mark_duplicates)
        md.finishAsync();
    adaptor.finishAsync();
    
    cerr << "done." << endl << "Clearing " << m_tempFilenames.size() << " temp files...";
    
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

bool MergeSortCommand::MergeSortCommandImplementation::RunSort(void) {
    // this does a single pass, chunking up the input file into smaller sorted temp files,
    // then write out using BamMultiReader to handle merging
    
    if ( GenerateSortedRuns() )
        return MergeSortedRuns();
    else
        return false;
}

template<class T> MergeSortCommand::MergeSortCommandImplementation::SortJob<T>::SortJob(typename vector<T>::iterator begin, typename vector<T>::iterator end, pthread_mutex_t & completion_lock, const MergeSortCommandImplementation & implementation)
: begin(begin)
, end(end)
, completion_lock(completion_lock)
, m_implementation(implementation)
{
}

template<class T> void MergeSortCommand::MergeSortCommandImplementation::SortJob<T>::runJob()
{
#ifdef __linux__
    prctl(PR_SET_NAME,"bt_tempfile_sort",0,0,0);
#endif
    if(m_implementation.sort_by_names)
        std::stable_sort( begin, end, Sort::ByName() );
    else
        std::stable_sort( begin, end, Sort::ByPosition() );
    
    pthread_mutex_unlock(&completion_lock);
}

//this function is designed to accept either BamAlignment or BamAlignment* as T
template<class T>
void MergeSortCommand::MergeSortCommandImplementation::SortBuffer(vector<T>& buffer) {
    
    if(command.nothreads)
    {
        if(sort_by_names)
            std::stable_sort( buffer.begin(), buffer.end(), Sort::ByName() );
        else
            std::stable_sort( buffer.begin(), buffer.end(), Sort::ByPosition() );
    } else {
        int divisions = buffer.size() / MERGESORT_MIN_SORT_SIZE;
        if(divisions > ThreadPool::availableCores())
            divisions = ThreadPool::availableCores();
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

            jobs[ctr] = new SortJob<T>(begin, end, locks[ctr], *this);
            sort_thread_pool->addJob(jobs[ctr]);
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
            
            if(sort_by_names)
                std::inplace_merge(buffer.begin(), midpoint, end, Sort::ByName() );
            else
                std::inplace_merge(buffer.begin(), midpoint, end, Sort::ByPosition() );
        }
        
        delete [] jobs;
        delete [] locks;
    }
}

bool MergeSortCommand::MergeSortCommandImplementation::WriteTempFile(const vector<BamAlignment>& buffer,
                                                        const string& tempFilename)
{
    // open temp file for writing
    BamWriter tempWriter;
    
    if(compresstempfiles)
        tempWriter.SetCompressionMode(BamWriter::Compressed);
    else
        tempWriter.SetCompressionMode(BamWriter::Uncompressed);
    
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

bool MergeSortCommand::MergeSortCommandImplementation::WriteTempFile(const vector<BamAlignment *>& buffer,
                                                        const string& tempFilename)
{
    // open temp file for writing
    BamWriter tempWriter;
    
    if(compresstempfiles)
        tempWriter.SetCompressionMode(BamWriter::Compressed);
    else
        tempWriter.SetCompressionMode(BamWriter::Uncompressed);
    
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

////////////////////////////
// Command interface

void MergeSortCommand::getOptions()
{
    
    options.add_options()
    ("out,o", po::value<string>()->default_value("stdout"), "Output filename. Omit for stdout.")
    ("compression,c", po::value<int>()->default_value(6), "Compression level of the output. Valid 0-9.")
    ("region,r", po::value<string>(), "Genomic region to use.")
    ("byname,b", "Sort by name. Otherwise, sorts by position.")
    ("n,n", po::value<int>()->default_value(5e5), "Alignments per temp file.")
    ("compresstempfiles,C", "Compress temp files. By default, uncompressed")
    ("markduplicates,d", "Mark duplicates after sorting.")
    ("removeduplicates,R", "Remove duplicates.")
    ;
}

int MergeSortCommand::runCommand()
{
    MergeSortCommandImplementation impl(*this);
    
    // run MergeSort, return success/fail
    if ( impl.Run() )
        return 0;
    else
        return 1;
}

