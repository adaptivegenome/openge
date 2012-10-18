/*********************************************************************
 *
 * read_sorter.cpp: Sort a stream of reads.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 20 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************
 *
 * This file is based on the algorithm in bamtools_sort.cpp, but has
 * been parallelized and refactored as part of OpenGE. Original authors
 * are Derek Barnett, Erik Garrison, Lee C. Baker
 * Marth Lab, Department of Biology, Boston College
 *********************************************************************/

#include "read_sorter.h"
#include "../util/thread_pool.h"
#include <cassert>
using namespace std;

#include <api/algorithms/Sort.h>

#include <api/SamConstants.h>

#include "mark_duplicates.h"

#include "../util/bam_serializer.h"
#include "../util/bgzf_output_stream.h"

#include <pthread.h>

using BamTools::SamHeader;
using namespace BamTools::Algorithms;

#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include <iostream>
using namespace std;

const unsigned int SORT_DEFAULT_MAX_BUFFER_COUNT  = 500000;  // max numberOfAlignments for buffer
const unsigned int SORT_DEFAULT_MAX_BUFFER_MEMORY = 1024;    // Mb
const unsigned int MERGESORT_MIN_SORT_SIZE = 30000;    //don't parallelize sort jobs smaller than this many alignments

//actual run creates threads for other 'run'commands
bool ReadSorter::Run(void)
{
    // options
    if(!isNothreads()) {
        thread_pool = new ThreadPool();
        sort_thread_pool = new ThreadPool();
    } else if(isVerbose())
        cerr << "Thread pool use disabled." << endl;
    
    m_header = getHeader();
    m_header.SortOrder = ( sort_order == SORT_NAME
                          ? BamTools::Constants::SAM_HD_SORTORDER_QUERYNAME
                          : BamTools::Constants::SAM_HD_SORTORDER_COORDINATE );

    RunSort();

    if(!isNothreads()) {
        delete thread_pool;
        delete sort_thread_pool;
    }
    
    return sort_retval && merge_retval;
}

// generates mutiple sorted temp BAM files from single unsorted BAM file
bool ReadSorter::GenerateSortedRuns(void) {
    
    if(isVerbose())
        cerr << "Generating sorted temp files." << endl;

    // get basic data that will be shared by all temp/output files
    m_header.SortOrder = ( sort_order == SORT_NAME
                          ? BamTools::Constants::SAM_HD_SORTORDER_QUERYNAME
                          : BamTools::Constants::SAM_HD_SORTORDER_COORDINATE );
    
    // set up alignments buffer
    vector<OGERead *> * buffer = new vector<OGERead *>();
    buffer->reserve( (size_t)(alignments_per_tempfile*1.1) );
    bool bufferFull = false;
    OGERead * al = NULL;
    
    // if sorting by name, we need to generate full char data
    // so can't use GetNextAlignmentCore()
    if ( sort_order == SORT_NAME ) {
        
        // iterate through file
        while (true) {
            al = getInputAlignment();
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
                buffer = new vector<OGERead *>();
                buffer->reserve( (size_t)(alignments_per_tempfile*1.1) );
                buffer->push_back(al);
            }
        }
    }
    
    // sorting by position, can take advantage of GNACore() speedup
    else {
        // iterate through file
        while (true) {
            al = getInputAlignment();
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
                buffer = new vector<OGERead *>();
                buffer->push_back(al);
            }
            if(read_count % 100000 == 0 && verbose)
                cerr << "\rRead " << read_count/1000 << "K reads." << flush;
        }
    }
    
    // handle any leftover buffer contents
    if ( !buffer->empty() )
        CreateSortedTempFile(buffer);

    //wait for all temp files to be created in other threads
    if(!isNothreads())
        thread_pool->waitForJobCompletion();
    
    if(isVerbose())
        cerr << "\rRead " << read_count/1000 << "K reads (done)." << endl;
    
    return true;
}

ReadSorter::TempFileWriteJob::TempFileWriteJob(ReadSorter * tool, vector<OGERead *> * buffer, string filename) :
filename(filename), buffer(buffer), tool(tool)
{
}

void ReadSorter::TempFileWriteJob::runJob()
{
    ogeNameThread("sort_tmp_sort");

    // do sorting
    tool->SortBuffer(*buffer);
    
    ogeNameThread("sort_tmp_write");

    // as noted in the comments of the original file, success is never
    // used so we never return it.
    bool success = tool->WriteTempFile( *buffer, filename );
    if(!success)
        cerr << "Problem writing out temporary file " << filename;
    

    ogeNameThread("sort_tmp_cleanup");

    for(size_t i = 0; i < buffer->size(); i++)
        delete (*buffer)[i];
    delete buffer;
}

bool ReadSorter::CreateSortedTempFile(vector<OGERead* > * buffer) {
    //make filename
    stringstream filename_ss;
    filename_ss << tmp_file_dir << m_tempFilenameStub << m_numberOfRuns;
    string filename = filename_ss.str();
    m_tempFilenames.push_back(filename);

    ++m_numberOfRuns;
    
    bool success = true;
    
    if(isNothreads()) {
        // do sorting
        SortBuffer(*buffer);
        
        // write sorted contents to temp file, store success/fail
        success = WriteTempFile( *buffer, filename );
        delete buffer;
        
    } else {
        TempFileWriteJob * job = new TempFileWriteJob(this, buffer, filename);
        thread_pool->addJob(job);
    }
    
    // return success/fail of writing to temp file
    // TODO: a failure returned here is not actually caught and handled anywhere
    return success;
}

// merges sorted temp BAM files into single sorted output BAM file
bool ReadSorter::MergeSortedRuns(void) {
    if(verbose)
        cerr << "Combining temp files for final output..." << endl;

    MultiReader readers;
    
    if(!readers.open(m_tempFilenames)) {
        cerr << "Error opening reader for tempfiles: " << endl;

        for (vector<string>::const_iterator tempIter = m_tempFilenames.begin() ; tempIter != m_tempFilenames.end(); ++tempIter )
            cerr << "   " << *tempIter << endl;
        exit(-1);
    }

    while(true) {
        OGERead * a = readers.read();
        
        if(!a)
            break;
        
        putOutputAlignment(a);

        if(write_count % 100000 == 0 && verbose)
            cerr << "\rCombined " << write_count/1000 << "K reads (" << 100 * write_count / read_count << "%)." << flush;
    }
    if(verbose && read_count)
        cerr << "\rCombined " << write_count/1000 << "K reads (" << 100 * write_count / read_count << "%)." << endl;

    if(verbose)
        cerr << "Clearing " << m_tempFilenames.size() << " temp files...";
    
    // delete all temp files
    vector<string>::const_iterator tempIter = m_tempFilenames.begin();
    vector<string>::const_iterator tempEnd  = m_tempFilenames.end();
    for ( ; tempIter != tempEnd; ++tempIter ) {
        remove(tempIter->c_str());
    }

    if(isVerbose())
        cerr << "done." << endl;

    // return success
    return true;
}

bool ReadSorter::RunSort(void) {
    // this does a single pass, chunking up the input file into smaller sorted temp files,
    // then merging in the results from multiple readers.
    
    if ( GenerateSortedRuns() )
        return MergeSortedRuns();
    else
        return false;
}

//this function is designed to accept either BamAlignment or BamAlignment* as T
template<class T>
void ReadSorter::SortBuffer(vector<T>& buffer) {
    if(isNothreads())
    {
        if(sort_order == SORT_NAME)
            std::stable_sort( buffer.begin(), buffer.end(), Sort::ByName() );
        else
            std::stable_sort( buffer.begin(), buffer.end(), Sort::ByPosition() );
    } else {
        if(sort_order == SORT_NAME)
            ogeSortMt( buffer.begin(), buffer.end(), Sort::ByName() );
        else
            ogeSortMt( buffer.begin(), buffer.end(), Sort::ByPosition() );
    }
}

bool ReadSorter::WriteTempFile(const vector<OGERead *>& buffer,
                                                                     const string& tempFilename)
{
    // open temp file for writing
    BamSerializer<BgzfOutputStream> tempWriter;
    
    if(compress_temp_files)
        tempWriter.getOutputStream().setCompressionLevel(6);
    else
        tempWriter.getOutputStream().setCompressionLevel(0);
    
    if ( !tempWriter.open(tempFilename, m_header) ) {
        cerr << "ReadSorter ERROR: could not open tempfile " << tempFilename
        << " for writing." << endl;
        exit(-1);
    }
    
    // write data
    vector<OGERead *>::const_iterator buffIter = buffer.begin();
    vector<OGERead *>::const_iterator buffEnd  = buffer.end();
    for ( ; buffIter != buffEnd; ++buffIter )  {
        const OGERead * al = (*buffIter);
        tempWriter.write(*al);
    }
    
    // close temp file & return success
    tempWriter.close();
    return true;
}

SamHeader ReadSorter::getHeader()
{
    SamHeader header = source->getHeader();

    header.SortOrder = ( sort_order == SORT_NAME
                          ? BamTools::Constants::SAM_HD_SORTORDER_QUERYNAME
                          : BamTools::Constants::SAM_HD_SORTORDER_COORDINATE );
    return header;
}

int ReadSorter::runInternal()
{
    // run MergeSort, return success/fail
    return Run();
}
