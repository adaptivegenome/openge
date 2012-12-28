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

#include "../util/bamtools/Sort.h"

#include "mark_duplicates.h"

#include "../util/bam_serializer.h"
#include "../util/bgzf_output_stream.h"

#include <pthread.h>

using namespace BamTools::Algorithms;

#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include <iostream>
using namespace std;

// generates mutiple sorted temp BAM files from single unsorted BAM file
bool ReadSorter::GenerateSortedRuns(void) {
    
    if(isVerbose())
        cerr << "Generating sorted temp files." << endl;
    
    // set up alignments buffer
    vector<OGERead *> * buffer = new vector<OGERead *>();
    buffer->reserve( (size_t)(alignments_per_tempfile*1.1) );

    // iterate through file
    while (true) {
        OGERead * al = getInputAlignment();
        if(!al)
            break;

        // store alignments until buffer is "full"
        if ( buffer->size() < alignments_per_tempfile )
            buffer->push_back(al);
        
        // if buffer is "full"
        else {
            // create a sorted temp file with current buffer contents
            // then push "al" into fresh buffer
            CreateSortedTempFile(buffer);
            buffer = new vector<OGERead *>();
            buffer->reserve( (size_t)(alignments_per_tempfile*1.1) );
            buffer->push_back(al);
        }
        if(read_count % 100000 == 0 && verbose)
            cerr << "\rRead " << read_count/1000 << "K reads." << flush;
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
    // do sorting
    tool->SortBuffer(*buffer);

    // as noted in the comments of the original file, success is never
    // used so we never return it.
    bool success = tool->WriteTempFile( *buffer, filename );
    if(!success) {
        cerr << "Problem writing out sorted temporary file. Aborting." << filename;
        exit(-1);
    }

    for(size_t i = 0; i < buffer->size(); i++)
        OGERead::deallocate((*buffer)[i]);
    delete buffer;
}

bool ReadSorter::CreateSortedTempFile(vector<OGERead* > * buffer) {
    //make filename
    stringstream filename_ss;
    filename_ss << m_tempFilenameStub << "_" << m_numberOfRuns << ".bam";
    string filename = filename_ss.str();
    m_tempFilenames.push_back(filename);

    ++m_numberOfRuns;
    
    bool success = true;
    
    if(isNothreads()) {
        // do sorting
        SortBuffer(*buffer);
        
        // write sorted contents to temp file, store success/fail
        success = WriteTempFile( *buffer, filename );
        for(size_t i = 0; i < buffer->size(); i++)
            OGERead::deallocate((*buffer)[i]);
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
    for (vector<string>::const_iterator tempIter = m_tempFilenames.begin() ; tempIter != m_tempFilenames.end(); ++tempIter ) {
        remove(tempIter->c_str());
    }

    if(isVerbose())
        cerr << "done." << endl;

    // return success
    return true;
}

//this function is designed to accept either OGERead or OGERead* as T
template<class T>
void ReadSorter::SortBuffer(vector<T>& buffer) {
    if(isNothreads())
    {
        if(sort_order == BamHeader::SORT_COORDINATE)
            std::stable_sort( buffer.begin(), buffer.end(), Sort::ByName() );
        else
            std::stable_sort( buffer.begin(), buffer.end(), Sort::ByPosition() );
    } else {
        if(sort_order == BamHeader::SORT_QUERYNAME)
            ogeSortMt( buffer.begin(), buffer.end(), Sort::ByName() );
        else
            ogeSortMt( buffer.begin(), buffer.end(), Sort::ByPosition() );
    }
}

bool ReadSorter::WriteTempFile(const vector<OGERead *>& buffer, const string& tempFilename)
{
    // open temp file for writing
    BamSerializer<BgzfOutputStream> tempWriter;
    
    tempWriter.getOutputStream().setCompressionLevel(compress_temp_files ? 6 : 0);
    
    m_header_access.lock();
    if ( !tempWriter.open(tempFilename, m_header) ) {
        cerr << "ReadSorter ERROR: could not open tempfile " << tempFilename
        << " for writing." << endl;
        exit(-1);
    }
    m_header_access.unlock();
    
    // write data
    for (vector<OGERead *>::const_iterator buffIter = buffer.begin() ; buffIter != buffer.end(); ++buffIter )  {
        tempWriter.write(**buffIter);
    }
    
    // close temp file & return success
    tempWriter.close();
    return true;
}

const BamHeader & ReadSorter::getHeader()
{
    while(true) {
        m_header_access.lock();
        bool ret = header_loaded;
        m_header_access.unlock();
        
        if(ret) break;
        usleep(10000);
    }

    return m_header;
}

int ReadSorter::runInternal()
{
    // options
    if(!isNothreads()) {
        thread_pool = new ThreadPool();
    } else if(isVerbose())
        cerr << "Thread pool use disabled." << endl;
    
    m_header_access.lock();
    m_header = AlgorithmModule::getHeader();
    m_header.setSortOrder( sort_order );
    
    header_loaded = true;
    m_header_access.unlock();
    
    bool retval = GenerateSortedRuns();
    
    if(retval)
        retval = MergeSortedRuns();
    
    if(!isNothreads()) {
        delete thread_pool;
    }
    
    return retval;
}
