#ifndef OGE_ALGO_SORT_H
#define OGE_ALGO_SORT_H
/*********************************************************************
 *
 * read_sorter.h: Sort a stream of reads.
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

#include "algorithm_module.h"
#include "../commands/commands.h"
#include "../util/read_stream_reader.h"

#include <vector>
#include <string>

typedef enum {
    SORT_NAME,
    SORT_POSITION
} sorting_t;

class ReadSorter : public AlgorithmModule
{
public:
    ReadSorter(const std::string & temp_directory)
    : m_tempFilenameStub("/oge_sort_")
    , m_numberOfRuns(0)
    , m_numberOfAlignments(0)
    , sort_order (SORT_POSITION)
    , compress_temp_files (false)
    , alignments_per_tempfile(200000)
    {
        char buffer[16];
        pid_t pid = getpid();
        sprintf(buffer, "%d", pid);
        m_tempFilenameStub = temp_directory + m_tempFilenameStub + buffer;
    }
    
    ~ReadSorter(void) { }
    
    // interface
protected:
    bool Run(void);
    bool RunSort(void);
    bool RunMerge(void);
    
    static void * RunMergeThread(void *);
    static void * RunSortThread(void *);
    
    //from SortTool:
    bool CreateSortedTempFile(std::vector<OGERead *> * buffer);
    bool GenerateSortedRuns(void);
    bool MergeSortedRuns(void);
    bool WriteTempFile(const std::vector<OGERead *> & buffer, const std::string& tempFilename);
    template<class T> void SortBuffer(std::vector<T>& buffer);

    // data members
private:
    std::string m_tempFilenameStub;
    int m_numberOfRuns;
    int64_t m_numberOfAlignments;
    BamTools::SamHeader m_header;
    BamTools::RefVector m_references;
    std::vector<std::string> m_tempFilenames;
    bool sort_retval, merge_retval;
    ThreadPool * thread_pool;
    
    //options:
    sorting_t sort_order;
    bool compress_temp_files;
    size_t alignments_per_tempfile;

public:
    class TempFileWriteJob : public ThreadJob
    {
    public:
        TempFileWriteJob(ReadSorter * tool, std::vector<OGERead *> * buffer, std::string filename);
        virtual void runJob();
    protected:
        std::string filename;
        std::vector<OGERead *> * buffer;
        ReadSorter * tool;
    };

protected:
    virtual int runInternal();
    virtual const BamTools::SamHeader & getHeader();
    
public:
    sorting_t getSortBy() { return sort_order; }
    void setSortBy(sorting_t sort_order) { this->sort_order = sort_order; }

    bool getCompressTempFiles() { return compress_temp_files; }
    void setCompressTempFiles(bool compress_temp_files) { this->compress_temp_files = compress_temp_files; }

    size_t getAlignmentsPerTempfile() { return alignments_per_tempfile; }
    void setAlignmentsPerTempfile(size_t alignments_per_tempfile) { this->alignments_per_tempfile = alignments_per_tempfile; }
};
#endif