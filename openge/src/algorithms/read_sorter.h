#ifndef OGE_ALGO_SORT_H
#define OGE_ALGO_SORT_H

#include "algorithm_module.h"
#include "../commands/commands.h"

#include <vector>
#include <string>

typedef enum {
    SORT_NAME,
    SORT_POSITION
} sorting_t;

class ReadSorter : public AlgorithmModule
{
public:
    ReadSorter()
    : m_tempFilenameStub("bam_mergesort_temp_")
    , m_numberOfAlignments(0)
    { }
    
    ~ReadSorter(void) { }
    
    // interface
protected:
    bool Run(void);
    bool RunSort(void);
    bool RunMerge(void);
    
    static void * RunMergeThread(void *);
    static void * RunSortThread(void *);
    
    //from SortTool:
    bool CreateSortedTempFile(std::vector<BamTools::BamAlignment *> * buffer);
    bool GenerateSortedRuns(void);
    bool MergeSortedRuns(void);
    bool WriteTempFile(const std::vector<BamTools::BamAlignment> & buffer, const std::string& tempFilename);
    bool WriteTempFile(const std::vector<BamTools::BamAlignment *> & buffer, const std::string& tempFilename);
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
    ThreadPool * sort_thread_pool;
    
    //options:
    bool compresstempfiles;
    size_t alignments_per_tempfile;

public:
    class TempFileWriteJob : public ThreadJob
    {
    public:
        TempFileWriteJob(ReadSorter * tool, std::vector<BamTools::BamAlignment *> * buffer, std::string filename);
        virtual void runJob();
    protected:
        std::string filename;
        std::vector<BamTools::BamAlignment *> * buffer;
        ReadSorter * tool;
    };
    
    template<class T>
    class SortJob : public ThreadJob
    {
    public:
        SortJob(typename std::vector<T>::iterator begin,typename std::vector<T>::iterator end, pthread_mutex_t & completion_lock, const ReadSorter & implementation);
        virtual void runJob();
    protected:
        typename std::vector<T>::iterator begin;
        typename std::vector<T>::iterator end;
        pthread_mutex_t & completion_lock;
        const ReadSorter & m_implementation;
    };
protected:
    virtual int runInternal();
    
    sorting_t sort_order;
    bool compress_temp_files;
public:
    class MergeSortCommandImplementation;
    MergeSortCommandImplementation * impl;

    sorting_t getSortBy() { return sort_order; }
    void setSortBy(sorting_t sort_order) { this->sort_order = sort_order; }

    bool getCompressTempFiles() { return compress_temp_files; }
    void setCompressTempFiles(bool compress_temp_files) { this->compress_temp_files = compress_temp_files; }

    size_t getAlignmentsPerTempfile() { return alignments_per_tempfile; }
    void setAlignmentsPerTempfile(size_t alignments_per_tempfile) { this->alignments_per_tempfile = alignments_per_tempfile; }
};
#endif