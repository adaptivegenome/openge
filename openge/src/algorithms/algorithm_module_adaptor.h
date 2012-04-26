
#ifndef ALGORITHM_MODULE_ADAPTOR_H
#define ALGORITHM_MODULE_ADAPTOR_H

#include "algorithm_module.h"
#include "api/BamAlignment.h"

#include <pthread.h>

// AlgorithmModuleAdaptor
//
// This class acts as a source of reads when chaining together modules, and is useful
// when using/creating a source of reads other than a BAM file. Ideally, all algorithms
// and file readers should be written in subclasses of AlgorithmModule, but in the mean
// time, we use this adapter class.

class AlgorithmModuleAdaptor : public AlgorithmModule
{
    
protected:
    virtual int runInternal();
    const BamTools::SamHeader header;
    const BamTools::RefVector reference_data;
    virtual BamTools::SamHeader getHeader() { return header;}
    virtual BamTools::RefVector getReferences() { return reference_data; }
public:
    AlgorithmModuleAdaptor(const BamTools::SamHeader & header, const BamTools::RefVector & reference_data);
    void putAlignment(BamTools::BamAlignment * read) { putOutputAlignment(read); }
    void done();

    pthread_mutex_t mutex;
};

#endif