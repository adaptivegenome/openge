#include "algorithm_module_adaptor.h"

AlgorithmModuleAdaptor::AlgorithmModuleAdaptor(const BamTools::SamHeader & header, const BamTools::RefVector & reference_data) 
: header(header)
, reference_data(reference_data)
{
    if(0 != pthread_mutex_init(&mutex, NULL))
        perror("Error creating AlgorithmModuleAdaptor mutex");
    if(0 != pthread_mutex_lock(&mutex))
        perror("Error locking AlgorithmModuleAdaptor mutex in initialization");
}

int AlgorithmModuleAdaptor::runInternal()
{
    if(0 != pthread_mutex_lock(&mutex))
        perror("Error locking AlgorithmModuleAdaptor mutex");
    if(0 != pthread_mutex_unlock(&mutex))
        perror("Error unlocking AlgorithmModuleAdaptor mutex");
    if(0 != pthread_mutex_destroy(&mutex))
        perror("Error destroying AlgorithmModuleAdaptor mutex");

    return 0;
}

void AlgorithmModuleAdaptor::done()
{
    if(0 != pthread_mutex_unlock(&mutex))
        perror("Error unlocking AlgorithmModuleAdaptor mutex on finish");
}