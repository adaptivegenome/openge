#include "algorithm_module.h"
#include <pthread.h>
using namespace BamTools;
using namespace std;

AlgorithmModule::AlgorithmModule()
: source(NULL)
, finished_execution(false)
{ }

void * AlgorithmModule::algorithm_module_run(void * in)
{
    AlgorithmModule * m = (AlgorithmModule *) in;

    m->run_return_value = m->runInternal();
    m->finished_execution = true;

    return 0;
}

int AlgorithmModule::runChildren()
{
    startAsync();
    
    for(set<AlgorithmModule *>::iterator i = sinks.begin(); i != sinks.end(); i++)
        (*i)->runChildren();
    
    finishAsync();
    return 0;
}

int AlgorithmModule::runChain()
{
    AlgorithmModule * root_module = this;
    while(root_module->source)
        root_module = root_module->source;

    root_module->runChildren();

    finished_execution = true;
    
    return 0;
}

int AlgorithmModule::startAsync()
{
    pthread_create(&thread, NULL, AlgorithmModule::algorithm_module_run, this);
    return 0;
}

int AlgorithmModule::finishAsync()
{
    pthread_join(thread, NULL);
    return run_return_value;
}

void AlgorithmModule::putInputAlignment(BamAlignment * read)
{
    while(input_queue.size() > 6000)
        usleep(10000);
    input_queue.push(read);
}

void AlgorithmModule::putOutputAlignment(BamAlignment * read)
{
    if(sinks.size() == 0)
        delete read;

    for(set<AlgorithmModule *>::iterator i = sinks.begin(); i != sinks.end(); i++) {
        if(sinks.begin() == i)
            (*i)->putInputAlignment(read);
        else
            (*i)->putInputAlignment(new BamAlignment(*read));
    }
}

BamAlignment * AlgorithmModule::getInputAlignment()
{
    while(input_queue.size() == 0) {
        if(source->finished_execution)
            return NULL;
        usleep(20000);
    }

    return input_queue.pop();
}

SamHeader AlgorithmModule::getHeader()
{
    return source->getHeader();
}

RefVector AlgorithmModule::getReferences()
{
    return source->getReferences();
}
