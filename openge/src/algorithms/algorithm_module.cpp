/*********************************************************************
*
* algorithm_module.cpp:  Algorithm module default implementations
* Open Genomics Engine
*
* Author: Lee C. Baker, VBI
* Last modified: 31 May 2012
*
*********************************************************************
*
* This file is released under the Virginia Tech Non-Commercial 
* Purpose License. A copy of this license has been provided in 
* the openge/ directory.
*
*********************************************************************/

#include "algorithm_module.h"
#include <pthread.h>
using namespace BamTools;
using namespace std;

int BlackHoleModule::runInternal()
{
    ogeNameThread("am_BlackHole");

    while(true) {
        BamTools::BamAlignment * r = getInputAlignment();
        if(!r)
            return 0;
        delete r;
    }
}

bool AlgorithmModule::verbose = false;
bool AlgorithmModule::nothreads = false;

void AlgorithmModule::setNothreads(bool set_nothreads) {
    nothreads = set_nothreads;
}

void AlgorithmModule::setVerbose(bool set_verbose) {
    verbose = set_verbose;
}

AlgorithmModule::AlgorithmModule()
: source(NULL)
, finished_execution(false)
, read_count(0)
, write_count(0)
{ }

void * AlgorithmModule::algorithm_module_run(void * in)
{
    AlgorithmModule * m = (AlgorithmModule *) in;

    m->run_return_value = m->runInternal();
    m->finished_execution = true;

    return 0;
}

bool AlgorithmModule::removeSink(AlgorithmModule * sink) 
{ 
    sink->setSource(NULL); 
    if(find(sinks.begin(), sinks.end(), sink) == sinks.end()) 
        return false; 

    sinks.erase(std::remove(sinks.begin(), sinks.end(), sink), sinks.end()); 

    return true;
}

int AlgorithmModule::runChildren()
{
    startAsync();
    
    for(vector<AlgorithmModule *>::iterator i = sinks.begin(); i != sinks.end(); i++)
        (*i)->runChildren();

    return 0;
}

int AlgorithmModule::waitForChildrenCompletion()
{
    for(vector<AlgorithmModule *>::iterator i = sinks.begin(); i != sinks.end(); i++)
        (*i)->waitForChildrenCompletion();
    
    finishAsync();
    return 0;
}

int AlgorithmModule::runChain()
{
    // TODO LCB- when running a chain, we should parse the tree and add BlackHoleModules to each leaf node
    // for performance reasons.
    AlgorithmModule * first_leaf = this;
    
    while(first_leaf->sinks.size() > 0)
        first_leaf = *(first_leaf->sinks.begin());
    
    BlackHoleModule * bh = new BlackHoleModule;
    first_leaf->addSink(bh);
    
    AlgorithmModule * root_module = this;
    while(root_module->source)
        root_module = root_module->source;

    root_module->runChildren();
    root_module->waitForChildrenCompletion();

    finished_execution = true;
    
    
    delete bh;
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
    pthread_detach(thread);
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
    write_count++;
    
    if(sinks.size() == 0)
        delete read;

    for(vector<AlgorithmModule *>::iterator i = sinks.begin(); i != sinks.end(); i++) {
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
        usleep(10000);
    }
    
    read_count++;

    return input_queue.pop();
}

SamHeader AlgorithmModule::getHeader()
{
    return source->getHeader();
}
