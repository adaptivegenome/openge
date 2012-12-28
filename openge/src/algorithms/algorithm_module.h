/*********************************************************************
 *
 * algorithm_module.h:  Algorithm module default implementations
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
 *********************************************************************
 *
 * The AlgorithmModule abstract class provides an interface into which
 * a stream of reads can be passed.
 *
 *********************************************************************/

#ifndef OGE_ALGO_MODULE_H
#define OGE_ALGO_MODULE_H

#include "../util/bam_header.h"
#include "../util/oge_read.h"

#include <vector>

#include "../util/thread_pool.h"
#include "../commands/commands.h"

class AlgorithmModule
{
public:
    AlgorithmModule();
    virtual ~AlgorithmModule() {};

    // Data flow management
    // Use the following functions to connect several modules together. These functions
    // should not be called after run() has been called (on any algorithm module).
    void addSink(AlgorithmModule * sink) { sinks.push_back(sink); sink->setSource(this);}
    bool removeSink(AlgorithmModule * sink);

    // Once the data flow between modules has been set up, call this function on any module in 
    // the chain to begin processing
    int runChain();

    // Start module in its own thread and return.
    int startAsync();
    
    // Wait for a module to complete execution, if it has been started by startAsync
    int finishAsync();
public:

    // The following function is implemented by the superclass, and should be used to 
    // communicate data between algorithm modules. There is a single input queue in each module,
    // which is fed by modules supplying data to this module. There is no actual output queue-
    // the putOutputAlignment function calls putInputAlignment on any modules set up to receive
    // data from this module (using addSink,etc).
    //
    // IMPORTANT: If a module passes a BamAlignment to another module, the destination module is
    // responsible for deleting it! This function should not by called by most algorithm
    // modules; only call if you are doing something creative with which modules you pass data
    // to, eg. SplitByChromosome.
    virtual void putInputAlignment(OGERead * read);

protected:
    // Use these functions when processing to get input data, and to pass the data to the next 
    // module in the chain.
    virtual void putOutputAlignment(OGERead * read);
    OGERead * getInputAlignment();

public:
    virtual const BamHeader & getHeader();
    
    bool isVerbose() const { return AlgorithmModule::verbose; }
    bool isNothreads() const { return AlgorithmModule::nothreads; }
    static void setNothreads(bool nothreads);
    static void setVerbose(bool verbose);
    
    size_t getReadCount() { return read_count; }
    size_t getWriteCount() { return write_count; }
    
protected:
    // All data processing should be performed in this function.
    virtual int runInternal() = 0;
    
    // Do not override this function.
    int runChildren();
    int waitForChildrenCompletion();
    
    // Internal
    static void * algorithm_module_run(void * in);
    void setSource(AlgorithmModule * src) { source = src; }

    std::vector<AlgorithmModule *> sinks;
    AlgorithmModule * source;
    SynchronizedQueue<OGERead *> input_queue;
    pthread_t thread;
    SynchronizedFlag finished_execution;
    int run_return_value;
    
    static bool verbose;
    static bool nothreads;
    size_t read_count, write_count;
};

// The black hole module exists to terminate a chain of other modules.
// Unfortunately, delete'ing the current BamAlignment takes forever,
// so this class helps by running the delete in a separate thread.
class BlackHoleModule : public AlgorithmModule
{
    virtual int runInternal();
};

#endif