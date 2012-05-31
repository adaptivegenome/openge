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

#include <api/SamHeader.h>
#include <api/BamAlignment.h>
#include <api/BamAux.h>

#include <set>

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
    void addSink(AlgorithmModule * sink) { sinks.insert(sink); sink->setSource(this);}
    bool removeSink(AlgorithmModule * sink) { sink->setSource(NULL); return sinks.erase(sink); }

    // Once the data flow between modules has been set up, call this function on any module in 
    // the chain to begin processing
    int runChain();

    // Start module in its own thread and return.
    int startAsync();
    
    // Wait for a module to complete execution, if it has been started by startAsync
    int finishAsync();
protected:

    // The following function is implemented by the superclass, and should be used to 
    // communicate data between algorithm modules. There is a single input queue in each module,
    // which is fed by modules supplying data to this module. There is no actual output queue-
    // the putOutputAlignment function calls putInputAlignment on any modules set up to receive
    // data from this module (using addSink,etc).
    //
    // IMPORTANT: If a module passes a BamAlignment to another module, the destination module is
    // responsible for deleting it!
    void putInputAlignment(BamTools::BamAlignment * read);
    
    // Use these functions when processing to get input data, and to pass the data to the next 
    // module in the chain.
    void putOutputAlignment(BamTools::BamAlignment * read);
    BamTools::BamAlignment * getInputAlignment();

public:
    virtual BamTools::SamHeader getHeader();
    virtual BamTools::RefVector getReferences();
    
    bool isVerbose() const { return AlgorithmModule::verbose; }
    bool isNothreads() const { return AlgorithmModule::nothreads; }
    static void setNothreads(bool nothreads);
    static void setVerbose(bool verbose);
    
protected:
    // All data processing should be performed in this function.
    virtual int runInternal() = 0;
    
    // Do not override this function.
    int runChildren();
    
    // Internal
    static void * algorithm_module_run(void * in);
    void setSource(AlgorithmModule * src) { source = src; }

    std::set<AlgorithmModule *> sinks;
    AlgorithmModule * source;
    SynchronizedQueue<BamTools::BamAlignment *> input_queue;
    pthread_t thread;
    bool finished_execution;
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