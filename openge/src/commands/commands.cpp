#include "commands.h"

#include <vector>
#include <string>
#include <iostream>
#include <boost/exception/all.hpp>
#include <api/BamParallelismSettings.h>

#include "../util/thread_pool.h"

using namespace std;
namespace po = boost::program_options;

int OpenGECommand::runWithParameters(int argc, const char ** argv)
{
    options_positional.add("in", -1);
    getOptions();
    
    try {
        po::store(po::command_line_parser(argc, argv).options(options).positional(options_positional).run(), vm) ;
        po::notify(vm);    
    } catch( boost::exception & e )
    {
        cerr << "Unrecognized parameters." << endl << endl << "Valid options are:" << endl << options << endl;
        return -1;
    }
    
    verbose = 0 < vm.count("verbose");

    if(vm.count("in") == 0)
        input_filenames.push_back("stdin");
    else
        input_filenames = vm["in"].as<vector<string> >();
    
    num_threads = vm["threads"].as<unsigned int>();
    BamParallelismSettings::setNumberThreads(num_threads);
    
    nothreads = vm.count("nothreads") != 0;
    
    if(nothreads) {
        if(verbose)
            cerr << "Multithreading disabled." << endl;

        BamParallelismSettings::disableMultithreading();
    } else {
        if(verbose)
            cerr << num_threads << " cores for use in thread pool." << endl;
        BamParallelismSettings::enableMultithreading();
    }
    
    return runCommand();
}

OpenGECommand::OpenGECommand()
{
    options.add_options()
    ("in,i", po::value<vector<string> >(),"Input files. If not specified, defaults to stdin. Can be specified without --in or -i")
    ("verbose,v" ,"Display detailed messages while processing")
    ("threads,t", po::value<unsigned int>()->default_value(ThreadPool::availableCores()), "Select the number of threads to be used in each threadpool")
    ("nothreads,d", "Disable use of thread pools for parallel processing.")
    ;
}

OpenGECommand * CommandMarshall::commandWithName(const string name) {
    const char * cname = name.c_str();
    
    if(!strcmp(cname, "convert"))
        return new ConvertCommand;
    else if(!strcmp(cname, "count"))
        return new CountCommand;
    else if(!strcmp(cname, "dedup"))
        return new DedupCommand;
    else if(!strcmp(cname, "head"))
        return new HeadCommand;
    else if(!strcmp(cname, "help"))
        return new HelpCommand;
    else if(!strcmp(cname, "mergesort"))
        return new MergeSortCommand;
    else if(!strcmp(cname, "stats"))
        return new StatsCommand;
    else if(!strcmp(cname, "view"))
        return new ViewCommand;
    else {
        cerr << "Unknown command " << name << "." << endl;
        return NULL;
    }
}