/*********************************************************************
 *
 * commands.cpp: Base class for an OpenGE command.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 21 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include "commands.h"

#include <vector>
#include <string>
#include <iostream>
#include <boost/exception/all.hpp>

#include <sys/time.h>
#include <sys/resource.h>

#include "../util/thread_pool.h"

#include "../algorithms/algorithm_module.h"

using namespace std;
namespace po = boost::program_options;

int OpenGECommand::runWithParameters(int argc, const char ** argv)
{
    command_line.append("openge ");
    for(int i = 0; i < argc; i++) {
        command_line.append(argv[i]);
        command_line.append(" ");
    }
    
    options.add(global_options);
    options.add(io_options);
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
    tmpdir = vm["tmpdir"].as<string>() + "/";

    if(vm.count("in") == 0)
        input_filenames.push_back("stdin");
    else
        input_filenames = vm["in"].as<vector<string> >();
    
    num_threads = vm["threads"].as<unsigned int>();
    OGEParallelismSettings::setNumberThreads(num_threads);
    
    nothreads = vm.count("nothreads") != 0;
    
    AlgorithmModule::setNothreads(nothreads);
    AlgorithmModule::setVerbose(verbose);
    
    if(nothreads) {
        if(verbose)
            cerr << "Multithreading disabled." << endl;

        OGEParallelismSettings::disableMultithreading();
    } else {
        if(verbose)
            cerr << num_threads << " cores for use in thread pool." << endl;
        OGEParallelismSettings::enableMultithreading();
    }
    
    timeval start_time;
    gettimeofday(&start_time, NULL);
    
    int ret = runCommand();
    
    if(verbose)
    {
        //print cpu usage and time
        rusage r;
        getrusage(RUSAGE_SELF, &r);
        timeval stop_time;
        gettimeofday(&stop_time, NULL);
        timeval real_time;
        real_time.tv_sec = stop_time.tv_sec - start_time.tv_sec;
        real_time.tv_usec = stop_time.tv_usec - start_time.tv_usec;
        fprintf(stderr, "Elapsed time: %3ldm%06.3fs\n", real_time.tv_sec /60, float(real_time.tv_sec %60) + (1.e-6 * real_time.tv_usec) );
        fprintf(stderr, "CPU time: %3ldm%06.3fs (user) / %3ldm%06.3fs (sys)\n", r.ru_utime.tv_sec /60, float(r.ru_utime.tv_sec %60) + (1.e-6 * r.ru_utime.tv_usec), r.ru_stime.tv_sec /60, float(r.ru_stime.tv_sec %60) + (1.e-6 * r.ru_stime.tv_usec));
#ifndef __linux__
        // on linux, the max ram consumption is in KB. On Mac, it is in bytes.
        r.ru_maxrss /= 1024;
#endif
        fprintf(stderr, "Max mem: %6ld MB\n", r.ru_maxrss /1024);
    }
    
    BamTools::BamAlignment::clearCachedAllocations();
    
    return ret;
}

OpenGECommand::OpenGECommand()
{
    io_options.add_options()
    ("in,i", po::value<vector<string> >(),"Input files. If not specified, defaults to stdin. Can be specified without --in or -i")
    ("format,F", po::value<string>(),"File output format")
    ("compression,c", po::value<int>()->default_value(6), "Compression level of the output. Valid 0-9.")
    ;
    global_options.add_options()
    ("verbose,v" ,"Display detailed messages while processing")
    ("threads,t", po::value<unsigned int>()->default_value(ThreadPool::availableCores()), "Select the number of threads to be used in each threadpool")
    ("nothreads,d", "Disable use of thread pools for parallel processing.")
    ("tmpdir,T", po::value<string>()->default_value("/tmp"), "Directory to use for temporary files")
    ("nosplit","Do not split by chromosome (for speed) when processing")
    ;
}

OpenGECommand * CommandMarshall::commandWithName(const string name) {
    const char * cname = name.c_str();
    
    if(!strcmp(cname, "compare"))
        return new CompareCommand;
    if(!strcmp(cname, "count"))
        return new CountCommand;
    if(!strcmp(cname, "coverage"))
        return new CoverageCommand;
    else if(!strcmp(cname, "dedup"))
        return new DedupCommand;
    else if(!strcmp(cname, "help"))
        return new HelpCommand;
    else if(!strcmp(cname, "history"))
        return new HistoryCommand;
    else if(!strcmp(cname, "localrealign"))
        return new LocalRealignCommand;
    else if(!strcmp(cname, "mergesort"))
        return new MergeSortCommand;
    else if(!strcmp(cname, "stats"))
        return new StatsCommand;
    else if(!strcmp(cname, "version"))
        return new VersionCommand;
    else if(!strcmp(cname, "view"))
        return new ViewCommand;
    else {
        cerr << "Unknown command " << name << "." << endl;
        return NULL;
    }
}