/*********************************************************************
 *
 * command_bpipe.cpp: Execute a BPIPE command script
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 12 September 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include "commands.h"

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include "../util/bpipe.h"
using namespace std;
namespace po = boost::program_options;


//environ is a list of all environment variables
#include <unistd.h>
#include <cstdlib>
extern char ** environ;

void BPipeCommand::getOptions()
{
    options.add_options()
    ("test,t", "Reads and checks a scripts' pipeline without actually running the commands.")
    ("print", "Print the commands that will be executed by the pipeline.")
    ("print_execution,x", "Print the execution structure of the pipeline.")
    ("define,p", po::value<vector<string> >(), "Define a variable var=value")
    ;
}

int BPipeCommand::runCommand()
{
    if(input_filenames.size() != 1 && input_filenames.size() != 2) {
        cerr << "One input execute script is required." << endl;
        exit(-1);
    }
    
    BPipe pipe;
    
    if(!pipe.load(input_filenames[0].c_str())) {
        cerr << "Error loading execute file " << input_filenames[0] << endl;
        exit(-1);
    }
    string input_filename;
    if(input_filenames.size() > 1)
        input_filename = input_filenames[1];
    
    //set all environment variables as bpipe variables    
    for(int i = 0; environ[i]; i++) {
        char * eq_loc = strchr(environ[i], '=');
        
        if(!*eq_loc)
            break;
        string var(environ[i], eq_loc);
        string value(eq_loc + 1);

        pipe.define(var, value);
    }

    //set define
    if(vm.count("define")) {
        const vector<string> defines = vm["define"].as<vector<string> >();
        for(int i = 0; i < defines.size(); i++) {
            size_t equals_position = defines[i].find("=");
            
            if(equals_position == string::npos) {
                cerr << "Error: variable definitions must be in \"-p var=value\" form. Exiting." << endl;
                exit(-1);
            }
            
            string var = string(&defines[i][0], equals_position);
            string value = &defines[i][equals_position + 1];
            
            pipe.define(var, value);
        }
    }

    if(!pipe.check(input_filename)) {
        cerr << "Parsing bpipe file " << input_filenames[0] << " failed." << endl;
        exit(-1);
    }
    
    if(!vm.count("test") && !vm.count("print") && !vm.count("print_execution")) {
        if(!pipe.execute()) {
            cerr << "Executing bpipe file " << input_filenames[0] << " failed." << endl;
            exit(-1);
        }
    }
    
    if(vm.count("print"))
        pipe.print();
    
    return 0;
}