/*********************************************************************
 *
 * command_help.cpp: Display help information.
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

#include <iostream>
using namespace std;

namespace po = boost::program_options;

int HelpCommand::runCommand()
{
    if(1 != vm.count("command")) {
        cerr << "Usage: openge help command" << endl;
        cerr << "Valid commands are: count coverage dedup help mergesort stats version view" << endl;
        return 0;
    }
    
    string command_name = vm["command"].as<string>();
    
    OpenGECommand * command = CommandMarshall::commandWithName(command_name);
    
    if(!command) {
        cerr << "Invalid command \"" << command << "\"." << endl;
        return 0;
    }
    
    command->getOptions();
    
    cerr << "Global options:" << endl;
    cerr << command->global_options << endl;
    cerr << "File I/O options:" << endl;
    cerr << command->io_options << endl;
    cerr << command_name << "-specific options:" << endl;
    cerr << command->options << endl;
    
    delete command;
    
    return 0;
}

void HelpCommand::getOptions()
{
    options_positional.add("command", 1);
    
    options.add_options()
    ("command", po::value<string>(), "Command to display help for");
}
