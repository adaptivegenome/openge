/*********************************************************************
 *
 * command_localrealignment.cpp: Perofrm local realignment.
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
using namespace std;
namespace po = boost::program_options;

#include "../algorithms/print_history.h"
#include "../algorithms/file_reader.h"

void HistoryCommand::getOptions()
{}

int HistoryCommand::runCommand()
{
    FileReader reader;
    PrintHistory history;
    
    reader.addFiles(input_filenames);
    
    reader.addSink(&history);
    
    reader.runChain();

    return 0;
}