/*********************************************************************
 *
 * command_stats.cpp: Display statistics of a file.
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

#include "../algorithms/statistics.h"
#include "../algorithms/file_reader.h"

namespace po = boost::program_options;

void StatsCommand::getOptions()
{
    options.add_options()
    ("inserts,I", "Show detailed insert statistics")
    ("lengths,L", "Show read length statistics")
    ;
}

int StatsCommand::runCommand()
{
    Statistics s;
    FileReader reader;
    
    reader.addFiles( input_filenames );
    
    s.showInsertSizeSummary(vm.count("inserts"));
    s.showReadLengthSummary(vm.count("lengths"));
    
    reader.addSink(&s);
    
    reader.runChain();
    return 0;
}
