/*********************************************************************
 *
 * command_count.cpp: Count reads in one or more files.
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
#include "../algorithms/file_reader.h"
namespace po = boost::program_options;

void CountCommand::getOptions()
{ }

int CountCommand::runCommand()
{
    FileReader reader;

    reader.setLoadStringData(false);

    reader.addFiles(input_filenames);
    reader.runChain();
    
    cout << reader.getCount() << endl;
    
    return 0;
}