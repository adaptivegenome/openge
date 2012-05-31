/*********************************************************************
*
* main.cpp: entry point to the OpenGE command line executable.
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

#include <iostream>

using namespace std;

#include "commands/commands.h"


int main(int argc, const char ** argv)
{

    if(argc == 1)
    {
        cerr << "Usage:" << endl << "    openge command [options]" << endl << endl << "Run 'openge help' for more details." << endl;
        return 0;
    }

    OpenGECommand * command = CommandMarshall::commandWithName(argv[1]);
    
    if(!command)
        return -1;

    return command->runWithParameters(argc-1, &argv[1]);
}
