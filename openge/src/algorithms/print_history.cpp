/*********************************************************************
 *
 * print_history.cpp: Print the history of commands from a file's header.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 2 August 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial
 * Purpose License. A copy of this license has been provided in
 * the openge/ directory.
 *
 *********************************************************************/

#include "print_history.h"
using namespace std;

int PrintHistory::runInternal()
{
    BamHeader h = getHeader();
    
    for( BamProgramRecords::const_iterator i = h.getPrograms().begin(); i < h.getPrograms().end(); i++ ) {
        if( i->getCommandLine().empty())
            cout << i->id << endl;
        else
            cout << i->getCommandLine() << endl;
    }
    
    while(true) {
        OGERead * a = getInputAlignment();
        if(!a) break;
        putOutputAlignment(a);
    }

    return 0;
}