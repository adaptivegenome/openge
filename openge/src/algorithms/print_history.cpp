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
#include <api/SamHeader.h>
using namespace std;
using BamTools::SamHeader;
using BamTools::SamProgramConstIterator;

int PrintHistory::runInternal()
{
    SamHeader h = getHeader();
    
    for( SamProgramConstIterator i = h.Programs.Begin(); i < h.Programs.End(); i++ ) {
        if( i->HasCommandLine())
            cout << i->CommandLine << endl;
        else
            cout << i->ID << endl;
    }
    
    while(true) {
        OGERead * a = getInputAlignment();
        if(!a) break;
        putOutputAlignment(a);
    }

    return 0;
}