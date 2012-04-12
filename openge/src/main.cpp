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
