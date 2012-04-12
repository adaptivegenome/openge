#include "commands.h"

#include <vector>
#include <string>
using namespace std;

#include <api/BamMultiReader.h>
using namespace BamTools;
namespace po = boost::program_options;

void CountCommand::getOptions()
{ }

int CountCommand::runCommand()
{
    BamMultiReader reader;
    
    if(!reader.Open(input_filenames)) {
        cerr << "Error opening bam files." << endl;
        return -1;
    }
    
    size_t count = 0;
    BamAlignment al;
    while(reader.GetNextAlignment(al))
        count++;
    
    reader.Close();
    
    cout << count << endl;
    
    return count;
}