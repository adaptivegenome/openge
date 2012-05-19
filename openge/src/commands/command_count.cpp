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