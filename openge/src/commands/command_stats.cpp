#include "commands.h"

#include <vector>
#include <string>
using namespace std;

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
using namespace BamTools;
namespace po = boost::program_options;

void StatsCommand::getOptions()
{
}

int StatsCommand::runCommand()
{
    return 0;
}
