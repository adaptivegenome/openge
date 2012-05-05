#include "commands.h"

#include "../algorithms/statistics.h"
#include "../algorithms/file_reader.h"

namespace po = boost::program_options;

void StatsCommand::getOptions()
{
    options.add_options()
    ("inserts,I", "Show detailed insert statistics")
    ;
}

int StatsCommand::runCommand()
{
    Statistics s;
    FileReader reader;
    
    reader.addFiles( input_filenames );
    
    s.showInsertSizeSummary(vm.count("inserts"));
    
    reader.addSink(&s);
    
    reader.runChain();
    return 0;
}
