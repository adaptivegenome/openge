#include "commands.h"

#include <vector>
#include <string>
using namespace std;

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
using namespace BamTools;
namespace po = boost::program_options;

void HeadCommand::getOptions()
{
    options.add_options()
    ("out,o", po::value<string>()->default_value("stdout"), "Output filename. Omit for stdout.")
    ("count,n", po::value<size_t>()->default_value(10), "Number of alignments to copy");
}

int HeadCommand::runCommand()
{
    BamMultiReader reader;
    
    size_t count_limit = vm["limit"].as<size_t>();
    
    if(!reader.Open(input_filenames)) {
        cerr << "Error opening bam files." << endl;
        return -1;
    }
    
    BamWriter writer;
    string filename("stdout");
    
    if(vm.count("out") != 0)
        filename = vm["out"].as<string>();
    
    writer.Open(filename, reader.GetHeader(), reader.GetReferenceData());
    
    size_t count = 0;
    BamAlignment al;
    while(reader.GetNextAlignment(al) && count < count_limit) {
        writer.SaveAlignment(al);
        count++;
    }
    
    reader.Close();
    writer.Close();
    
    if(verbose)
      cerr << count << " alignments written to " << filename << endl;
    
    return count;
}
