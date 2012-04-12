#include "commands.h"

using namespace std;

#include <api/SamReader.h>
#include <api/BamWriter.h>
using namespace BamTools;
namespace po = boost::program_options;

void ConvertCommand::getOptions()
{
    options.add_options()
    ("out,o", po::value<string>()->default_value("stdout"), "Output filename. Omit for stdout.")
    ;
}

int ConvertCommand::runCommand()
{
    string filename_in = input_filenames[0];
    string filename_out = vm["out"].as<string>();
    
    if(input_filenames.size() > 1)
        cerr << "More than one input filename provided - only using " << filename_in << "." << endl;
    
    if(verbose)
        cerr << "Converting " << filename_in << " from SAM to BAMfile " << filename_out << "." << endl;
    
    SamReader reader;
    
    if(!reader.Open(filename_in)) {
        cerr << "Error opening " << filename_in << " for reading" << endl;
        return -1;
    }
    
    BamWriter writer;
    
    if(!writer.Open(filename_out, reader.GetHeader(), reader.GetRefData()))
    {
        cerr << "Error opening " << filename_out << " for writing" << endl;
        return -1;
    }
    
    size_t count = 0;
    while(true)
    {
        BamAlignment al;
        if(!reader.GetNextAlignment(al))
            break;
        writer.SaveAlignment(al);
        count++;
    }
    
    reader.Close();
    writer.Close();
    
    if(verbose)
        cout << count << " alignments written to " << filename_out << endl;
    
    return count;
}