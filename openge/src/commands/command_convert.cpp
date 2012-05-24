#include "commands.h"

#include "../algorithms/file_reader.h"
#include "../algorithms/file_writer.h"

namespace po = boost::program_options;
using namespace std;

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

    if(verbose)
        cerr << "Converting " << filename_in << " from SAM to BAMfile " << filename_out << "." << endl;
    
    FileReader reader;
    FileWriter writer;
    reader.addSink(&writer);
    
    reader.setFormat(FileReader::FORMAT_SAM);
    reader.addFiles(input_filenames);
    writer.setFilename(filename_out);
    
    reader.runChain();
    
    if(verbose)
        cout << writer.getCount() << " alignments written to " << filename_out << endl;
    
    return 0;
}