/*********************************************************************
 *
 * command_view.cpp: Filter a file
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 21 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include "commands.h"

#include <vector>
#include <string>
using namespace std;

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif

#include "../algorithms/file_reader.h"
#include "../algorithms/file_writer.h"
#include "../algorithms/filter.h"
using namespace BamTools;
namespace po = boost::program_options;

void ViewCommand::getOptions()
{
    options.add_options()
    ("out,o", po::value<string>()->default_value("stdout"), "Output filename. Omit for stdout.")
    ("count,n", po::value<size_t>(), "Number of alignments to copy")
    ("mapq,q", po::value<int>(), "Minimum map quality")
    ("length,l", po::value<string>(), "Range of acceptable read lengths")
    ("trimbegin,B", po::value<int>(), "Trim the beginning of read by arg bases")
    ("trimend,E", po::value<int>(), "Trim the beginning of read by end bases")
    ("region,r", po::value<string>(), "Genomic region to use.");
}

int ViewCommand::runCommand()
{
    string filename_in = input_filenames[0];
    string filename_out = vm["out"].as<string>();
    
    if(input_filenames.size() > 1)
        cerr << "More than one input filename provided - only using " << filename_in << "." << endl;

    FileReader reader;
    Filter filter;
    FileWriter writer;
    reader.addSink(&filter);
    filter.addSink(&writer);
    
    if(vm.count("count")) {
        size_t count_limit = vm["count"].as<size_t>();
        filter.setCountLimit(count_limit);
    }

    if(vm.count("mapq")) {
        int mapq = vm["mapq"].as<int>();
        filter.setQualityLimit(mapq);
    }

    if(vm.count("format"))
        writer.setFormat(vm["format"].as<string>());

    if(vm.count("length"))
        filter.setReadLengths(vm["length"].as<string>());
    
    if(vm.count("trimbegin") > 0)
        filter.setTrimBeginLength(vm["trimbegin"].as<int>());

    if(vm.count("trimend") > 0)
        filter.setTrimEndLength(vm["trimend"].as<int>());

    bool hasregion = vm.count("region") != 0;

    if(hasregion) {
        std::string region_string = vm["region"].as<string>();
        filter.setRegion(region_string);
    }
    
    reader.addFiles(input_filenames);
    writer.setFilename(filename_out);
    writer.addProgramLine(command_line);

    if((vm.count("trimbegin") > 0 || vm.count("trimend") > 0) && writer.getFileFormat() != FORMAT_FASTQ) {
        cerr << "Trimming reads is only supported for the FASTQ format at this time. Aborting." << endl;
        exit(-1);
    }
    
    if(filename_out == "stdout")
        writer.setDefaultFormat(FORMAT_SAM);
    
    reader.runChain();
    
    return 0;
}
