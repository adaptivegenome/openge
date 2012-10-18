/*********************************************************************
 *
 * command_dedup.cpp: Remove duplicates in a file.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 14 April 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include "commands.h"
#include "../algorithms/file_reader.h"
#include "../algorithms/measure_coverage.h"
#include <iostream>

namespace po = boost::program_options;
using namespace std;

void CoverageCommand::getOptions()
{
    options.add_options()
    ("out,o", po::value<string>()->default_value("stdout"), "Output tsv filename. Omit for stdout.")
    ("verifymapping,V", "Verify mapping with read name (see docs)")
    ("omituncoveredbases", "Do not output info on bases with no coverage")
    ("binsize,b",po::value<int>()->default_value(100),"Bin size of output data")
    ("strict,S","When verifying mapping, use both coordinates (see docs)")
    ;
}

int CoverageCommand::runCommand()
{
    MeasureCoverage coverage;
    FileReader reader;
    
    coverage.setOutputFile(vm["out"].as<string>());
    coverage.setVerifyCorrectMapping(vm.count("verifymapping") > 0);
    coverage.setPrintZeroCoverageBases(vm.count("omituncoveredbases") == 0);
    coverage.setBinSize(vm["binsize"].as<int>());
    coverage.setStrict(vm.count("strict") > 0);
    
    reader.addSink(&coverage);

    reader.addFiles(input_filenames);

    reader.runChain();

    return 0;
}
