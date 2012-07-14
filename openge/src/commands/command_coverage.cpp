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

using namespace BamTools;
namespace po = boost::program_options;
using namespace std;

void CoverageCommand::getOptions()
{
    options.add_options()
    ("out,o", po::value<string>()->default_value("stdout"), "Output csv filename. Omit for stdout.")
    ("verifymapping,V", "Verify mapping with read name (see docs)")
    ("omituncoveredbases", "Do not output info on bases with no coverage")
    ;
}

int CoverageCommand::runCommand()
{
    MeasureCoverage coverage;
    FileReader reader;
    
    coverage.setOutputFile(vm["out"].as<string>());
    coverage.setVerifyCorrectMapping(vm.count("verifymapping") > 0);
    coverage.setPrintZeroCoverageBases(vm.count("omituncoveredbases") == 0);
    
    reader.addSink(&coverage);

    reader.addFiles(input_filenames);

    reader.runChain();

    return 0;
}
