/*********************************************************************
 *
 * command_stats.cpp: Display statistics of a file.
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
#include <string>

#include "../algorithms/file_reader.h"
#include "../algorithms/repeatseq.h"

using namespace std;
namespace po = boost::program_options;

void RepeatseqCommand::getOptions() {
    options.add_options()
    ("reference,R", po::value<string>(), "FASTA reference")
    ("intervals,L", po::value<string>(), "Regions file")
    ("length,l", po::value<string>(), "use only a specific read length or range of read lengths (e.g. LENGTH or MIN:MAX)")
    ("before,B", po::value<int>(), "required number of reference matching bases BEFORE the repeat (default 3)")
    ("after,A", po::value<int>(), "required number of reference matching bases AFTER the repeat (default 3)")
    ("quality,Q", po::value<int>(), "minimum mapping quality for a read to be used for allele determination")
    ("multi,M", "exclude reads flagged with the XT:A:R tag")
    ("properlypaired,P", "exclude reads that are not properly paired (for PE reads only)")
    ("haploid,H", "assume a haploid rather than diploid genome")
    ("repeatseq", "write .repeatseq file containing additional information about reads")
    ("calls", "write .calls file")
    ("somaticA", po::value<string>(),"")
    ("somaticB", po::value<string>(),"")
    ("tag", po::value<string>(), "include user-defined tag in the output filename")
    ("flank", po::value<int>()->default_value(8), "number of flanking bases to output from each read")
    ("outfile", po::value<string>(), "When using stdin as input, set output filename")
    ;
    
    /*
     cout << endl << "-----------------------------------------------------------\n\n";
     cout << "RepeatSeq v" << VERSION << "\n\n";
     cout << "Usage:\t repeatseq [options] <in.bam> <in.fasta> <in.regions>\n\n";
     cout << "Options:";
     cout << " -r\t\tuse only a specific read length or range of read lengths (e.g. LENGTH or MIN:MAX)";
     cout << "\n\t -L\t\trequired number of reference matching bases BEFORE the repeat [3]";
     cout << "\n\t -R\t\trequired number of reference matching bases AFTER the repeat [3]";
     cout << "\n\t -M\t\tminimum mapping quality for a read to be used for allele determination";
     cout << "\n\t -multi\t\texclude reads flagged with the XT:A:R tag";
     cout << "\n\t -pp\t\texclude reads that are not properly paired (for PE reads only)";
     cout << "\n";
     cout << "\n\t -haploid\tassume a haploid rather than diploid genome";
     cout << "\n";
     cout << "\n\t -repeatseq\twrite .repeatseq file containing additional information about reads";
     cout << "\n\t -calls\t\twrite .calls file";
     cout << "\n\t -t\t\tinclude user-defined tag in the output filename";
     cout << "\n\t -o\t\tnumber of flanking bases to output from each read";
     cout << "\n";
     cout << endl << "-----------------------------------------------------------" << endl;
     */
}

int RepeatseqCommand::runCommand() {
    
    Repeatseq repeatseq;
    FileReader reader;
    reader.addSink(&repeatseq);

    //SETTINGS:
    string paramString;
    
    if(input_filenames.size() == 1 && input_filenames[0] == "stdin") {
        if(!vm.count("outfile")) {
            cerr << "Error. When using stdin as input, the --outfile parameter is required to name the output files." << endl;
            exit(-1);
        }
        reader.addFile("stdin");
        repeatseq.setOutputFilename(vm["outfile"].as<string>());
    } else {
        if(vm.count("outfile")) {
            cerr << "Warning. The --outfile parameter is only used when using stdin for input." << endl;
        }
        reader.addFiles(input_filenames);
        repeatseq.setOutputFilename(input_filenames[0]);
    }
    
    if(vm.count("somaticA") || vm.count("somaticB")) {
        if(vm.count("calls") || vm.count("repeatseq")) {
            cerr << "Somatic modes don't support calls or VCF output format. Quitting." << endl;
            exit(-1);
        }
        
        repeatseq.setMakeVcfFile(false);
        repeatseq.setMakeRepeatseqFile(true);
        if(vm.count("somaticA")) {
            repeatseq.setErrorModel(Repeatseq::ERROR_SOMATIC_A);
            repeatseq.setSomaticInput(vm["somaticA"].as<string>());
        } else {
            repeatseq.setErrorModel(Repeatseq::ERROR_SOMATIC_B);
            repeatseq.setSomaticInput(vm["somaticB"].as<string>());
        }
    }
    
    if(vm.count("intervals")) {
        repeatseq.setIntervalsFilename(vm["intervals"].as<string>());
    }else {
        cerr << "One intervals file is required. Quitting." << endl;
        exit(-1);
    }

    if(vm.count("reference")) {
        repeatseq.setFastaFilename(vm["reference"].as<string>());
    }else {
        cerr << "One FASTA reference file is required. Quitting." << endl;
        exit(-1);
    }

    if (vm.count("length")) {
        //read length select setting
        string range = vm["length"].as<string>();
        size_t firstColon = range.find(":");
        if (firstColon == string::npos){
            //no colon was found:
            int len = atoi(range.c_str());
            repeatseq.setLength(len, len);
            paramString += ".readlength";
            paramString += range.c_str();
        }
        else {
            // a colon was found:
            string first = range.substr(0, firstColon);
            string second = range.substr(firstColon+1,string::npos);
            int min_len = atoi(first.c_str());
            int max_len = atoi(second.c_str());
            repeatseq.setLength(min_len, max_len);
            paramString += ".readlength";
            paramString += first;
            paramString += "-";
            paramString += second;
        }
    }

    if (vm.count("tag")) {
        //append tag to paramString
        paramString += ".";
        paramString += vm["tag"].as<string>();
    }

    if (vm.count("flank")) {
        //set number of L/R chars to print
        repeatseq.setLRCharsToPrint(vm["flank"].as<int>());
    }

    if (vm.count("haploid")) {
        if(vm.count("somaticA") || vm.count("somaticB")) {
            cerr << "Haploid genomes aren't allowed with the somatic error model. Quitting." << endl;
            exit(-1);
        }
        repeatseq.setHaploid(true);
    }

    //FILTERS:
    if (vm.count("properlypaired")) {
        //Properly Paired Filter
        repeatseq.setProperlyPaired(true);
        paramString += ".ProperlyPaired";
    }

    if (vm.count("before")) {
        //Characters that much consecutively match to the LEFT for a read to be used
        repeatseq.setLeftFlank(vm["before"].as<int>());
        paramString += ".L";
        paramString += vm["before"].as<int>();
    }

    if (vm.count("after")) {
        //Characters that much consecutively match to the RIGHT for a read to be used
        repeatseq.setRightFlank(vm["after"].as<int>());
        paramString += ".R";
        paramString += vm["after"].as<int>();
    }

    if (vm.count("quality")) {
        //MINIMUM MapQuality Score
        repeatseq.setMinimumQuality(vm["quality"].as<int>());
        paramString += ".M";
        paramString += vm["quality"].as<int>();
    }

    if (vm.count("multi")) {
        //MULTI Filter (exclude read if XT:A:R tag is present)
        paramString += ".multi";
        repeatseq.setMulti(true);
    }

    if (vm.count("repeatseq"))
        repeatseq.setMakeRepeatseqFile(true);

    if (vm.count("calls"))
        repeatseq.setMakeCallsFile(true);
    
    repeatseq.setParamString(paramString);
    
    return reader.runChain();
}