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
#include "../util/picard_structures.h"
#include "../algorithms/mark_duplicates.h"
#include "../algorithms/file_writer.h"
#include "../algorithms/file_reader.h"
#include "../algorithms/sorted_merge.h"
#include "../algorithms/split_by_chromosome.h"
#include <iostream>

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
using namespace BamTools;
namespace po = boost::program_options;
using namespace std;

void DedupCommand::getOptions()
{
    options.add_options()
    ("out,o", po::value<string>()->default_value("stdout"), "Output filename. Omit for stdout.")
    ("remove,r", "Remove duplicates")
    ("compression,c", po::value<int>()->default_value(6), "Compression level of the output. Valid 0-9.")
    ("nosplit","Do not split by chromosome (for speed) when processing")
    ;
}

int DedupCommand::runCommand()
{
    bool do_remove_duplicates = vm.count("remove") != 0;
    bool no_split = vm.count("nosplit") != 0;
    int compression_level = vm["compression"].as<int>();

    if(no_split && verbose)
        cerr << "Disabling split-by-chromosome." << endl;

    int num_chains = min(12,ThreadPool::availableCores()/2);

    if(nothreads || no_split || num_chains <= 1)
    {
        FileReader reader;
        MarkDuplicates mark_duplicates;
        FileWriter writer;

        reader.setLoadStringData(false);

        reader.addSink(&mark_duplicates);

        mark_duplicates.addSink(&writer);
        mark_duplicates.removeDuplicates = do_remove_duplicates;

        reader.addFiles(input_filenames);
        writer.setFilename(vm["out"].as<string>());
        writer.setCompressionLevel(compression_level);
        
        return writer.runChain();
    } else {
        FileReader reader;
        SortedMerge merge;
        SplitByChromosome split;
        FileWriter writer;

        vector<MarkDuplicates *> duplicate_markers;
        
        //read-filter-split
        reader.setLoadStringData(false);
        reader.addSink(&split);
        
        //merge-write
        merge.addSink(&writer);
        
        // split-sort-dedup-merge
        // each iteration of this loop forms one chain inside the split
        for(int ctr = 0; ctr < num_chains; ctr++)
        {  
            MarkDuplicates * mark_duplicates = new MarkDuplicates;
            duplicate_markers.push_back(mark_duplicates);
            merge.addSource(mark_duplicates);
            mark_duplicates->removeDuplicates = do_remove_duplicates;
            
            char filename[48];
            sprintf(filename, "/dedup_%8x.bam",  (uint32_t)(0xffffffff & (uint64_t)mark_duplicates));
            mark_duplicates->setBufferFileName(tmpdir + string(filename));

            
            split.addSink(mark_duplicates);
        }

        reader.addFiles(input_filenames);
        writer.setFilename(vm["out"].as<string>());
        writer.setCompressionLevel(compression_level);
        
        int ret = writer.runChain();
        
        //clean up allocated objects
        for(int ctr = 0; ctr < num_chains; ctr++)
            delete duplicate_markers[ctr];
        
        return ret;
    }
}
