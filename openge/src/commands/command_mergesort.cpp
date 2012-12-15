/*********************************************************************
 *
 * command_mergesort.cpp: Perform merging, sorting, and mark duplicates.
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
namespace po = boost::program_options;

#include "../algorithms/file_reader.h"
#include "../algorithms/read_sorter.h"
#include "../algorithms/mark_duplicates.h"
#include "../algorithms/filter.h"
#include "../algorithms/file_writer.h"
#include "../algorithms/split_by_chromosome.h"
#include "../algorithms/sorted_merge.h"

#include <string>
using namespace std;


////////////////////////////
// Command interface

void MergeSortCommand::getOptions()
{
    
    options.add_options()
    ("out,o", po::value<string>()->default_value("stdout"), "Output filename. Omit for stdout.")
    ("region,r", po::value<string>(), "Genomic region to use.")
    ("mapq,q", po::value<int>(), "Minimum map quality allowed in reads")
    ("byname,b", "Sort by name. Otherwise, sorts by position.")
    ("n,n", po::value<int>()->default_value(5e5), "Alignments per temp file.")
    ("compresstempfiles,C", "Compress temp files. By default, uncompressed")
    ("markduplicates,M", "Mark duplicates after sorting.")
    ("removeduplicates,R", "Remove duplicates.")
    ;
}

int MergeSortCommand::runCommand()
{
    bool do_mark_duplicates = vm.count("markduplicates") != 0;
    bool do_remove_duplicates = vm.count("removeduplicates") != 0;
    bool no_split = vm.count("nosplit") != 0;
    if(do_remove_duplicates)
        do_mark_duplicates = true;

    bool sort_by_names = vm.count("byname") != 0;
    int compression_level = vm["compression"].as<int>();
    bool compresstempfiles = vm.count("compresstempfiles") != 0;
    int alignments_per_tempfile = vm["n"].as<int>();
    
    if(no_split && verbose)
        cerr << "Disabling split-by-chromosome." << endl;
    
    int num_chains = min(12,OGEParallelismSettings::getNumberThreads ()/2);
    
    if(nothreads || no_split || !do_mark_duplicates || num_chains <= 1)
    {
        //The chain for this command goes something like this:
        //Reader->Filter->Sort->MarkDuplicates->Writer (->BlackHole)
        //
        // Filter is omitted if there is no region
        // MarkDuplicates is omitted if it isn't required
        // BlackHole is automatically included when run.
        
        FileReader reader;
        Filter filter;
        ReadSorter sort_reads(tmpdir);
        MarkDuplicates mark_duplicates(tmpdir);
        FileWriter writer;
        
        reader.setLoadStringData(false);
        
        if(vm.count("region") || vm.count("mapq")) {
            if(vm.count("region"))
                filter.setRegion(vm["region"].as<string>());
            if(vm.count("mapq"))
                filter.setQualityLimit(vm["mapq"].as<int>());
            reader.addSink(&filter);
            filter.addSink(&sort_reads);
        } else {
            reader.addSink(&sort_reads);
        }

        if(do_mark_duplicates) {
            sort_reads.addSink(&mark_duplicates);
            mark_duplicates.addSink(&writer);
            mark_duplicates.removeDuplicates = do_remove_duplicates;
        }
        else {
            sort_reads.addSink(&writer);
        }

        sort_reads.setSortBy(sort_by_names ? BamHeader::SORT_QUERYNAME : BamHeader::SORT_COORDINATE);
        sort_reads.setCompressTempFiles(compresstempfiles);
        sort_reads.setAlignmentsPerTempfile(alignments_per_tempfile);

        reader.addFiles(input_filenames);
        writer.setFilename(vm["out"].as<string>());
        if(!vm.count("nopg"))
            writer.addProgramLine(command_line);
        writer.setCompressionLevel(compression_level);
        if(vm.count("format"))
            writer.setFormat(vm["format"].as<string>());
        
        return writer.runChain();
    } else {
        //The chain for this command goes something like this:
        //Reader->Filter->Sort->Split->MarkDuplicates(multiple)->Merge->Writer (->BlackHole)
        //
        // Filter is omitted if there is no region
        // BlackHole is automatically included when run.
        
        FileReader reader;
        Filter filter;
        SortedMerge merge;
        SplitByChromosome split;
        FileWriter writer;
        ReadSorter sort_reads(tmpdir);

        vector<MarkDuplicates *> duplicate_markers;
        
        //read-filter-split
        reader.setLoadStringData(false);
        if(vm.count("region")) {
            string region = vm["region"].as<string>();
            filter.setRegion(region);
            reader.addSink(&filter);
            filter.addSink(&sort_reads);
        } else {
            reader.addSink(&sort_reads);
        }
        
        sort_reads.addSink(&split);

        //merge-write
        merge.addSink(&writer);

        // split-sort-dedup-merge
        // each iteration of this loop forms one chain inside the split
        for(int ctr = 0; ctr < num_chains; ctr++)
        {  
            MarkDuplicates * mark_duplicates = new MarkDuplicates(tmpdir);
            duplicate_markers.push_back(mark_duplicates);
            merge.addSource(mark_duplicates);
            split.addSink(mark_duplicates);

            mark_duplicates->removeDuplicates = do_remove_duplicates;
        }

        sort_reads.setSortBy(sort_by_names ? BamHeader::SORT_QUERYNAME : BamHeader::SORT_COORDINATE);
        sort_reads.setCompressTempFiles(compresstempfiles);
        sort_reads.setAlignmentsPerTempfile(alignments_per_tempfile);
        
        reader.addFiles(input_filenames);
        writer.setFilename(vm["out"].as<string>());
        if(!vm.count("nopg"))
            writer.addProgramLine(command_line);
        writer.setCompressionLevel(compression_level);
        
        int ret = writer.runChain();
        
        //clean up allocated objects
        for(int ctr = 0; ctr < num_chains; ctr++)
            delete duplicate_markers[ctr];

        return ret;
    }
}

