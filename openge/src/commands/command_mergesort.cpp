#include "commands.h"
namespace po = boost::program_options;

#include "../algorithms/file_reader.h"
#include "../algorithms/read_sorter.h"
#include "../algorithms/mark_duplicates.h"
#include "../algorithms/filter.h"
#include "../algorithms/file_writer.h"

#include <string>
using namespace std;


////////////////////////////
// Command interface

void MergeSortCommand::getOptions()
{
    
    options.add_options()
    ("out,o", po::value<string>()->default_value("stdout"), "Output filename. Omit for stdout.")
    ("compression,c", po::value<int>()->default_value(6), "Compression level of the output. Valid 0-9.")
    ("region,r", po::value<string>(), "Genomic region to use.")
    ("byname,b", "Sort by name. Otherwise, sorts by position.")
    ("n,n", po::value<int>()->default_value(5e5), "Alignments per temp file.")
    ("compresstempfiles,C", "Compress temp files. By default, uncompressed")
    ("markduplicates,d", "Mark duplicates after sorting.")
    ("removeduplicates,R", "Remove duplicates.")
    ;
}

int MergeSortCommand::runCommand()
{
    bool do_mark_duplicates = vm.count("markduplicates") != 0;
    bool do_remove_duplicates = vm.count("removeduplicates") != 0;
    if(do_remove_duplicates)
        do_mark_duplicates = true;

    bool sort_by_names = vm.count("byname") != 0;
    int compression_level = vm["compression"].as<int>();
    bool compresstempfiles = vm.count("compresstempfiles") != 0;
    int alignments_per_tempfile = vm["n"].as<int>();
    
    
    //The chain for this command goes something like this:
    //Reader->Filter->Sort->MarkDuplicates->Writer (->BlackHole)
    //
    // Filter is omitted if there is no region
    // MarkDuplicates is omitted if it isn't required
    // BlackHole is automatically included when run.
    
    FileReader reader;
    Filter filter;
    ReadSorter sort_reads;
    MarkDuplicates mark_duplicates;
    FileWriter writer;
    
    if(vm.count("region")) {
        string region = vm["region"].as<string>();
        filter.setRegion(region);
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

    sort_reads.setSortBy(sort_by_names ? SORT_NAME : SORT_POSITION);
    sort_reads.setCompressTempFiles(compresstempfiles);
    sort_reads.setAlignmentsPerTempfile(alignments_per_tempfile);

    reader.addFiles(input_filenames);
    writer.setFilename(vm["out"].as<string>());
    writer.setCompressionLevel(compression_level);
    
    return writer.runChain();
}

