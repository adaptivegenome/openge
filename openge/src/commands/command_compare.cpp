/*********************************************************************
 *
 * command_compare.cpp: Compare the reads in multiple files, checking 
 *                      for differences
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 15 Aug 2012
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
#include <iomanip>
using namespace std;

#include <api/BamAlignment.h>
#include <api/BamReader.h>
using namespace BamTools;

#include "../algorithms/filter.h"
namespace po = boost::program_options;

class CompareElement {
    int chr;
    int position;
    vector<CigarOp> & cigar;
public:
    CompareElement(int chr, int position, vector<CigarOp> & cigar)
    : chr(chr)
    , position(position)
    , cigar(cigar)
    {}
    
    const CompareElement & operator=(const CompareElement & c) {
        chr = c.chr;
        position = c.position;
        cigar = c.cigar;
        return *this;
    }
    
    bool operator<(const CompareElement & e) const
    {
        if(chr != e.chr)
            return chr < e.chr;
        if(position != e.position)
            return position < e.position;
        if(cigar.size() != e.cigar.size())
            return cigar.size() < e.cigar.size();
        for(int i = 0; i < cigar.size(); i++) {
            if(cigar[i].Length != e.cigar[i].Length)
                return cigar[i].Length < e.cigar[i].Length;
            if(cigar[i].Type != e.cigar[i].Type)
                return cigar[i].Type < e.cigar[i].Type;
        }
        return false;
    }
};
            
bool regionContainsRead(const BamRegion & region, const BamAlignment & read) {
    return
        (read.RefID >= region.LeftRefID && read.RefID <= region.RightRefID) &&
    (read.Position >= region.LeftPosition && read.Position <= region.RightPosition);
}

void CompareCommand::getOptions()
{
    options.add_options()
    ("region,r", po::value<vector<string> >(), "Genomic regions to use.");
    ;
}

int CompareCommand::runCommand()
{
    //open ref file
    BamReader ref_reader;
    if(!ref_reader.Open(input_filenames[0])) {
        cerr << "Couldn't open input file" << input_filenames[0] << ". Aborting." << endl;
        exit(-1);
    }
    vector<BamRegion> regions;
    vector<string> region_strings;
    const SamSequenceDictionary sequences = ref_reader.GetHeader().Sequences;

    size_t longest_region_string = 0;
    if(vm.count("region")) {
        region_strings = vm["region"].as<vector<string> >();
        regions = vector<BamRegion>(region_strings.size());

        for(int region = 0; region < regions.size(); region++) {
            bool successful_parse = Filter::ParseRegionString(region_strings[region], regions[region], sequences);
            
            if(!successful_parse) {
                cerr << "Couldn't understand region string " << region_strings[region] << ". Exiting." << endl;
                exit(-1);
            }
            longest_region_string = max(longest_region_string, region_strings[region].size());
        }
    } else {
        region_strings.push_back("");
        regions.push_back(BamRegion(0, 0, INT_MAX, INT_MAX ));
    }
    
    vector<vector<CompareElement> > * ref_data_p = NULL;
    
    //we should only print filenames when not shooting results to wc -l or something.
    bool print_file_names = input_filenames.size() > 2;

    for(int file = 0; file < input_filenames.size(); file++) {
        
        vector<vector<CompareElement> > * compare_data_p = new vector<vector<CompareElement> >(regions.size());   //compare_data[region][element]
        if(file == 0)
            ref_data_p = compare_data_p;
        vector<vector<CompareElement> > & compare_data = *compare_data_p;
        vector<vector<CompareElement> > & ref_data = *ref_data_p;

        BamReader reader;
        if(!reader.Open(input_filenames[file])) {
            cerr << "Couldn't open input file " << input_filenames[file] << ". Aborting." << endl;
            exit(-1);
        }
        
        const SamSequenceDictionary compare_sequences = ref_reader.GetHeader().Sequences;

        if(compare_sequences != sequences) {
            cerr << "Sequence dictionaries between " << input_filenames[0] << " and " << input_filenames[file] << " differ. This may produce inconsistent results." << endl;
        }

        while(true) {
            BamAlignment * read = reader.GetNextAlignment();
            if(!read)
                break;
    
            for(int region = 0; region < regions.size(); region++) {
                if(regionContainsRead(regions[region], *read))
                    compare_data[region].push_back(CompareElement(read->RefID, read->Position, read->CigarData));
            }
        }

        //sort the vectors of each file
        for(int region = 0; region < regions.size(); region++)
            sort(compare_data[region].begin(), compare_data[region].end());
        
        //check for differences
        if(file != 0) {
            if(print_file_names)
                cerr << input_filenames[file] << ":\n";

            for(int region = 0; region < regions.size(); region++) {
                vector<CompareElement> common_elements;
                std::set_intersection(  compare_data[region].begin(), compare_data[region].end(), ref_data[region].begin(), ref_data[region].end(), std::back_inserter( common_elements )  );
                int additions = compare_data[region].size() - common_elements.size();
                int deletions = ref_data[region].size() - common_elements.size();
                if(additions || deletions) {
                    if(print_file_names)
                        cout << "   ";
                    cout << setw(longest_region_string) << region_strings[region] << "\t" << setw(8) << additions << " added\t" << setw(8) << deletions << " removed" << endl;
                }
            }
            if(file != 0)
                delete compare_data_p;
        }
    }
    delete ref_data_p;
    
    return 0;
}
