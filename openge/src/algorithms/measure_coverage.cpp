/*********************************************************************
 *
 * filter.cpp:  Algorithm module that filters a stream of reads.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 22 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************
 *
 * Filter algorithm module class. Filters a stream of reads by the 
 * region of the genome, or by count, or both. Some code based on 
 * Bamtools' merge code, from bamtools_merge.cpp.
 *
 *********************************************************************/

#include "measure_coverage.h"

#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cassert>

using namespace BamTools;
using namespace std;

MeasureCoverage::MeasureCoverage()
: verify_mapping(false)
, print_zero_cover_bases(false)
, strict(false)
, binsize(100)
{ }

int MeasureCoverage::runInternal()
{
    SamHeader header = this->getHeader();
    int total_length = 0;
    int num_correct_maps = 0;
    int num_skipped_reads = 0;
    bool overflow = false;
    try {
        if(verbose)
            cerr << "Setting up coverage counting structures" << endl;
        for(SamSequenceConstIterator i = header.Sequences.Begin(); i != header.Sequences.End(); i++) {
            size_t length = (atol(i->Length.c_str()) + binsize - 1) / binsize;
            total_length += length;
            
            coverage_map[i->Name].reserve(length);
            coverage_map[i->Name].insert(coverage_map[i->Name].begin(), length, 0);
            
            if(verify_mapping) {
                correctness_map[i->Name].reserve(length);
                correctness_map[i->Name].insert(correctness_map[i->Name].begin(), length, 0);
            }
        }
    } catch (bad_alloc& ba) {
        cerr << "Ran out of memory allocating tracking structures. Try using a higher value for --binsize." << endl;
    }
    
    if(verbose) 
        cerr << "Measuring coverage" << endl;
    while(true) {
        BamAlignment * al = getInputAlignment();
        
        if(!al)
            break;
        
        if(al->RefID == -1 || al->Position == -1) {
            num_skipped_reads++;
            putOutputAlignment(al);
            continue;
        }

        string & chr = header.Sequences[al->RefID].Name;
        assert(coverage_map.count(chr) > 0);
        vector<unsigned int> & chr_vector = coverage_map[chr];

        for(int i = al->Position; i <= al->Position + al->Length; i++) {
            chr_vector[i/binsize]++;
            
            if(chr_vector[i/binsize] == UINT_MAX)
                overflow = true;
        }

        // we are scanning read names that look like
        // "@chr1_80429992_80429614_1_0_0_0_0:0:0_0:0:0_2ed4"
        // and define correctness as chromosome name matching, and 
        // position being within 5 bases of either of the two numbers
        // in the name (8042xxxxx in this example
        if(verify_mapping) {
            char name[32] = {0};
            int n1, n2;
            const char * name_buffer = al->Name.c_str();
            const char * first_underscore = strchr(al->Name.c_str(), '_');
            strncpy(name, name_buffer, first_underscore - name_buffer);
            if( first_underscore ) {
                int ret = sscanf(first_underscore, "_%d_%d", &n1, &n2);
                bool n1_match = (5 >= abs(n1 - al->Position));
                bool n2_match = (5 >= abs(n2 - al->Position));
                bool n_match = strict ? n1_match : (n1_match || n2_match);
                if( 0 == strcmp(name, chr.c_str()) && ret == 2
                   && n_match
               ) {
                
                    assert(correctness_map.count(chr) > 0);
                    vector<unsigned int> & correct_vector = correctness_map[chr];
                    for(int i = al->Position; i <= al->Position + al->Length; i++) {
                        correct_vector[i/binsize]++;
                        
                        if(correct_vector[i/binsize] == UINT_MAX)
                            overflow = true;
                    }
                
                    num_correct_maps++;
                }
            }
        }
        
        putOutputAlignment(al);
    }

    if(num_skipped_reads)
        cerr << "Skipped " << num_skipped_reads << " unmapped reads." << endl;
    
    //now print coverage:
    cerr << "Average coverage:" << endl;
    for(map<string, vector<unsigned int> >::const_iterator vec = coverage_map.begin(); vec != coverage_map.end(); vec++) {
        int64_t total_coverage = 0;
        for(int i = 0; i < vec->second.size(); i++)
            total_coverage += vec->second[i];
        cerr << "   " << setw(20) << vec->first << ": " << setw(8) << double( total_coverage) / atoi(header.Sequences[vec->first].Length.c_str()) << "x" << endl;
    }
    
    if(verbose && verify_mapping)
        cerr << "Found " << num_correct_maps << " / " << write_count << " (" << 100. * num_correct_maps / write_count << " %) reads were correctly mapped." << endl;

    
    if(out_filename != "stdout") {
        if(verbose) {
            if(print_zero_cover_bases)
                cerr << "Writing coverage file (expect " << total_length/binsize << " lines)" << endl;
            else
                cerr << "Writing coverage file" << endl;
        }

        ofstream outfile(out_filename.c_str());
        
        if(verify_mapping)
            outfile << "chromosome\tposition\tcoverage\tcorrect_maps\n";
        else
            outfile << "chromosome\tposition\tcoverage\n";
        int write_ct = 0;
        for(map<string, vector<unsigned int> >::const_iterator vec = coverage_map.begin(); vec != coverage_map.end(); vec++) {
            const unsigned int * count_data = &(vec->second[0]);
            const unsigned int * correctness_data = &(correctness_map[vec->first][0]);
            
            for(int i = 0; i < vec->second.size(); i++) {
                if(print_zero_cover_bases || count_data[i] != 0) {
                    if(verify_mapping)
                        outfile << vec->first << "\t" << binsize * i+1 << "\t" << count_data[i] << "\t" << correctness_data[i] << "\n";
                    else
                        outfile << vec->first << "\t" << binsize * i+1 << "\t" << count_data[i] << "\n";
                }

                write_ct += binsize;
                if(verbose && write_ct % 5000000 == 0)
                    cerr << "\rWriting " << 100. * i / vec->second.size() << "% done";
            }
        }
        if(verbose)
            cerr << "\rWriting 100% done" << endl;
    }
    
    if(overflow)
        cerr << "Error: at least one overflow occurred when measuring coverage- try reducing binsizes." << endl;

    return 0;
}
