#include "commands.h"

#include <vector>
#include <string>
using namespace std;

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
using namespace BamTools;
namespace po = boost::program_options;

void ViewCommand::getOptions()
{
    options.add_options()
    ("out,o", po::value<string>()->default_value("stdout"), "Output filename. Omit for stdout.")
    ("count,n", po::value<size_t>(), "Number of alignments to copy")
    ("region,r", po::value<string>(), "Genomic region to use.");
}

// this has been copied from bamtools utilities, since it isn't in the API. Original file is bamtools_utilities.cpp .
bool ParseRegionString(const string& regionString, const BamMultiReader& reader, BamRegion& region)
{
    // -------------------------------
    // parse region string
    
    // check first for empty string
    if ( regionString.empty() ) 
        return false;   
    
    // non-empty string, look for a colom
    size_t foundFirstColon = regionString.find(':');
    
    // store chrom strings, and numeric positions
    string startChrom;
    string stopChrom;
    int startPos;
    int stopPos;
    
    // no colon found
    // going to use entire contents of requested chromosome 
    // just store entire region string as startChrom name
    // use BamReader methods to check if its valid for current BAM file
    if ( foundFirstColon == string::npos ) {
        startChrom = regionString;
        startPos   = 0;
        stopChrom  = regionString;
        stopPos    = -1;
    }
    
    // colon found, so we at least have some sort of startPos requested
    else {
        
        // store start chrom from beginning to first colon
        startChrom = regionString.substr(0,foundFirstColon);
        
        // look for ".." after the colon
        size_t foundRangeDots = regionString.find("..", foundFirstColon+1);
        
        // no dots found
        // so we have a startPos but no range
        // store contents before colon as startChrom, after as startPos
        if ( foundRangeDots == string::npos ) {
            startPos   = atoi( regionString.substr(foundFirstColon+1).c_str() ); 
            stopChrom  = startChrom;
            stopPos    = -1;
        } 
        
        // ".." found, so we have some sort of range selected
        else {
            
            // store startPos between first colon and range dots ".."
            startPos = atoi( regionString.substr(foundFirstColon+1, foundRangeDots-foundFirstColon-1).c_str() );
            
            // look for second colon
            size_t foundSecondColon = regionString.find(':', foundRangeDots+1);
            
            // no second colon found
            // so we have a "standard" chrom:start..stop input format (on single chrom)
            if ( foundSecondColon == string::npos ) {
                stopChrom  = startChrom;
                stopPos    = atoi( regionString.substr(foundRangeDots+2).c_str() );
            }
            
            // second colon found
            // so we have a range requested across 2 chrom's
            else {
                stopChrom  = regionString.substr(foundRangeDots+2, foundSecondColon-(foundRangeDots+2));
                stopPos    = atoi( regionString.substr(foundSecondColon+1).c_str() );
            }
        }
    }
    
    // -------------------------------
    // validate reference IDs & genomic positions
    
    const RefVector references = reader.GetReferenceData();
    
    // if startRefID not found, return false
    int startRefID = reader.GetReferenceID(startChrom);
    if ( startRefID == -1 ) return false;
    
    // startPos cannot be greater than or equal to reference length
    const RefData& startReference = references.at(startRefID);
    if ( startPos >= startReference.RefLength ) return false;
    
    // if stopRefID not found, return false
    int stopRefID = reader.GetReferenceID(stopChrom);
    if ( stopRefID == -1 ) return false;
    
    // stopPosition cannot be larger than reference length
    const RefData& stopReference = references.at(stopRefID);
    if ( stopPos > stopReference.RefLength ) return false;
    
    // if no stopPosition specified, set to reference end
    if ( stopPos == -1 ) stopPos = stopReference.RefLength;
    
    // -------------------------------
    // set up Region struct & return
    
    region.LeftRefID     = startRefID;
    region.LeftPosition  = startPos;
    region.RightRefID    = stopRefID;;
    region.RightPosition = stopPos;
    return true;
}

int ViewCommand::runCommand()
{
    BamMultiReader reader;
    
    size_t count_limit = SIZE_T_MAX;
    
    if(vm.count("count"))
        count_limit = vm["count"].as<size_t>();
    
    if(!reader.Open(input_filenames)) {
        cerr << "Error opening bam files." << endl;
        return -1;
    }
    
    BamWriter writer;
    string filename("stdout");
    
    if(vm.count("out") != 0)
        filename = vm["out"].as<string>();
    
    writer.Open(filename, reader.GetHeader(), reader.GetReferenceData());
    
    size_t count = 0;
    
    bool hasregion = vm.count("region") != 0;
    
    // if no region specified, store entire contents of file(s)
    if ( !hasregion ) {
        BamAlignment al;
        while ( reader.GetNextAlignmentCore(al) && count < count_limit ) {
            writer.SaveAlignment(al);
            count++;
        }
    }
    
    // otherwise attempt to use region as constraint
    else {
        string region_string = vm["region"].as<string>();
        
        if(verbose)
            cerr << "Filtering to region " << region_string << endl;
        
        // if region string parses OK
        BamRegion region;
        if (ParseRegionString(region_string, reader, region) ) {
            
            // attempt to find index files
            reader.LocateIndexes();
            
            // if index data available for all BAM files, we can use SetRegion
            if ( reader.HasIndexes() ) {
                
                // attempt to use SetRegion(), if failed report error
                if ( !reader.SetRegion(region.LeftRefID,
                                       region.LeftPosition,
                                       region.RightRefID,
                                       region.RightPosition) )
                {
                    cerr << "bamtools merge ERROR: set region failed. Check that REGION describes a valid range"
                    << endl;
                    reader.Close();
                    return false;
                }
                
                // everything checks out, just iterate through specified region, storing alignments
                BamAlignment al;
                while (reader.GetNextAlignmentCore() && count < count_limit ) {
                    writer.SaveAlignment(al);
                    count++;
                }
            }
            
            // no index data available, we have to iterate through until we
            // find overlapping alignments
            else {
                BamAlignment al;
                while ( reader.GetNextAlignmentCore(al)  && count < count_limit) {
                    if ( (al.RefID >= region.LeftRefID)  && ( (al.Position + al.Length) >= region.LeftPosition ) &&
                        (al.RefID <= region.RightRefID) && ( al.Position <= region.RightPosition) ) {
                        writer.SaveAlignment(al);
                        count++;
                    }
                }
            }
        }
        
        // error parsing REGION string
        else {
            cerr << "ERROR: could not parse region'" << region_string << "'" << endl;
            cerr << "Check that region description is in valid format (see documentation) and that the coordinates are valid"
            << endl;
            reader.Close();
            return false;
        }
    }
    
    reader.Close();
    writer.Close();
    
    if(verbose)
        cerr << count << " alignments written to " << filename << endl;
    
    return 0;
}
