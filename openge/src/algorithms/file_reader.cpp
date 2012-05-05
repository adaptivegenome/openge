#include "file_reader.h"

#include "api/BamMultiReader.h"
#include "api/SamReader.h"
using namespace BamTools;


int FileReader::runInternal()
{
    if(format == FORMAT_BAM)
    {
        BamMultiReader reader;
        
        if(!reader.Open(filenames)) {
            cerr << "Error opening BAM files." << endl;
            return -1;
        }
        
        header = reader.GetHeader();
        references = reader.GetReferenceData();
        open = true;
        
        BamAlignment * al;
        
        while(true)
        {
            al = reader.GetNextAlignment();
            if(!al)
                break;
            
            putOutputAlignment(al);
            count++;
        }
        
        reader.Close();
    } else if(format == FORMAT_SAM) {
        SamReader reader;
        
        if(filenames.size() != 1) {
            cerr << "Too many input files. SamReader only supports a single file input" << endl;
            return -1;
        }
        
        if(!reader.Open(filenames[0])) {
            cerr << "Error opening BAM files." << endl;
            return -1;
        }
        
        header = reader.GetHeader();
        references = reader.GetReferenceData();
        open = true;

        while(true)
        {
            BamAlignment al;
            if(!reader.GetNextAlignment(al))
                break;
            
            if(!sinks.empty())
                putOutputAlignment(new BamAlignment(al));
            count++;
        }
        
        reader.Close();
    } else 
        cerr << "Unrecognized output format." << endl;
    
    return 0;
}

void FileReader::addFile(std::string filename) 
{
    filenames.push_back(filename); 
}

void FileReader::addFiles(std::vector<std::string> filenames) 
{
    this->filenames.insert(this->filenames.end(), filenames.begin(), filenames.end()); 
}