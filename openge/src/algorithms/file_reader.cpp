#include "file_reader.h"

#include "api/BamMultiReader.h"
using namespace BamTools;


int FileReader::runInternal()
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
    }
    
    reader.Close();
    
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