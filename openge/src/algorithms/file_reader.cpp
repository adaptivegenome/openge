#include "file_reader.h"

#include "api/BamMultiReader.h"
#include "../util/sam_reader.h"
using namespace BamTools;
using namespace std;

FileReader::file_format_t FileReader::deduceFileFormat()
{
    FILE * fp = fopen(filenames[0].c_str(), "rb");
    
    if(!fp) {
        cerr << "Couldn't open file " << filenames[0] << endl;
        return FORMAT_UNKNOWN;
    }

    unsigned char data[2];
    if(2 != fread(data, 2, 1, fp)) {
        cerr << "Couldn't read from file " << filenames[0] << endl;
        return FORMAT_UNKNOWN;
    }

    fclose(fp);
    
    if(data[0] == '@')
        return FORMAT_SAM;
    
    if(data[0] == 31 && data[1] == 139)
        return FORMAT_BAM;
    
    return FORMAT_UNKNOWN;
}

int FileReader::runInternal()
{
    if(!format_specified)
        format = deduceFileFormat();

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
            if(load_string_data)
                al = reader.GetNextAlignment();
            else
                al = reader.GetNextAlignmentCore();

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
            cerr << "Error opening SAM file." << endl;
            return -1;
        }
        
        header = reader.GetHeader();
        references = reader.GetReferenceData();
        open = true;

        
        BamAlignment * al = NULL;
        while(true)
        {
            al = reader.GetNextAlignment();
            
            if(NULL == al)
                break;
            
            if(!sinks.empty())
                putOutputAlignment(al);
            count++;
        }

        reader.Close();
    } else {
        cerr << "FileReader couldn't detect file format. Aborting." << endl;
        return -1;
    }

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