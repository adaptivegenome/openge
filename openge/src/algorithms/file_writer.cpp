#include "file_writer.h"

#include "api/BamWriter.h"
using namespace BamTools;


int FileWriter::runInternal()
{
    BamWriter writer;
    
    if(!writer.Open(filename, getHeader(), getReferences())) {
        cerr << "Error opening BAM file to write." << endl;
        return -1;
    }
    
    writer.SetCompressionLevel(compression_level);
    
    BamAlignment * al;
    
    while(true)
    {
        al = getInputAlignment();
        if(!al)
            break;

        writer.SaveAlignment(*al);
        putOutputAlignment(al);
    }
    
    writer.Close();
    
    return 0;
}