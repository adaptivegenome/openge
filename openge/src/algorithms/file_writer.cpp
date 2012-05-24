#include "file_writer.h"

#include "api/BamWriter.h"
using namespace BamTools;
using namespace std;

#ifdef __linux__
#include <sys/prctl.h>
#endif

int FileWriter::runInternal()
{
#ifdef __linux__
    prctl(PR_SET_NAME,"am_FileWriter",0,0,0);
#endif
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
        count++;
    }
    
    writer.Close();
    
    if(isVerbose())
        cerr << "Wrote " << count << " files to " << filename << endl;
    
    return 0;
}