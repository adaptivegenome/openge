/*********************************************************************
 *
 * file_writer.cpp:  Algorithm module writes a BAM or SAM file.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 31 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************
 *
 * File writer algorithm module. Writes a stream of reads to a BAM
 * file. Eventually this will be extended to support SAM or CRAM
 * formats.
 *
 *********************************************************************/

#include "file_writer.h"
#include "openge_constants.h"

#include "../util/sam_writer.h"
#include "../util/fastq_writer.h"

#include "api/BamWriter.h"
using namespace BamTools;
using namespace std;

//from http://stackoverflow.com/questions/874134/find-if-string-endswith-another-string-in-c
bool hasEnding (std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

void FileWriter::setFormat(const std::string & format_name)
{
    file_format_t ret = detectFileFormatFromFilename(format_name);
    
    if(ret == FORMAT_UNKNOWN) {
        cerr << "Unknown file format specified: " << format_name << ". Aborting." << endl;
        exit(-1);
    }
    
    setFormat(ret);
}

file_format_t FileWriter::getFileFormat() {
    file_format_t determined_format = file_format;
    
    if(determined_format == FORMAT_UNKNOWN) {
        determined_format = detectFileFormatFromFilename(filename);
        if(determined_format == FORMAT_UNKNOWN)
            determined_format = default_file_format;
    }
    
    return determined_format;
}

int FileWriter::runInternal()
{
    ogeNameThread("am_FileWriter");
    
    file_format = getFileFormat();
    SamHeader header = getHeader();
    
    if(command_line_options.size() > 0) {
        SamProgram pg;
        pg.ID = string("openge");
	pg.Version = string(OPENGE_VERSION_STRING);

        for(int i = 2; header.Programs.Contains( pg.ID); i++) {
            stringstream s;
            s << "openge-" << i;
            pg.ID = s.str();
        }

        pg.CommandLine = command_line_options;
        header.Programs.Add(pg);
    }

    switch(file_format) {
        case FORMAT_SAM:
            {
                SamWriter writer;
                
                if(!writer.Open(filename, header)) {
                    cerr << "Error opening BAM file to write." << endl;
                    exit(-1);
                }
                
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
                
            }
            break;
        case FORMAT_FASTQ:
            {
                FastqWriter writer;
                
                if(!writer.Open(filename, header)) {
                    cerr << "Error opening FASTQ file to write." << endl;
                    exit(-1);
                }
                
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
                
            }
            break;
        case FORMAT_BAM:
            {
                BamWriter writer;

                writer.SetCompressionMode(BamWriter::Compressed);
                writer.SetCompressionLevel(compression_level);
                
                RefVector references;
                
                for(SamSequenceConstIterator i = header.Sequences.Begin(); i != header.Sequences.End(); i++) {
                    RefData d;
                    d.RefName = i->Name;
                    d.RefLength = atoi(i->Length.c_str());
                }

                if(!writer.Open(filename, header, references)) {
                    cerr << "Error opening BAM file to write." << endl;
                    exit(-1);
                }

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
            }
            break;
        default:
            cerr << "Unsupported output file format selected. Aborting." << endl;
            exit(-1);
            break;
    }
    
    if(isVerbose())
        cerr << "Wrote " << write_count << " reads to " << filename << endl;
    
    return 0;
}