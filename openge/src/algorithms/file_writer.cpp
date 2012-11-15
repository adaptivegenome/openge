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

#include <fstream>

#include "../util/sam_writer.h"
#include "../util/fastq_writer.h"
#include "../util/bam_serializer.h"
#include "../util/bgzf_output_stream.h"

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
    BamHeader header = getHeader();
    
    if(command_line_options.size() > 0) {
        BamProgramRecord pg;
        pg.id = string("openge");
        pg.version = string(OPENGE_VERSION_STRING);

        for(int i = 2; header.getPrograms().contains( pg.id); i++) {
            stringstream s;
            s << "openge-" << i;
            pg.id = s.str();
        }

        pg.commandLine = command_line_options;
        header.getPrograms().add(pg);
    }

    switch(file_format) {
        case FORMAT_SAM:
            {
                SamWriter writer;
                
                if(!writer.open(filename, header)) {
                    cerr << "Error opening BAM file to write." << endl;
                    exit(-1);
                }
                
                OGERead * al;
                
                while(true)
                {
                    al = getInputAlignment();
                    if(!al)
                        break;
                    
                    writer.write(*al);
                    putOutputAlignment(al);
                }
                
                writer.close();
                
            }
            break;
        case FORMAT_FASTQ:
            {
                FastqWriter writer;
                
                if(!writer.open(filename, header)) {
                    cerr << "Error opening FASTQ file to write." << endl;
                    exit(-1);
                }
                
                OGERead * al;
                
                while(true)
                {
                    al = getInputAlignment();
                    if(!al)
                        break;
                    
                    writer.write(*al);
                    putOutputAlignment(al);
                }
                
                writer.close();
                
            }
            break;
        case FORMAT_BAM:
            {
                BamSerializer<BgzfOutputStream> writer;

                writer.getOutputStream().setCompressionLevel(compression_level);

                if(!writer.open(filename, header)) {
                    cerr << "Error opening BAM file to write." << endl;
                    exit(-1);
                }

                OGERead * al;
                
                while(true)
                {
                    al = getInputAlignment();
                    if(!al)
                        break;

                    writer.write(*al);
                    putOutputAlignment(al);
                }
                
                writer.close();
            }
            break;
        case FORMAT_RAWBAM:
        {
            BamSerializer<ofstream> writer;
            
            if(!writer.open(filename, header)) {
                cerr << "Error opening RAWBAM file to write." << endl;
                exit(-1);
            }
            
            OGERead * al;
            
            while(true)
            {
                al = getInputAlignment();
                if(!al)
                    break;
                
                writer.write(*al);
                putOutputAlignment(al);
            }
            
            writer.close();
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