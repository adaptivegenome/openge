/*********************************************************************
 *
 * command_help.cpp: A high performance FASTA reader.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 7 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <string>
#include <sstream>
#include <cassert>
#include <cstring>
#include <fcntl.h>

#include "fasta_reader.h"
#include <sys/mman.h>

#include <algorithm>
using namespace std;
using namespace BamTools;

FastaReader::FastaReader()
: is_open(false)
{
    
}

FastaReader::~FastaReader()
{
    if(is_open)
        Close();
}

SamSequenceDictionary FastaReader::getSequenceDictionary()
{
    SamSequenceDictionary d;
    
    for(vector<fasta_sequence_t>::const_iterator i = ordered_sequences.begin(); i != ordered_sequences.end(); i++)
        d.Add(i->name, i->length);

    return d;
}

bool FastaReader::Open(const string filename)
{
    //first determine file length
    file = open(filename.c_str(), O_RDONLY);
    if(-1 == file) {
        cerr << "Could not open file: " << filename << endl;
        perror("Reason:");
        return false;
    }
    
    const int SEEK_BEGIN = 0;

    file_length = lseek(file, 0, SEEK_END);
    lseek(file, 0, SEEK_BEGIN);

    file_data = (char *) mmap(NULL, file_length, PROT_READ, MAP_SHARED, file, 0);

    if(file_data == MAP_FAILED) {
        cerr << "Couldn't mmap file " << filename << endl;
        perror("Reason:");
        return false;
    }
    
    // For now, don't support index files- always generate our own index, then load it.
    // Unfortunately this is slow, but we can revisit this later.
    cerr << "Generating FASTA index (using existing index is not implemented yet)...";
    string index = generateFastaIndex();
    cerr << "Index looks like this: " << endl << index << endl;

    is_open = true;
    cerr << "done." << endl;
    return true;
}

string FastaReader::getSubsequenceAt(string name, size_t start, size_t stop)
{
    return readSequence(name, start, stop - start);
}

string FastaReader::readSequence(string name, size_t start, size_t length)
{
    string read;
    
    map<string, fasta_sequence_t>::iterator i = sequences.find(name);
    
    if(i == sequences.end()) {
        cerr << "Sequence " << name << " not found in FASTA. Aborting." << endl;
        exit(-1);
    }
    
    fasta_sequence_t & seq = i->second;
    
    if(start + length >= seq.length) {
        cerr << "Requested FASTA read is beyond end of sequence. Aborting." << endl;
        assert(0);
        exit(-1);
    }
    
    size_t lines_to_skip = start / seq.line_data_length;

    char * data_p = seq.sequence_start + lines_to_skip * seq.line_length;
    size_t data_ix = lines_to_skip * seq.line_data_length;
    
    // first line to be read doesn't necessarily start from the beginning of the line,
    // or end at the end of the line.
    {
        char * read_start = data_p + (start - data_ix);
        size_t len_minus_start = seq.line_data_length - (start - data_ix);
        int read_len = min(length, len_minus_start);
        
        read.append(read_start, read_len);
        
        if(read.find('\n') != read.npos)
            assert(0);
        
        data_p += seq.line_length;
        data_ix += seq.line_data_length;
    }
    
    
    //Read lines, starting from the beginning of the line.
    while(data_ix < start + length) {
        
        int read_len = min((size_t)seq.line_data_length, start+length - data_ix);

        read.append(data_p, read_len);
        
        if(read.find('\n') != read.npos)
            assert(0);
        
        data_p += seq.line_length;
        data_ix += seq.line_data_length;
    }
    
    if(read.find('\n') != read.npos)
        assert(0);

    return read;
}

// Generates an index from the file. Can be used to create indexes. Maybe we should
// create and save indexes automatically? TODO LCB
string FastaReader::generateFastaIndex()
{
    char * p = file_data;
    while(p)
    {
        //find beginning of sequence.
        p = (char *) memchr(p, '>', file_length - (p - file_data));
        if(p == NULL)
            break;
        p++;
        
        fasta_sequence_t sequence;

        //parse out name
        char * name_whitespace = p + strcspn(p, " \t\n");
        sequence.name = string(p, name_whitespace);

        //find beginning of next line- start of sequence
        p = (char *) memchr(p, '\n', file_length - (p - file_data));
        sequence.sequence_start = ++p;
        char * line_start = p;

        //find end of data on this line
        char * line_data_end = line_start + strspn(line_start, "abcdefghiklmnopqrstuvwyzxABCDEFGHIKLMNOPQRSTUVWYZX*-");
        sequence.line_data_length = line_data_end - line_start;
        
        //find end of line (data stride)
        char * line_end = (char *) memchr(p, '\n', file_length - (p - file_data));
        sequence.line_length = line_end - line_start + 1;

        p = line_end + 1;

        int data_lines = 0;
        char * line = p;

        //find end of sequence by searching for either:
        // 1. a line that ends before we think it should by:
        //  a. end of file
        //  b. no newline where we expect one
        // 2. or a line starting with > where we expect data- ie immediately after the newline.
        
        while(true) {   //stay in this while loop until we are for sure on the last line of the sequence
            data_lines++;
            char * expected_line_end = line + sequence.line_length-1;
            if(expected_line_end >= file_data + file_length)  //if file ends, this is the last line
                break;
            if(expected_line_end[0] != '\n')  //no newline where expected
                break;
            if(expected_line_end[1] == '>') //next line starts a new sequence.
                break;

            line = expected_line_end + 1;
        }
        
        //so now 'line' is for sure the last line of data in this sequence
        line_end = line + strspn(line, "abcdefghiklmnopqrstuvwyzxABCDEFGHIKLMNOPQRSTUVWYZX*-");
        size_t data_on_this_line = line_end - line;
        sequence.length = data_on_this_line + data_lines * sequence.line_data_length;
        
        sequences[sequence.name] = sequence;
        ordered_sequences.push_back(sequence);
    }
    
    //now, generate index data.
    stringstream ss("");
    for(vector<fasta_sequence_t>::const_iterator i = ordered_sequences.begin(); i != ordered_sequences.end(); i++)
        ss << i->name << "\t" << i->length << "\t" << (size_t)(i->sequence_start - file_data) << "\t" << i->line_data_length << "\t" << i->line_length << endl;
    
    return ss.str();
}

void FastaReader::Close()
{
    //remove file mapping
    munmap(file_data, file_length);

    close(file);
    file = 0;

    is_open = false;
}
