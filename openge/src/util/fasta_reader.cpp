/*********************************************************************
 *
 * FastaReader.cpp: A high performance FASTA reader.
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
#include <fstream>
#include <cassert>
#include <cstring>
#include <fcntl.h>
#include <errno.h>

#include "fasta_reader.h"
#include <sys/mman.h>

#include <algorithm>
using namespace std;
using BamTools::SamSequenceDictionary;

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

bool hasSuffix (std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

bool fileExists(const std::string & filename)
{
    return ifstream(filename.c_str()).good();
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
    
    string index_filename;
    if(fileExists(filename + ".fai"))
        index_filename = filename + ".fai"; //append .fai
    else if(hasSuffix(filename, "fa") && fileExists(filename + "i"))
        index_filename = filename + "i"; // .fa to .fai
    else if(hasSuffix(filename, "fasta") && fileExists(filename.substr(0, filename.size() - 3) + "i") )
        index_filename = filename.substr(0, filename.size() - 3) + "i"; // .fasta to .fai

    bool index_exists = true;
    FILE *file;
    if ((file = fopen(index_filename.c_str(), "r")) == NULL) {
        if (errno == ENOENT) {
            index_exists = false;
        } else {
            // Check for other errors too, like EACCES and EISDIR
        }
    } else {
        fclose(file);
    }
    
    if(!index_exists) { //create index
        ofstream ix_file(index_filename.c_str(),ios_base::out);
        ix_file << generateFastaIndex();
        cerr << "Generated FASTA index at " << index_filename << endl;
    } else
        readFastaIndex(index_filename);

    is_open = true;
    return true;
}

string FastaReader::getSubsequenceAt(string name, size_t start, size_t stop) const
{
    return readSequence(name, start, stop - start);
}

string FastaReader::readSequence(string name, size_t start, size_t length) const
{
    string read;
    
    map<string, fasta_sequence_t>::const_iterator i = sequences.find(name);
    
    if(i == sequences.end()) {
        cerr << "Sequence " << name << " not found in FASTA. Aborting." << endl;
        exit(-1);
    }
    
    const fasta_sequence_t & seq = i->second;
    
    if(start + length > seq.length) {
        cerr << "Requested FASTA read " << name << ":" << start << "-" << (start + length - 1) << " is beyond end of sequence(" << seq.length << "). Aborting." << endl;
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
    return writeFastaIndex();
}

string FastaReader::writeFastaIndex() const
{
    stringstream ss("");
    for(vector<fasta_sequence_t>::const_iterator i = ordered_sequences.begin(); i != ordered_sequences.end(); i++)
        ss << i->name << "\t" << i->length << "\t" << (size_t)(i->sequence_start - file_data) << "\t" << i->line_data_length << "\t" << i->line_length << endl;
    
    return ss.str();
}

bool FastaReader::readFastaIndex(string filename)
{
    ifstream file(filename.c_str());

    while(true) {
        string line;
        getline(file, line);
        if(file.fail())
            break;
        
        stringstream ss(line);
        fasta_sequence_t seq = {0};
        
        int64_t start_offset;
        ss >> seq.name >> seq.length >> start_offset >> seq.line_data_length >> seq.line_length;
        seq.sequence_start = &file_data[start_offset];
        sequences[seq.name] = seq;
        ordered_sequences.push_back(seq);
    }
    
    //check that FASTA entries are where we expect them
    if(ordered_sequences.size() == 0) {
        cerr << "Error- no sequences from FASTA " << filename << " were loaded. Quitting." << endl;
    } else {
        char * data_p = file_data;
        for(int seq_ctr = 0; seq_ctr < ordered_sequences.size(); seq_ctr++) {
            
            // check for the header of the next sequence
            if(*data_p != '>') {
                cerr << "Error: FASTA index does not match FASTA file. Make sure the FASTA and FAI are properly\nformatted. Aborting." << endl;
                exit(-1);
            }
            fasta_sequence_t & seq = ordered_sequences[seq_ctr];
            int full_lines = seq.length / seq.line_data_length;
            int last_line_length = seq.length % seq.line_data_length;
            if(last_line_length != 0)
                last_line_length++;
            
            // check for the beginning of the next sequence.
            while(*(++data_p) != '\n');
            
            if(data_p + 1 != seq.sequence_start) {
                cerr << "Error: FASTA index does not match FASTA file. Make sure the FASTA and FAI are properly\nformatted. Aborting." << endl;
                exit(-1);
            }

            if(seq_ctr + 1 < ordered_sequences.size()) {
                data_p = ordered_sequences[seq_ctr].sequence_start + full_lines * seq.line_length + last_line_length ;
            }
        }
    }
    
    return true;
}

void FastaReader::Close()
{
    //remove file mapping
    munmap(file_data, file_length);

    close(file);
    file = 0;

    is_open = false;
}
