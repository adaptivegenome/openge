//
//  FastaReader.h
//  OpenGE
//
//  Created by Lee Baker on 5/7/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#ifndef OpenGE_FastaReader_h
#define OpenGE_FastaReader_h

#include <string>
#include <map>
#include <vector>
#include "api/SamSequenceDictionary.h"

class FastaReader {
public:
protected:
    typedef struct {
        char * sequence_start;
        std::string name;
        int line_length;
        int line_data_length;
        size_t length;
    } fasta_sequence_t;
    std::map<std::string, fasta_sequence_t> sequences;
    std::vector<fasta_sequence_t> ordered_sequences;
    bool is_open;
    char * file_data;
    size_t file_length;
    int file;
public:
    FastaReader();
    ~FastaReader();
    BamTools::SamSequenceDictionary getSequenceDictionary();
    bool Open(const std::string filename);
    std::string readSequence(std::string name, size_t start, size_t length);
    std::string getSubsequenceAt(std::string name, size_t start, size_t stop);
    std::string generateFastaIndex();
    void Close();
protected:
};

#endif
