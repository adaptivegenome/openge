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
    bool open(const std::string filename);
    size_t getSequenceLength(const std::string & name) const;
    std::string readSequence(const std::string & name, size_t start, size_t length) const;
    std::string getSubsequenceAt(const std::string & name, size_t start, size_t stop) const;
    std::string generateFastaIndex();
    bool readFastaIndex(std::string filename);
    std::string writeFastaIndex() const;
    void close();
protected:
};

#endif
