//
//  Header.h
//  OpenGE
//
//  Created by Lee Baker on 5/6/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#ifndef OpenGE_GenomeLocParser_h
#define OpenGE_GenomeLocParser_h

#include "api/BamAlignment.h"
#include "api/BamAux.h"
#include "api/SamSequenceDictionary.h"
#include "GenomeLoc.h"

#include <string>

class GenomeLocParser {
private:
    BamTools::SamSequenceDictionary * contigInfo;
public:
    GenomeLocParser( BamTools::SamSequenceDictionary seqDict);
    bool contigIsInDictionary(const std::string contig) const;
    bool indexIsInDictionary(const int index) const;

    BamTools::SamSequence getContigInfo(std::string contig);
    int getContigIndex(std::string contig) const;
    std::string getContig(const int index) const;
    GenomeLoc parseGenomeLoc(const std::string str) const;

    GenomeLoc createGenomeLoc(const BamTools::BamAlignment & read) const;

    GenomeLoc createGenomeLoc(const std::string contig, const int start, const int stop) const;    
    GenomeLoc createGenomeLoc(const std::string & contig, const int start, const int stop, bool mustBeOnReference) const;
    GenomeLoc createGenomeLoc(const std::string & contig, int index, const int start, const int stop) const;
    GenomeLoc createGenomeLoc(const std::string & contig, int index, const int start, const int stop, bool mustBeOnReference) const;
protected:
    int getContigIndexWithoutException(std::string contig) const;
private:
    bool validateGenomeLoc(const std::string contig, const int contigIndex, const int start, const int stop, const bool mustBeOnReference, const const bool exceptOnError) const;
};

#endif
