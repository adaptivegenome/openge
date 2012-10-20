/*********************************************************************
 *
 * GenomeLocParser.h: Port of GATK's GenomeLocParser.
 * Open GenomeLocParser Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 6 May 2012
 *
 *********************************************************************
 *
 * This file has been ported from GATK's implementation in Java, and
 * is released under the Virginia Tech Non-Commercial Purpose License.
 * A copy of this license has been provided in  the openge/ directory.
 * 
 * The original file, GenomeLocParser.java, was released 
 * under the following license:
 *
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef OpenGE_GenomeLocParser_h
#define OpenGE_GenomeLocParser_h

#include "api/BamAlignment.h"
#include "api/BamAux.h"
#include "api/SamSequenceDictionary.h"
#include "GenomeLoc.h"
#include "../oge_read.h"

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

    GenomeLoc createGenomeLoc(const OGERead & read) const;

    GenomeLoc createGenomeLoc(const std::string contig, const int start, const int stop) const;    
    GenomeLoc createGenomeLoc(const std::string & contig, const int start, const int stop, bool mustBeOnReference) const;
    GenomeLoc createGenomeLoc(const std::string & contig, int index, const int start, const int stop) const;
    GenomeLoc createGenomeLoc(const std::string & contig, int index, const int start, const int stop, bool mustBeOnReference) const;
protected:
    int getContigIndexWithoutException(std::string contig) const;
private:
    bool validateGenomeLoc(const std::string contig, const int contigIndex, const int start, const int stop, const bool mustBeOnReference, const bool exceptOnError) const;
};

#endif
