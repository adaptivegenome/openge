/*********************************************************************
 *
 * AlignmentUtils.cpp: Port of GATK's AlignmentUtils.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 9 May 2012
 *
 *********************************************************************
 *
 * This file has been ported from GATK's implementation in Java, and
 * is released under the Virginia Tech Non-Commercial Purpose License.
 * A copy of this license has been provided in  the openge/ directory.
 * 
 * The original file, AlignmentUtils.java, was released 
 * under the following license:
 *
 * Copyright (c) 2010 The Broad Institute. 
 * Ported to C++ by Lee C. Baker, Virginia Bioinformatics Institute
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

#ifndef OpenGE_AlignmentUtils_h
#define OpenGE_AlignmentUtils_h

#include "../oge_read.h"
#include <string>
#include <vector>

class AlignmentUtils
{
public:
    class MismatchCount {
    public:
        int numMismatches;
        long mismatchQualities;
        MismatchCount() 
        : numMismatches(0)
        , mismatchQualities(0)
        {}
    };
    
    static long mismatchingQualities(const OGERead * r, std::string refSeq, int refIndex);
    static MismatchCount getMismatchCount(const OGERead * r, std::string refSeq, int refIndex);
    static MismatchCount getMismatchCount(const OGERead * r, std::string refSeq, int refIndex, int startOnRead, int nReadBases) ;
    
    static std::vector<CigarOp> leftAlignIndel( std::vector<CigarOp> cigar, const std::string refSeq, const std::string readSeq, const int refIndex, const int readIndex);
    
private:
    static std::vector<CigarOp> moveCigarLeft(const std::vector<CigarOp> & cigar, int indexOfIndel);
    static bool cigarHasZeroSizeElement(const std::vector<CigarOp> & c);
    static std::vector<CigarOp> cleanUpCigar(const std::vector<CigarOp> & c);
    static std::string createIndelString(const std::vector<CigarOp> & cigar, const int indexOfIndel, const std::string refSeq, const std::string readSeq, int refIndex, int readIndex);
};

#endif
