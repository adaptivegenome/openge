//
//  AlignmentUtils.h
//  OpenGE
//
//  Created by Lee Baker on 5/9/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#ifndef OpenGE_AlignmentUtils_h
#define OpenGE_AlignmentUtils_h

#include "api/BamAlignment.h"
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
    
    static long mismatchingQualities(const BamTools::BamAlignment * r, std::string refSeq, int refIndex);    
    static MismatchCount getMismatchCount(const BamTools::BamAlignment * r, std::string refSeq, int refIndex);
    static MismatchCount getMismatchCount(const BamTools::BamAlignment * r, std::string refSeq, int refIndex, int startOnRead, int nReadBases) ;
    
    static std::vector<BamTools::CigarOp> leftAlignIndel( std::vector<BamTools::CigarOp> cigar, const std::string refSeq, const std::string readSeq, const int refIndex, const int readIndex);
    
private:
    static std::vector<BamTools::CigarOp> moveCigarLeft(const std::vector<BamTools::CigarOp> & cigar, int indexOfIndel);
    static bool cigarHasZeroSizeElement(const std::vector<BamTools::CigarOp> & c);
    static std::vector<BamTools::CigarOp> cleanUpCigar(const std::vector<BamTools::CigarOp> & c);
    static std::string * createIndelString(const std::vector<BamTools::CigarOp> & cigar, const int indexOfIndel, const std::string refSeq, const std::string readSeq, int refIndex, int readIndex);
};

#endif
