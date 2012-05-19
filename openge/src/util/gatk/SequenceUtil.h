//
//  SequenceUtil.h
//  OpenGE
//
//  Created by Lee Baker on 5/9/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#ifndef OpenGE_SequenceUtil_h
#define OpenGE_SequenceUtil_h

#include <string>
#include "api/BamAlignment.h"

class SequenceUtil {
public:
    static bool basesEqual(char lhs, char rhs);
    static bool bisulfiteBasesEqual(const bool negativeStrand, const char read, const char reference);
    static int countMismatches(const BamTools::BamAlignment * read, const std::string referenceBases, const int referenceOffset = 0, const bool bisulfiteSequence = false);
    static int calculateSamNmTag(const BamTools::BamAlignment * read, const std::string referenceBases, const int referenceOffset = 0, const bool bisulfiteSequence = false);
    static int sumQualitiesOfMismatches(const BamTools::BamAlignment * read, const std::string referenceBases, const int referenceOffset = 0, const bool bisulfiteSequence = false);
};

#endif
