/*********************************************************************
 *
 * SequenceUtil.h: Port of GATK's SequenceUtil.
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
 * The original file, SequenceUtil.java, was released 
 * under the following license:
 *
 * Copyright (c) 2009 The Broad Institute
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
