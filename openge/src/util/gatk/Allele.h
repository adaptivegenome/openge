/*********************************************************************
 *
 * Allele.h: Port of GATK's Allele class.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 8 May 2012
 *
 *********************************************************************
 *
 * This file has been ported from GATK's implementation in Java, and
 * is released under the Virginia Tech Non-Commercial Purpose License.
 * A copy of this license has been provided in  the openge/ directory.
 * 
 * The original file, Allele.java, was released 
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

#ifndef OpenGE_Allele_h
#define OpenGE_Allele_h

#include <string>
#include <vector>

class Allele{
private:
    static const char * EMPTY_ALLELE_BASES;
    
    bool m_isRef;
    bool m_isNull;
    bool m_isNoCall;
    bool m_isSymbolic;
    
    std::string bases;
    
    const static Allele REF_A;
    const static Allele ALT_A;
    const static Allele REF_C;
    const static Allele ALT_C;
    const static Allele REF_G;
    const static Allele ALT_G;
    const static Allele REF_T;
    const static Allele ALT_T;
    const static Allele REF_N;
    const static Allele ALT_N;
    const static Allele REF_NULL;
    const static Allele ALT_NULL;
    
public:
    const static std::string NULL_ALLELE_STRING;
    const static std::string NO_CALL_STRING;
    const static Allele NO_CALL;
    
private:
    Allele(std::string bases, bool isRef);
    
public:
    static bool wouldBeNullAllele(std::string bases);
    static bool wouldBeNoCallAllele(std::string bases);
    static bool acceptableAlleleBases(std::string bases);
    static bool wouldBeSymbolicAllele(std::string bases);

    bool isNull() const             { return m_isNull; }
    bool isNonNull() const          { return ! isNull(); }
    
    bool isNoCall() const           { return m_isNoCall; }
    bool isCalled() const           { return ! isNoCall(); }
    
    bool isReference() const        { return m_isRef; }
    bool isNonReference() const     { return ! isReference(); }
    
    bool isSymbolic() const         { return m_isSymbolic; }
    
    std::string getBases();
    int length() const;
    bool operator==(const Allele & other) const;
};

#endif
