/*********************************************************************
 *
 * VariantContext.h: Port of GATK's VariantContext.
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
 * The original file, VariantContext.java, was released 
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

#ifndef OpenGE_VariantContext_h
#define OpenGE_VariantContext_h

#include <vector>
#include "Allele.h"
#include "GATKFeature.h"

class VariantContext : public GATKFeature{
public:
    typedef enum {
        NO_VARIATION,
        SNP,
        MNP,    // a multi-nucleotide polymorphism
        INDEL,
        SYMBOLIC,
        MIXED,
    } Type;

private:
    // a fast cached access point to the ref / alt alleles for biallelic case
    Allele * REF;
    
    // set to the alt allele when biallelic, otherwise == null
    Allele * ALT;
    
    /** A set of the alleles segregating in this context */
    std::vector<Allele> alleles;
    /** The type (cached for performance reasons) of this context */
    Type * type;
public:
    VariantContext()
    : REF(NULL)
    , ALT(NULL)
    , type(NULL)
    {}

    bool isSimpleInsertion();
    bool isSimpleDeletion();
    Allele getReference();
    bool isBiallelic();
    bool isComplexIndel();
    bool isSymbolic();
    bool isMNP();
    bool isIndel();
    
    int getNAlleles();
    std::vector<Allele> getAlternateAlleles();
    Allele getAlternateAllele(int i);
    Type getType();
    std::string getChr();
    int getStart();
    int getEnd();private:
    void determineType();
private:
    void determinePolymorphicType();
    static VariantContext::Type typeOfBiallelicVariant(const Allele & ref, const Allele & allele);
};

#endif
