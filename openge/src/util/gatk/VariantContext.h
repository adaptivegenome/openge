//
//  VariantContext.h
//  OpenGE
//
//  Created by Lee Baker on 5/8/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#ifndef OpenGE_VariantContext_h
#define OpenGE_VariantContext_h

#include <vector>
#include "Allele.h"

class VariantContext {
public:
    typedef enum {
        NO_VARIATION,
        SNP,
        MNP,    // a multi-nucleotide polymorphism
        INDEL,
        SYMBOLIC,
        MIXED,
    } Type;
protected:   
    std::string contig;
    long start;
    long stop;
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
