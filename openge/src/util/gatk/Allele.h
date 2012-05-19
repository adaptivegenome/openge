//
//  Allele.h
//  OpenGE
//
//  Created by Lee Baker on 5/8/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

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
