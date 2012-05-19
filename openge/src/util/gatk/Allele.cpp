//
//  Allele.cpp
//  OpenGE
//
//  Created by Lee Baker on 5/8/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#include "Allele.h"
#include <iostream>
#include <string>
using namespace std;
/*
package org.broadinstitute.sting.utils.variantcontext;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
 */

/**
 * Immutable representation of an allele
 *
 * Types of alleles:
 *
 * Ref: a t C g a // C is the reference base
 *
 *    : a t G g a // C base is a G in some individuals
 *
 *    : a t - g a // C base is deleted w.r.t. the reference
 *
 *    : a t CAg a // A base is inserted w.r.t. the reference sequence
 *
 * In these cases, where are the alleles?
 *
 * SNP polymorphism of C/G  -> { C , G } -> C is the reference allele
 * 1 base deletion of C     -> { C , - } -> C is the reference allele
 * 1 base insertion of A    -> { - ; A } -> Null is the reference allele
 *
 * Suppose I see a the following in the population:
 *
 * Ref: a t C g a // C is the reference base
 *    : a t G g a // C base is a G in some individuals
 *    : a t - g a // C base is deleted w.r.t. the reference
 *
 * How do I represent this?  There are three segregating alleles:
 *
 *  { C , G , - }
 *
 * Now suppose I have this more complex example:
 *
 * Ref: a t C g a // C is the reference base
 *    : a t - g a
 *    : a t - - a
 *    : a t CAg a
 *
 * There are actually four segregating alleles:
 *
 *   { C g , - g, - -, and CAg } over bases 2-4
 *
 * However, the molecular equivalence explicitly listed above is usually discarded, so the actual
 * segregating alleles are:
 *
 *   { C g, g, -, C a g }
 *
 * Critically, it should be possible to apply an allele to a reference sequence to create the
 * correct haplotype sequence:
 *
 * Allele + reference => haplotype
 *
 * For convenience, we are going to create Alleles where the GenomeLoc of the allele is stored outside of the
 * Allele object itself.  So there's an idea of an A/C polymorphism independent of it's surrounding context.
 *
 * Given list of alleles it's possible to determine the "type" of the variation
 *
 *      A / C @ loc => SNP with
 *      - / A => INDEL
 *
 * If you know where allele is the reference, you can determine whether the variant is an insertion or deletion.
 *
 * Alelle also supports is concept of a NO_CALL allele.  This Allele represents a haplotype that couldn't be
 * determined. This is usually represented by a '.' allele.
 *
 * Note that Alleles store all bases as bytes, in **UPPER CASE**.  So 'atc' == 'ATC' from the perspective of an
 * Allele.
 
 * @author ebanks, depristo
 */

const char * Allele::EMPTY_ALLELE_BASES = "";

const string Allele::NULL_ALLELE_STRING("-");
const string Allele::NO_CALL_STRING(".");

const Allele Allele::REF_A("A", true);
const Allele Allele::ALT_A("A", false);
const Allele Allele::REF_C("C", true);
const Allele Allele::ALT_C("C", false);
const Allele Allele::REF_G("G", true);
const Allele Allele::ALT_G("G", false);
const Allele Allele::REF_T("T", true);
const Allele Allele::ALT_T("T", false);
const Allele Allele::REF_N("N", true);
const Allele Allele::ALT_N("N", false);
const Allele Allele::REF_NULL(NULL_ALLELE_STRING, true);
const Allele Allele::ALT_NULL(NULL_ALLELE_STRING, false);
const Allele Allele::NO_CALL(NO_CALL_STRING, false);

#if 0
class Allele {
    private static final byte[] EMPTY_ALLELE_BASES = new byte[0];
    
    private boolean isRef = false;
    private boolean isNull = false;
    private boolean isNoCall = false;
    private boolean isSymbolic = false;
    
    private byte[] bases = null;
    
    public final static String NULL_ALLELE_STRING = "-";
    public final static String NO_CALL_STRING = ".";
    /** A generic static NO_CALL allele for use */
    
#endif
    // no public way to create an allele
    Allele::Allele(string bases, bool isRef)
    : m_isRef(false)
    , m_isNull(false)
    , m_isNoCall(false)
    , m_isSymbolic(false)
    {
        // standardize our representation of null allele and bases
        if ( wouldBeNullAllele(bases) ) {
            bases = EMPTY_ALLELE_BASES;
            m_isNull = true;
        } else if ( wouldBeNoCallAllele(bases) ) {
            bases = EMPTY_ALLELE_BASES;
            m_isNoCall = true;
            if ( isRef ) {
                cerr << "Cannot tag a NoCall allele as the reference allele" << endl;
                abort();
            }
        } else if ( wouldBeSymbolicAllele(bases) ) {
            m_isSymbolic = true;
            if ( isRef ) {
                cerr << "Cannot tag a symbolic allele as the reference allele" << endl;
                abort();
            }
        }
        //        else
        //            bases = new String(bases).toUpperCase().getBytes(); // todo -- slow performance
        
        this->m_isRef = isRef;
        this->bases = bases;
        
        if ( ! acceptableAlleleBases(bases) ) {
            cerr << "Unexpected base in allele bases \'" << bases << "\'" << endl;
            abort();
        }
    }
#if 0
    
    
    private final static Allele REF_A = new Allele("A", true);
    private final static Allele ALT_A = new Allele("A", false);
    private final static Allele REF_C = new Allele("C", true);
    private final static Allele ALT_C = new Allele("C", false);
    private final static Allele REF_G = new Allele("G", true);
    private final static Allele ALT_G = new Allele("G", false);
    private final static Allele REF_T = new Allele("T", true);
    private final static Allele ALT_T = new Allele("T", false);
    private final static Allele REF_N = new Allele("N", true);
    private final static Allele ALT_N = new Allele("N", false);
    private final static Allele REF_NULL = new Allele(NULL_ALLELE_STRING, true);
    private final static Allele ALT_NULL = new Allele(NULL_ALLELE_STRING, false);
    public final static Allele NO_CALL = new Allele(NO_CALL_STRING, false);
    
    // ---------------------------------------------------------------------------------------------------------
    //
    // creation routines
    //
    // ---------------------------------------------------------------------------------------------------------
    
    /**
     * Create a new Allele that includes bases and if tagged as the reference allele if isRef == true.  If bases
     * == '-', a Null allele is created.  If bases ==  '.', a no call Allele is created.
     *
     * @param bases the DNA sequence of this variation, '-', of '.'
     * @param isRef should we make this a reference allele?
     * @throws IllegalArgumentException if bases contains illegal characters or is otherwise malformated
     */
    public static Allele create(byte[] bases, boolean isRef) {
        if ( bases == null )
            throw new IllegalArgumentException("create: the Allele base string cannot be null; use new Allele() or new Allele(\"\") to create a Null allele");
        
        if ( bases.length == 1 ) {
            // optimization to return a static constant Allele for each single base object
            switch (bases[0]) {
                case '.':
                    if ( isRef ) throw new IllegalArgumentException("Cannot tag a NoCall allele as the reference allele");
                    return NO_CALL;
                case '-': return isRef ? REF_NULL : ALT_NULL;
                case 'A': case 'a' : return isRef ? REF_A : ALT_A;
                case 'C': case 'c' : return isRef ? REF_C : ALT_C;
                case 'G': case 'g' : return isRef ? REF_G : ALT_G;
                case 'T': case 't' : return isRef ? REF_T : ALT_T;
                case 'N': case 'n' : return isRef ? REF_N : ALT_N;
                default: throw new IllegalArgumentException("Illegal base: " + (char)bases[0]);
            }
        } else {
            return new Allele(bases, isRef);
        }
    }
    
    public static Allele create(byte base, boolean isRef) {
        //    public Allele(byte base, boolean isRef) {
        return create( new byte[]{ base }, isRef);
    }
    
    public static Allele create(byte base) {
        return create( base, false );
    }
    
    public static Allele extend(Allele left, byte[] right) {
        if (left.isSymbolic())
            throw new IllegalArgumentException("Cannot extend a symbolic allele");
        byte[] bases = null;
        if ( left.length() == 0 )
            bases = right;
        else {
            bases = new byte[left.length() + right.length];
            System.arraycopy(left.getBases(), 0, bases, 0, left.length());
            System.arraycopy(right, 0, bases, left.length(), right.length);
        }
        
        return create(bases, left.isReference());
    }
#endif
    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent the null allele
     */
    bool Allele::wouldBeNullAllele(string bases) {
        return (bases.size() == 1 && bases[0] == '-') || bases.size() == 0;
    }
    
    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent the NO_CALL allele
     */
    bool Allele::wouldBeNoCallAllele(string bases) {
        return bases.size() == 1 && bases[0] == '.';
    }
    
    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent a symbolic allele
     */
    bool Allele::wouldBeSymbolicAllele(string bases) {
        if ( bases.size() <= 2 )
            return false;
        else {
            const string strBases = bases;
            return (bases[0] == '<' && bases[bases.size()-1] == '>') ||
            (string::npos != strBases.find("[") || string::npos != strBases.find("]"));
        }
    }
    
    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent the well formatted allele
     */
    bool Allele::acceptableAlleleBases(string bases) {
        if ( wouldBeNullAllele(bases) || wouldBeNoCallAllele(bases) || wouldBeSymbolicAllele(bases) )
            return true;
        
        for ( int i = 0; i < bases.size(); i++ ) {
            switch (bases[i]) {
                case 'A': case 'C': case 'G': case 'T': case 'N' : case 'a': case 'c': case 'g': case 't': case 'n' :
                    break;
                default:
                    return false;
            }
        }
        
        return true;
    }
#if 0
    
    /**
     * @see Allele(byte[], boolean)
     *
     * @param bases  bases representing an allele
     * @param isRef  is this the reference allele?
     */
    public static Allele create(String bases, boolean isRef) {
        //public Allele(String bases, boolean isRef) {
        return create(bases.getBytes(), isRef);
    }
    
    
    /**
     * Creates a non-Ref allele.  @see Allele(byte[], boolean) for full information
     *
     * @param bases  bases representing an allele
     */
    public static Allele create(String bases) {
        return create(bases, false);
    }
    
    /**
     * Creates a non-Ref allele.  @see Allele(byte[], boolean) for full information
     *
     * @param bases  bases representing an allele
     */
    public static Allele create(byte[] bases) {
        return create(bases, false);
        //this(bases, false);
    }
    
    // ---------------------------------------------------------------------------------------------------------
    //
    // accessor routines
    //
    // ---------------------------------------------------------------------------------------------------------
    
    // Returns a nice string representation of this object
    public String toString() {
        return (isNull() ? NULL_ALLELE_STRING : ( isNoCall() ? NO_CALL_STRING : getDisplayString() )) + (isReference() ? "*" : "");
    }
#endif
    
    /**
     * Return the DNA bases segregating in this allele.  Note this isn't reference polarized,
     * so the Null allele is represented by a vector of length 0
     *
     * @return the segregating bases
     */
    string Allele::getBases() { return m_isSymbolic ? EMPTY_ALLELE_BASES : bases; }
#if 0
    /**
     * Return the DNA bases segregating in this allele in String format.
     * This is useful, because toString() adds a '*' to reference alleles and getBases() returns garbage when you call toString() on it.
     *
     * @return the segregating bases
     */
    public String getBaseString() { return new String(getBases()); }
    
    /**
     * Return the printed representation of this allele.
     * Same as getBaseString(), except for symbolic alleles.
     * For symbolic alleles, the base string is empty while the display string contains <TAG>.
     *
     * @return the allele string representation
     */
    public String getDisplayString() { return new String(bases); }
    
    /**
     * @param other  the other allele
     *
     * @return true if these alleles are equal
     */
    public boolean equals(Object other) {
        return ( ! (other instanceof Allele) ? false : equals((Allele)other, false) );
    }
    
    /**
     * @return hash code
     */
    public int hashCode() {
        int hash = 1;
        for (int i = 0; i < bases.length; i++)
            hash += (i+1) * bases[i];
        return hash;
    }
#endif
    
    /**
     * Returns true if this and other are equal.  If ignoreRefState is true, then doesn't require both alleles has the
     * same ref tag
     *
     * @param other            allele to compare to
     * @param ignoreRefState   if true, ignore ref state in comparison
     * @return true if this and other are equal
     */
    bool Allele::operator==(const Allele & other) const {
        bool ignoreRefState = false;
        return this == &other 
          || ((m_isRef == other.m_isRef || ignoreRefState)
          && m_isNull == other.m_isNull 
          && m_isNoCall == other.m_isNoCall && bases == other.bases);
    }
#if 0
    
    /**
     * @param test  bases to test against
     *
     * @return  true if this Alelle contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
     */
    public boolean basesMatch(byte[] test) { return !isSymbolic && (bases == test || Arrays.equals(bases, test)); }
    
    /**
     * @param test  bases to test against
     *
     * @return  true if this Alelle contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
     */
    public boolean basesMatch(String test) { return basesMatch(test.toUpperCase().getBytes()); }
    
    /**
     * @param test  allele to test against
     *
     * @return  true if this Alelle contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
     */
    public boolean basesMatch(Allele test) { return basesMatch(test.getBases()); }
    
#endif
    /**
     * @return the length of this allele.  Null and NO_CALL alleles have 0 length.
     */
    int Allele::length() const {
        return m_isSymbolic ? 0 : bases.size();
    }
#if 0
    
    // ---------------------------------------------------------------------------------------------------------
    //
    // useful static functions
    //
    // ---------------------------------------------------------------------------------------------------------
    
    public static Allele getMatchingAllele(Collection<Allele> allAlleles, String alleleBases) {
        return getMatchingAllele(allAlleles, alleleBases.getBytes());
    }
    
    public static Allele getMatchingAllele(Collection<Allele> allAlleles, byte[] alleleBases) {
        for ( Allele a : allAlleles ) {
            if ( a.basesMatch(alleleBases) ) {
                return a;
            }
        }
        
        if ( wouldBeNoCallAllele(alleleBases) )
            return NO_CALL;
        else
            return null;    // couldn't find anything
    }
    
    public static List<Allele> resolveAlleles(List<Allele> possibleAlleles, List<String> alleleStrings) {
        List<Allele> myAlleles = new ArrayList<Allele>(alleleStrings.size());
        
        for ( String alleleString : alleleStrings ) {
            Allele allele = getMatchingAllele(possibleAlleles, alleleString);
            
            if ( allele == null ) {
                if ( Allele.wouldBeNoCallAllele(alleleString.getBytes()) ) {
                    allele = create(alleleString);
                } else {
                    throw new IllegalArgumentException("Allele " + alleleString + " not present in the list of alleles " + possibleAlleles);
                }
            }
            
            myAlleles.add(allele);
        }
        
        return myAlleles;
    }
    
    public int compareTo(Allele other) {
        if ( isReference() && other.isNonReference() )
            return -1;
        else if ( isNonReference() && other.isReference() ) 
            return 1;
        else
            return getBaseString().compareTo(other.getBaseString()); // todo -- potential performance issue
    }
    
    public static boolean oneIsPrefixOfOther(Allele a1, Allele a2) {
        if ( a1.isNull() || a2.isNull() )
            return true;
        
        if ( a2.length() >= a1.length() )
            return firstIsPrefixOfSecond(a1, a2);
        else
            return firstIsPrefixOfSecond(a2, a1);
    }
    
    private static boolean firstIsPrefixOfSecond(Allele a1, Allele a2) {
        String a1String = a1.getBaseString();
        return a2.getBaseString().substring(0, a1String.length()).equals(a1String);
    }
}
#endif
