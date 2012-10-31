/*********************************************************************
 *
 * GenomeLoc.cpp: Port of GATK's GenomeLoc.
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
 * The original file, GenomeLoc.java, was released 
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

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Mar 2, 2009
 * Time: 8:50:11 AM
 *
 * Genome location representation.  It is *** 1 *** based closed.  Note that GenomeLocs start and stop values
 * can be any positive or negative number, by design.  Bound validation is a feature of the GenomeLocParser,
 * and not a fundamental constraint of the GenomeLoc
 */
#include "GenomeLoc.h"

#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cstdio>
using namespace std;

// --------------------------------------------------------------------------------------------------------------
//
// constructors
//
// --------------------------------------------------------------------------------------------------------------

//@Requires({
//    "contig != null",
//    "contigIndex >= 0", // I believe we aren't allowed to create GenomeLocs without a valid contigIndex
//    "start <= stop"})
GenomeLoc::GenomeLoc( const string contig, int contigIndex, int start, int stop )
    : contigName (contig)
    , contigIndex(contigIndex)
    , start(start)
    , stop(stop)
    {}

/** Unsafe constructor for special constant genome locs */
GenomeLoc::GenomeLoc( const string contig )
    : contigName(contig)
    , contigIndex(-1)
    , start(0)
    , stop(0)
{}


//@Ensures("result != null")
// toString method - 
string GenomeLoc::toString() const  {
    char start_str [33], stop_str[33];
    

    if(GenomeLoc::isUnmapped(*this)) return "unmapped";
    if ( throughEndOfContigP() && atBeginningOfContigP() )
        return getContig();
    else if ( throughEndOfContigP() || getStart() == getStop() ) {
        sprintf(start_str, "%d", getStart()+1);
        return string("") + getContig() + string(":") + start_str;
    }
    else {
        sprintf(start_str, "%d", getStart()+1);
        sprintf(stop_str, "%d", getStop()+1);
        return string("") + getContig() + string(":") + start_str + string("-") + stop_str;
    }
}

/**
 * Returns a new GenomeLoc that represents the entire span of this and that.  Requires that
 * this and that GenomeLoc are contiguous and both mapped
 */
//@Requires({
//    "that != null",
//    "isUnmapped(this) == isUnmapped(that)"})
//@Ensures("result != null")
GenomeLoc GenomeLoc::merge( const GenomeLoc & that ) { //throws ReviewedStingException {
    if(isUnmapped(*this) || isUnmapped(that)) {
        if(! isUnmapped(*this) || !isUnmapped(that))
            cerr << "Tried to merge a mapped and an unmapped genome loc" << endl;
            return UNMAPPED;
    }
    
    assert(contiguousP(that));
        // throw new ReviewedStingException("The two genome loc's need to be contigous");

    return GenomeLoc(getContig(), contigIndex,
                         min( getStart(), that.getStart() ),
                         max( getStop(), that.getStop()) );
}

/**
 * Returns a new GenomeLoc that represents the region between the endpoints of this and that. Requires that
 * this and that GenomeLoc are both mapped.
 */
//@Requires({"that != null", "isUnmapped(this) == isUnmapped(that)"})
//@Ensures("result != null")
GenomeLoc GenomeLoc::endpointSpan(const GenomeLoc & that) {// throws ReviewedStingException {
    if(GenomeLoc::isUnmapped(*this) || GenomeLoc::isUnmapped(that)) {
        cerr << "Cannot get endpoint span for unmerged genome locs" << endl;
        assert(0);
    }
    
    if ( getContig() != that.getContig() ) {
        cerr << "Cannot get endpoint span for genome locs on different contigs" << endl;
        assert(0);
    }
    
    return GenomeLoc(getContig(),contigIndex,min(getStart(),that.getStart()),max(getStop(),that.getStop()));
}

/**
 * Splits the contig into to regions: [start,split point) and [split point, end].
 * @param splitPoint The point at which to split the contig.  Must be contained in the given interval.
 * @return A two element array consisting of the genome loc before the split and the one after.
 */
pair<GenomeLoc, GenomeLoc> GenomeLoc::split( int splitPoint) const {
    if(splitPoint < getStart() || splitPoint > getStop()) {
        cerr << "Unable to split contig " << this->toString() << " at split point " << splitPoint << "; split point is not contained in region." << endl;
        assert(0);
    }
    return pair<GenomeLoc, GenomeLoc> ( GenomeLoc(getContig(),contigIndex,getStart(),splitPoint-1), GenomeLoc(getContig(),contigIndex,splitPoint,getStop()));
}

//@Requires("that != null")
//@Ensures("result != null")
GenomeLoc GenomeLoc::intersect( const GenomeLoc & that ) const { //throws ReviewedStingException {
    if(isUnmapped(*this) || isUnmapped(that)) {
        if(! isUnmapped(*this) || !isUnmapped(that)) {
            cerr << "Tried to intersect a mapped and an unmapped genome loc" << endl;
            assert(0);
        }
        return UNMAPPED;
    }
    
    if (!(this->overlapsP(that))) {
        cerr << "GenomeLoc::intersect(): The two genome loc's need to overlap" << endl;
        assert(0);
    }
    
    return GenomeLoc(getContig(), contigIndex,
                         max(getStart(), that.getStart()),
                         min( getStop(), that.getStop()) );
}

//@Requires("that != null")
vector<GenomeLoc> GenomeLoc::subtract( const GenomeLoc & that ) {
    vector<GenomeLoc> ret;
    if(isUnmapped(*this) || isUnmapped(that)) {
        if(! isUnmapped(*this) || !isUnmapped(that)) {
            cerr << "Tried to intersect a mapped and an unmapped genome loc" << endl;
            assert(0);
        }
        ret.push_back(UNMAPPED);
        return ret;
    }
    
    if (!(overlapsP(that))) {
        cerr << "GenomeLoc::minus(): The two genome loc's need to overlap" << endl;
        assert(0);
    }
    
    if (*this == that) {
        return ret;
    } else if (containsP(that)) {
        /**
         * we have to create two new region, one for the before part, one for the after
         * The old region:
         * |----------------- old region (g) -------------|
         *        |----- to delete (e) ------|
         *
         * product (two new regions):
         * |------|  + |--------|
         *
         */
        int afterStop = getStop(), afterStart = that.getStop() + 1;
        int beforeStop = that.getStart() - 1, beforeStart = getStart();
        if (afterStop - afterStart >= 0) {
            GenomeLoc after(getContig(), getContigIndex(), afterStart, afterStop);
            ret.push_back(after);
        }
        if (beforeStop - beforeStart >= 0) {
            GenomeLoc before(getContig(), getContigIndex(), beforeStart, beforeStop);
            ret.push_back(before);
        }
        
        return ret;
    } else if (that.containsP(*this)) {
        /**
         * e completely contains g, delete g, but keep looking, there may be more regions
         * i.e.:
         *   |--------------------- e --------------------|
         *       |--- g ---|    |---- others ----|
         */
        return ret;   // don't need to do anything
    } else {
        /**
         * otherwise e overlaps some part of g
         *
         * figure out which region occurs first on the genome.  I.e., is it:
         * |------------- g ----------|
         *       |------------- e ----------|
         *
         * or:
         *       |------------- g ----------|
         * |------------ e -----------|
         *
         */

        if (that.getStart() < getStart()) {
            ret.push_back( GenomeLoc(getContig(), getContigIndex(), that.getStop() + 1, getStop()));
        } else {
            ret.push_back( GenomeLoc(getContig(), getContigIndex(), getStart(), that.getStart() - 1));
        }
        
        // replace g with the new region
        return ret;
    }
}

//@Requires("that != null")
bool GenomeLoc::containsP(const GenomeLoc & that) const {
    return onSameContig(that) && getStart() <= that.getStart() && getStop() >= that.getStop();
}

//@Requires("that != null")
bool GenomeLoc::onSameContig(const GenomeLoc & that) const {
    return (contigIndex == that.contigIndex);
}

//@Requires("that != null")
//@Ensures("result >= 0")
int GenomeLoc::distance( const GenomeLoc & that ) const {
    if ( onSameContig(that) )
        return abs(getStart() - that.getStart());
    else
        return INT_MAX;
}

//@Requires({"left != null", "right != null"})
bool GenomeLoc::isBetween( const GenomeLoc & left, const GenomeLoc & right ) {
    return compareTo(left) > -1 && compareTo(right) < 1;
}

/**
 * Tests whether this contig is completely before contig 'that'.
 * @param that Contig to test against.
 * @return true if this contig ends before 'that' starts; false if this is completely after or overlaps 'that'.
 */
//@Requires("that != null")
bool GenomeLoc::isBefore( GenomeLoc that ) const {
    int comparison = compareContigs(that);
    return ( comparison == -1 || ( comparison == 0 && getStop() < that.getStart() ));        
}

/**
 * Tests whether any portion of this contig is before that contig.
 * @param that Other contig to test.
 * @return True if the start of this contig is before the start of the that contig.
 */
//@Requires("that != null")
bool GenomeLoc::startsBefore(const GenomeLoc & that) {
    int comparison = compareContigs(that);
    return ( comparison == -1 || ( comparison == 0 && getStart() < that.getStart() ));
}

/**
 * Tests whether this contig is completely after contig 'that'.
 * @param that Contig to test against.
 * @return true if this contig starts after 'that' ends; false if this is completely before or overlaps 'that'.
 */
//@Requires("that != null")
bool GenomeLoc::isPast( const GenomeLoc & that ) {
    int comparison = compareContigs(that);
    return ( comparison == 1 || ( comparison == 0 && getStart() > that.getStop() ));
}

/**
 * Return the minimum distance between any pair of bases in this and that GenomeLocs:
 */
//@Requires("that != null")
//@Ensures("result >= 0")
int GenomeLoc::minDistance( const GenomeLoc & that ) {
    if (!onSameContig(that))
        return INT_MAX;
    
    int minDistance;
    if (isBefore(that))
        minDistance = distanceFirstStopToSecondStart(*this, that);
    else if (that.isBefore(*this))
        minDistance = distanceFirstStopToSecondStart(that, *this);
    else // this and that overlap [and possibly one contains the other]:
        minDistance = 0;
    
    return minDistance;
}

//@Requires({
//    "locFirst != null",
//    "locSecond != null",
//    "locSecond.isPast(locFirst)"
//})
//@Ensures("result >= 0")
int GenomeLoc::distanceFirstStopToSecondStart(GenomeLoc locFirst, GenomeLoc locSecond) {
    return locSecond.getStart() - locFirst.getStop();
}

//@Override
int GenomeLoc::hashCode() {
    return start << 16 | stop << 4 | contigIndex;
}


/**
 * conpare this genomeLoc's contig to another genome loc
 * @param that the genome loc to compare contigs with
 * @return 0 if equal, -1 if that.contig is greater, 1 if contig is greater
 */
//@Requires("that != null")
//@Ensures("result == 0 || result == 1 || result == -1")
int GenomeLoc::compareContigs( GenomeLoc that ) const {
    if (contigIndex == that.contigIndex)
        return 0;
    else if (contigIndex > that.contigIndex)
        return 1;
    return -1;
}

//@Requires("that != null")
//@Ensures("result == 0 || result == 1 || result == -1")
int GenomeLoc::compareTo( const GenomeLoc & that ) const {
    int result = 0;
    
    if ( *this == that ) {
        result = 0;
    }
    else if(isUnmapped(*this))
        result = 1;
    else if(isUnmapped(that))
        result = -1;
    else {
        int cmpContig = compareContigs(that);
        
        if ( cmpContig != 0 ) {
            result = cmpContig;
        } else {
            if ( getStart() < that.getStart() ) result = -1;
            if ( getStart() > that.getStart() ) result = 1;
        }
    }
    
    return result;
}

/**
 * reciprocialOverlap: what is the min. percent of gl1 and gl2 covered by both
 *
 * gl1.s ---------- gk1.e
 * gl2.s ---------- gl2.e
 * 100%
 *
 * gl1.s ---------- gk1.e
 *      gl2.s ---------- gl2.e
 * 50%
 *
 * gl1.s ---------- gk1.e
 *      gl2.s -------------------- gl2.e
 * 25% (50% for gl1 but only 25% for gl2)
 */
double GenomeLoc::reciprocialOverlapFraction(const GenomeLoc & o) {
    if ( overlapsP(o) )
        return min(overlapPercent(*this, o), overlapPercent(o, *this));
    else
        return 0.0;
}

long GenomeLoc::sizeOfOverlap( const GenomeLoc & that ) const {
    return ( overlapsP(that) ? min( getStop(), that.getStop() ) - max( getStart(), that.getStart() ) + 1L : 0L );
}


GenomeLoc GenomeLoc::UNMAPPED("");
GenomeLoc GenomeLoc::WHOLE_GENOME("all");
