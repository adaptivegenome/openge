/*********************************************************************
 *
 * GenomeLocParser.cpp: Port of GATK's GenomeLocParser.
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
 * The original file, GenomeLocParser.java, was released 
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

#include <iostream>
#include "GenomeLocParser.h"
#include <cassert>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>

#include <inttypes.h>
#include <errno.h>
using namespace BamTools;
using namespace std;

/**
 * Factory class for creating GenomeLocs
 */
//@Invariant({
//    "logger != null",
//    "contigInfo != null"})

#if 0

/**
 * set our internal reference contig order
 * @param refFile the reference file
 */
@Requires("refFile != null")
public GenomeLocParser(ReferenceSequenceFile refFile) {
    this(refFile.getSequenceDictionary());
}
#endif
GenomeLocParser::GenomeLocParser( BamSequenceRecords seqDict) {
    contigInfo = new BamSequenceRecords(seqDict);
    /*
    cerr << "Prepared reference sequence contig dictionary" << endl;
    for (SamSequenceConstIterator contig = contigInfo->Begin(); contig != contigInfo->End(); contig++) {
        fprintf(stderr, "\t%25s (%8d bp)\n", contig->Name.c_str(), atoi(contig->Length.c_str()));
    }
     */
}

/**
 * Determines whether the given contig is valid with respect to the sequence dictionary
 * already installed in the GenomeLoc.
 *
 * @return True if the contig is valid.  False otherwise.
 */
bool GenomeLocParser::contigIsInDictionary(const string contig) const {
    return contigInfo->contains(contig);
}

bool GenomeLocParser::indexIsInDictionary(const int index) const {
    return index >= 0 && index < contigInfo->size();
}
string GenomeLocParser::getContig(const int index) const
{
    assert(index >= 0 && index < contigInfo->size());
    return (*contigInfo)[index].name;
}

/**
 * get the contig's SAMSequenceRecord
 *
 * @param contig the string name of the contig
 *
 * @return the sam sequence record
 */
//@Ensures("result != null")
//@ThrowEnsures({"UserException.MalformedGenomeLoc", "!contigIsInDictionary(contig) || contig == null"})
BamSequenceRecord GenomeLocParser::getContigInfo(string contig) {
    if ( ! contigIsInDictionary(contig) )
        fprintf(stderr, "Contig %s given as location, but this contig isn't present in the Fasta sequence dictionary", contig.c_str());
    return (*contigInfo)[contig];
}

/**
 * Returns the contig index of a specified string version of the contig
 *
 * @param contig the contig string
 *
 * @return the contig index, -1 if not found
 */
//@Ensures("result >= 0")
//@ThrowEnsures({"UserException.MalformedGenomeLoc", "!contigIsInDictionary(contig) || contig == null"})
int GenomeLocParser::getContigIndex(string contig) const {
    return contigInfo->indexOfString(contig);
}

int GenomeLocParser::getContigIndexWithoutException(string contig) const {
    if ( ! contigInfo->contains(contig) )
        return -1;
    return contigInfo->indexOfString(contig);
}

// --------------------------------------------------------------------------------------------------------------
//
// Low-level creation functions
//
// --------------------------------------------------------------------------------------------------------------

/**
 * create a genome loc, given the contig name, start, and stop
 *
 * @param contig the contig name
 * @param start  the starting position
 * @param stop   the stop position
 *
 * @return a new genome loc
 */
//@Ensures("result != null")
//@ThrowEnsures({"UserException.MalformedGenomeLoc", "!isValidGenomeLoc(contig, start, stop)"})
GenomeLoc GenomeLocParser::createGenomeLoc(const string contig, const int start, const int stop) const {
    return createGenomeLoc(contig, getContigIndex(contig), start, stop);
}

GenomeLoc GenomeLocParser::createGenomeLoc(const string & contig, const int start, const int stop, bool mustBeOnReference) const {
    return createGenomeLoc(contig, getContigIndex(contig), start, stop, mustBeOnReference);
}

//@ThrowEnsures({"UserException.MalformedGenomeLoc", "!isValidGenomeLoc(contig, start, stop, false)"})
GenomeLoc GenomeLocParser::createGenomeLoc(const string & contig, int index, const int start, const int stop) const {
    return createGenomeLoc(contig, index, start, stop, false);
}

//@ThrowEnsures({"UserException.MalformedGenomeLoc", "!isValidGenomeLoc(contig, start, stop,mustBeOnReference)"})
GenomeLoc GenomeLocParser::createGenomeLoc(const string & contig, int index, const int start, const int stop, bool mustBeOnReference) const {
    //cerr << "CGL: " << contig << " ix: " << index << " " << start << ":" << stop << endl;
    validateGenomeLoc(contig, index, start, stop, mustBeOnReference, true);
    return GenomeLoc(contig, index, start, stop);
}
/**
 * validate a position or interval on the genome as valid
 *
 * Requires that contig exist in the master sequence dictionary, and that contig index be valid as well.  Requires
 * that start <= stop.
 *
 * if mustBeOnReference is true,
 * performs boundary validation for genome loc INTERVALS:
 * start and stop are on contig and start <= stop
 *
 * @param contig the contig name
 * @param start  the start position
 * @param stop   the stop position
 *
 * @return true if it's valid, false otherwise.  If exceptOnError, then throws a UserException if invalid
 */
bool GenomeLocParser::validateGenomeLoc(const string contig, const int contigIndex, const int start, const int stop, const bool mustBeOnReference, const bool exceptOnError) const {
    if ( ! contigInfo->contains(contig) ) {
        cerr << "validateGenomeLoc: Unknown contig " << contig << endl;
        assert(0);
        return false;
    }
    
    if (stop < start) {
        cerr << "validateGenomeLoc: The stop position " << stop << " is less than start " << start << " in contig " << contig << endl;
        assert(0);
        return false;
    }
    
    if (contigIndex < 0){
        cerr << "validateGenomeLoc: The contig index " << contigIndex << " is less than 0" << endl;
        assert(0);
        return false;
    }
    
    if (contigIndex >= contigInfo->size()) {
        cerr << "validateGenomeLoc: The contig index " << contigIndex << " is greater than the stored sequence count (" << contigInfo->size() << ")" << endl;
        assert(0);
        return false;
    }
    
    if ( mustBeOnReference ) {
        if (start < 0) {
            cerr << "validateGenomeLoc: The start position " << start << " is less than 0" << endl;
            assert(0);
            return false;
        }

        if (stop < 0) {
            cerr << "validateGenomeLoc: The stop position " << stop << " is less than 0" << endl;
            assert(0);
            return false;
        }

        int contigSize = (*contigInfo)[contigIndex].getLength();
        if (start > contigSize || stop > contigSize) {
            cerr << "validateGenomeLoc: The genome loc coordinates " << start << "-" << stop << " exceed the contig size (" << contigSize << ")" << endl;
            assert(0);
            return false;
        }
    }
    
    // we passed
    return true;
}

#if 0

public bool isValidGenomeLoc(string contig, int start, int stop, bool mustBeOnReference ) {
    return validateGenomeLoc(contig, getContigIndexWithoutException(contig), start, stop, mustBeOnReference, false);
}

public bool isValidGenomeLoc(string contig, int start, int stop ) {
    return validateGenomeLoc(contig, getContigIndexWithoutException(contig), start, stop, true, false);
}

private bool vglHelper(bool exceptOnError, string msg) {
    if ( exceptOnError )
        throw new UserException.MalformedGenomeLoc(string("Parameters to GenomeLocParser are incorrect:") + msg);
    else
        return false;
}

// --------------------------------------------------------------------------------------------------------------
//
// Parsing genome locs
//
// --------------------------------------------------------------------------------------------------------------
#endif
/**
 * parse a genome interval, from a location string
 *
 * Performs interval-style validation:
 *
 * contig is valid; start and stop less than the end; start <= stop, and start/stop are on the contig
 * @param str the string to parse
 *
 * @return a GenomeLoc representing the string
 *
 */
//@Requires("str != null")
//@Ensures("result != null")
GenomeLoc GenomeLocParser::parseGenomeLoc(const string str) const {
    // 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000'
    //System.out.printf("Parsing location '%s'%n", str);
    
    char * line = new char[str.size() + 1];
    strncpy(line, str.c_str(), str.size() + 1);
    
    char * colon_loc = strchr(line, ':');
    char * hyphen_loc = strchr(line, '-');
    
    if(!colon_loc) {  //parsing error, abort
        cerr << "Could not find colon in interval file: " << str << endl;
        exit(-1);
    }
    
    *colon_loc = 0;
    
    if(hyphen_loc) {  //must be a single
        *hyphen_loc = 0;
    }
    
    int start = strtol(colon_loc + 1, NULL, 10) - 1;
    
    if(start == 0 && errno == EINVAL) {
        cerr << "Could not find start in interval file: " << str << endl;
        perror("Reason:");
        exit(-1);
    }
    
    int stop = start;   //for the chr1:1 case
    if(hyphen_loc) {    // chr1:1-4 case
        stop = strtol(hyphen_loc + 1, NULL, 10) - 1;
        
        if(stop == 0 && errno == EINVAL) {
            cerr << "Could not find stop in interval file: " << str << endl;
            perror("Reason:");
            exit(-1);
        }
    }
    
    string contig(line);

    // is the contig valid?
    if (!contigIsInDictionary(contig)) {
        cerr << "Contig '" << contig << "' does not match any contig in the sequence dictionary derived from the reference; are you sure you are using the correct reference fasta file?" << endl;
        assert(0);
    }
    
    if (stop == INT_MAX)
        // lookup the actually stop position!
        stop = (*contigInfo)[contig].getLength();
    delete [] line;
    
    return createGenomeLoc(contig, getContigIndex(contig), start, stop, true);
}
#if 0

/**
 * Parses a number like 1,000,000 into a long.
 * @param pos
 */
@Requires("pos != null")
@Ensures("result >= 0")
private int parsePosition(string pos) {
    if(pos.indexOf('-') != -1) {
        throw new NumberFormatException("Position: '" + pos + "' can't contain '-'." );
    }
    
    if(strchr( pos.indexOf(',') != -1) {
        stringBuilder buffer = new stringBuilder();
        for(int i = 0; i < pos.length(); i++) {
            char c = pos.charAt(i);
            
            if(c == ',') {
                continue;
            } else if(c < '0' || c > '9') {
                throw new NumberFormatException("Position: '" + pos + "' contains invalid chars." );
            } else {
                buffer.append(c);
            }
        }
        return Integer.parseInt(buffer.tostring());
    } else {
        return Integer.parseInt(pos);
    }
}

// --------------------------------------------------------------------------------------------------------------
//
// Parsing string representations
//
// --------------------------------------------------------------------------------------------------------------
#endif
/**
 * create a genome loc, given a read. If the read is unmapped, *and* yet the read has a contig and start position,
 * then a GenomeLoc is returned for contig:start-start, otherwise and UNMAPPED GenomeLoc is returned.
 *
 * @param read
 *
 * @return
 */
//@Requires("read != null")
//@Ensures("result != null")
GenomeLoc GenomeLocParser::createGenomeLoc(const OGERead & read) const {
    if ( !read.IsMapped() && read.getRefID() == -1 )
        // read is unmapped and not placed anywhere on the genome
        return GenomeLoc::UNMAPPED;
    else {
        // Use max to ensure that end >= start (Picard assigns the end to reads that are entirely within an insertion as start-1)
        int length = read.getQueryBases().size();
        for( vector <CigarOp>::const_iterator i = read.getCigarData().begin(); i != read.getCigarData().end(); i++) {
            if(i->type == 'D')
                length += i->length;
            if(i->type == 'I')
                length -= i->length;
            if(i->type == 'S' || i->type == 'H')
                length -= i->length;
        }

        int end = read.IsMapped() ? max(read.getPosition() + length, read.getPosition()):read.getPosition() ;
        //cerr << "CGL READ: " << getContig(read.RefID) << " Pos: " << read.Position << ":" << end << " len:" << length << " mapped: " << (int) read.IsMapped() << "read " << read.QueryBases[0] << endl;

        return createGenomeLoc(getContig(read.getRefID()), read.getPosition(), max(end-1, read.getPosition()));
    }
}
#if 0
/**
 * Creates a GenomeLoc from a Tribble feature
 * @param feature
 * @return
 */
public GenomeLoc createGenomeLoc(Feature feature) {
    return createGenomeLoc(feature.getChr(), feature.getStart(), feature.getEnd());
}

/**
 * Creates a GenomeLoc corresponding to the variant context vc.  If includeSymbolicEndIfPossible
 * is true, and VC is a symbolic allele the end of the created genome loc will be the value
 * of the END info field key, if it exists, or vc.getEnd() if not.
 *
 * @param vc
 * @param includeSymbolicEndIfPossible
 * @return
 */
public GenomeLoc createGenomeLoc(VariantContext vc, bool includeSymbolicEndIfPossible) {
    if ( includeSymbolicEndIfPossible && vc.isSymbolic() ) {
        int end = vc.getAttributeAsInt(VCFConstants.END_KEY, vc.getEnd());
        return createGenomeLoc(vc.getChr(), vc.getStart(), end);
    }
    else
        return createGenomeLoc(vc.getChr(), vc.getStart(), vc.getEnd());
}

public GenomeLoc createGenomeLoc(VariantContext vc) {
    return createGenomeLoc(vc, false);
}

/**
 * create a new genome loc, given the contig name, and a single position. Must be on the reference
 *
 * @param contig the contig name
 * @param pos    the postion
 *
 * @return a genome loc representing a single base at the specified postion on the contig
 */
@Ensures("result != null")
@ThrowEnsures({"UserException.MalformedGenomeLoc", "!isValidGenomeLoc(contig, pos, pos, true)"})
public GenomeLoc createGenomeLoc(string contig, int pos) {
    return createGenomeLoc(contig, getContigIndex(contig), pos, pos);
}

/**
 * create a new genome loc from an existing loc, with a new start position
 * Note that this function will NOT explicitly check the ending offset, in case someone wants to
 * set the start of a new GenomeLoc pertaining to a read that goes off the end of the contig.
 *
 * @param loc   the old location
 * @param start a new start position
 *
 * @return the newly created genome loc
 */
public GenomeLoc setStart(GenomeLoc loc, int start) {
    return createGenomeLoc(loc.getContig(), loc.getContigIndex(), start, loc.getStop());
}

/**
 * create a new genome loc from an existing loc, with a new stop position
 * Note that this function will NOT explicitly check the ending offset, in case someone wants to
 * set the stop of a new GenomeLoc pertaining to a read that goes off the end of the contig.
 *
 * @param loc  the old location
 * @param stop a new stop position
 *
 * @return
 */
public GenomeLoc setStop(GenomeLoc loc, int stop) {
    return createGenomeLoc(loc.getContig(), loc.getContigIndex(), loc.start, stop);
}

/**
 * return a new genome loc, with an incremented position
 *
 * @param loc the old location
 *
 * @return a new genome loc
 */
public GenomeLoc incPos(GenomeLoc loc) {
    return incPos(loc, 1);
}

/**
 * return a new genome loc, with an incremented position
 *
 * @param loc the old location
 * @param by  how much to move the start and stop by
 *
 * @return a new genome loc
 */
public GenomeLoc incPos(GenomeLoc loc, int by) {
    return createGenomeLoc(loc.getContig(), loc.getContigIndex(), loc.start + by, loc.stop + by);
}

/**
 * Creates a GenomeLoc than spans the entire contig.
 * @param contigName Name of the contig.
 * @return A locus spanning the entire contig.
 */
@Requires("contigName != null")
@Ensures("result != null")
public GenomeLoc createOverEntireContig(string contigName) {
    SAMSequenceRecord contig = contigInfo.getSequence(contigName);
    return createGenomeLoc(contigName,contig.getSequenceIndex(),1,contig.getSequenceLength(), true);
}

/**
 * Creates a loc to the left (starting at the loc start + 1) of maxBasePairs size.
 * @param loc The original loc
 * @param maxBasePairs The maximum number of basePairs
 * @return The contiguous loc of up to maxBasePairs length or null if the loc is already at the start of the contig.
 */
@Requires({"loc != null", "maxBasePairs > 0"})
public GenomeLoc createGenomeLocAtStart(GenomeLoc loc, int maxBasePairs) {
    if (GenomeLoc.isUnmapped(loc))
        return null;
    string contigName = loc.getContig();
    SAMSequenceRecord contig = contigInfo.getSequence(contigName);
    int contigIndex = contig.getSequenceIndex();
    
    int start = loc.getStart() - maxBasePairs;
    int stop = loc.getStart() - 1;
    
    if (start < 1)
        start = 1;
    if (stop < 1)
        return null;
    
    return createGenomeLoc(contigName, contigIndex, start, stop, true);
}

/**
 * Creates a loc to the right (starting at the loc stop + 1) of maxBasePairs size.
 * @param loc The original loc
 * @param maxBasePairs The maximum number of basePairs
 * @return The contiguous loc of up to maxBasePairs length or null if the loc is already at the end of the contig.
 */
@Requires({"loc != null", "maxBasePairs > 0"})
public GenomeLoc createGenomeLocAtStop(GenomeLoc loc, int maxBasePairs) {
    if (GenomeLoc.isUnmapped(loc))
        return null;
    string contigName = loc.getContig();
    SAMSequenceRecord contig = contigInfo.getSequence(contigName);
    int contigIndex = contig.getSequenceIndex();
    int contigLength = contig.getSequenceLength();
    
    int start = loc.getStop() + 1;
    int stop = loc.getStop() + maxBasePairs;
    
    if (start > contigLength)
        return null;
    if (stop > contigLength)
        stop = contigLength;
    
    return createGenomeLoc(contigName, contigIndex, start, stop, true);
}
#endif
