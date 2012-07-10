#ifndef OpenGE_GATKFeature_h
#define OpenGE_GATKFeature_h

/*********************************************************************
 *
 * AlignmentUtils.cpp: Port of GATK's AlignmentUtils.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 13 May 2012
 *
 *********************************************************************
 *
 * This file has been ported from GATK's implementation in Java, and
 * is released under the Virginia Tech Non-Commercial Purpose License.
 * A copy of this license has been provided in  the openge/ directory.
 * 
 * The original file, AlignmentUtils.java, was released 
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
/*
package org.broadinstitute.sting.gatk.refdata.utils;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.HasGenomeLocation;
*/

/**
 * 
 * @author aaron 
 * 
 * Class GATKFeature
 *
 * This wraps a Tribble feature or a RODatum so that both present the same interface: a genome loc for position and a
 * way of retrieving the track name.
 */

#include "GenomeLoc.h"
#include "GenomeLocParser.h"

class GATKFeature { // implements Feature, HasGenomeLocation {
protected:
    std::string name;
    const GenomeLocParser * genomeLocParser;
    GenomeLoc * position;

    std::string chr;
    int start;
    int end;
    
protected:
    void setName(std::string name) {
        this->name = name;
    }
    
public:
    GATKFeature()
    : name("")
    , genomeLocParser(NULL)
    , position(NULL)
    {}
    
    GATKFeature(const GenomeLocParser * genomeLocParser, std::string name) 
    : name(name)
    , genomeLocParser(genomeLocParser)
    , position(NULL)
    { }
    
    std::string getName() {
        return name;
    }
    
    GenomeLoc getLocation() {
        if (position == NULL) 
            position = new GenomeLoc(genomeLocParser->createGenomeLoc(chr, start, end));
        return *position;
    }
    
    /** Return the features reference sequence name, e.g chromosome or contig */
    std::string getChr() {
        return chr;
    }
    
    /** Return the start position in 1-based coordinates (first base is 1) */
    int getStart() {
        return start;
    }
    
    /**
     * Return the end position following 1-based fully closed conventions.  The length of a feature is
     * end - start + 1;
     */
    int getEnd() {
        return end;
    }
    
    bool operator<(const GATKFeature & o) const {
        if(chr != o.chr)
            return chr < o.chr;
        if(start != o.start)
            return start < o.start;
        if(end != o.end)
            return end < o.end;
        if(*position != *(o.position))
            return *position < *position;
        return name < o.name;
    }
};

#endif
