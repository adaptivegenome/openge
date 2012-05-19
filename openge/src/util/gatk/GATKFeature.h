//
//  GATKFeature.h
//  OpenGE
//
//  Created by Lee Baker on 5/13/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#ifndef OpenGE_GATKFeature_h
#define OpenGE_GATKFeature_h

/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
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

class Feature 
{
public:
    std::string chr;
    int start;
    int end;
    
    bool operator<(const Feature & o) const {
        if(start != o.start)
            return start < o.start;
        if(end != o.end)
            return end < o.end;
        return chr < o.chr;
    }
    
    bool operator==(const Feature & o) const {
        return start == o.start && end == o.end && chr == o.chr;
    }
    bool operator!=(const Feature & o) const { return !(*this == o); }
};

//this is an abstract class- what actually implements it though?
class GATKFeature { // implements Feature, HasGenomeLocation {
protected:
    std::string name;
    const GenomeLocParser * genomeLocParser;
    Feature feature;
    GenomeLoc * position;
    
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
    
    GATKFeature(const GenomeLocParser * genomeLocParser,Feature f, std::string name) 
    : name(name)
    , genomeLocParser(genomeLocParser)
    , feature(f)
    , position(NULL)
    { }
    
    std::string getName() {
        return name;
    }
    
    GenomeLoc getLocation() {
        if (position == NULL) 
            position = new GenomeLoc(genomeLocParser->createGenomeLoc(feature.chr, feature.start, feature.end));
        return *position;
    }
    
    /** Return the features reference sequence name, e.g chromosome or contig */
    std::string getChr() {
        return feature.chr;
    }
    
    /** Return the start position in 1-based coordinates (first base is 1) */
    int getStart() {
        return feature.start;
    }
    
    /**
     * Return the end position following 1-based fully closed conventions.  The length of a feature is
     * end - start + 1;
     */
    int getEnd() {
        return feature.end;
    }
    
    // TODO: this should be a Feature, actually
    Feature & getUnderlyingObject() {
        return feature;
    }
    
    bool operator<(const GATKFeature & o) const {
        if(feature != o.feature)
            return feature < feature;
        if(*position != *(o.position))
            return *position < *position;
        return name < o.name;
    }
};

#endif
