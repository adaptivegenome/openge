//
//  ReferenceContext.cpp
//  OpenGE
//
//  Created by Lee Baker on 5/9/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//
#include "ReferenceContext.h"

#include <iostream>
#include <string>
using namespace std;

/*
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
/*
package org.broadinstitute.sting.gatk.contexts;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.util.StringUtil;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
*/
/**
 * The section of the reference that overlaps with the given
 * read / locus. 
 *
 * @author hanna
 * @version 0.1
 */


const bool ReferenceContext::UPPERCASE_REFERENCE = true;

/**
 * Contructor for a simple, windowless reference context.
 * @param locus locus of interest.
 * @param base reference base at that locus.
 */
ReferenceContext::ReferenceContext( GenomeLocParser genomeLocParser, GenomeLoc locus, char base ) 
: genomeLocParser(genomeLocParser)
, locus(locus)
, window(NULL)
, basesCache(NULL)
, basesProvider(ForwardingProvider(base))
{}
ReferenceContext::ReferenceContext( GenomeLocParser genomeLocParser, GenomeLoc locus, GenomeLoc window, string bases ) 
: genomeLocParser(genomeLocParser)
, locus(locus)
, window(new GenomeLoc(window))
, basesCache(NULL)
, basesProvider(ForwardingProvider(bases))
{}

ReferenceContext::ReferenceContext( GenomeLocParser genomeLocParser, GenomeLoc locus, GenomeLoc window, ForwardingProvider basesProvider )
: genomeLocParser(genomeLocParser)
, locus(locus)
, window(new GenomeLoc(window))
, basesCache(NULL)
, basesProvider(basesProvider)
{}

void ReferenceContext::fetchBasesFromProvider() {
    if ( basesCache == NULL ) {
        basesCache = new string(basesProvider.getBases());
        if (UPPERCASE_REFERENCE) 
            std::transform(basesCache->begin(), basesCache->end(), basesCache->begin(), (int (*)(int))std::toupper); 
    }
}

/**
 * @return The genome loc parser associated with this reference context
 */
GenomeLocParser ReferenceContext::getGenomeLocParser() {
    return genomeLocParser;
}

/**
 * The locus currently being examined.
 * @return The current locus.
 */
GenomeLoc ReferenceContext::getLocus() {
    return locus;
}

GenomeLoc ReferenceContext::getWindow() {
    return *window;
}

/**
 * Get the base at the given locus.
 * @return The base at the given locus from the reference.
 */
char ReferenceContext::getBase() {
    return getBases()[(int)(locus.getStart() - window->getStart())];
}

/**
 * All the bases in the window currently being examined.
 * @return All bases available.  If the window is of size [0,0], the array will
 *         contain only the base at the given locus.
 */
string ReferenceContext::getBases() {
    fetchBasesFromProvider();
    return *basesCache;
}
