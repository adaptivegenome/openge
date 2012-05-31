/*********************************************************************
 *
 * ReferenceContext.h: Port of GATK's ReferenceContext.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 9 May 2012
 *
 *********************************************************************
 *
 * This file has been ported from GATK's implementation in Java, and
 * is released under the Virginia Tech Non-Commercial Purpose License.
 * A copy of this license has been provided in  the openge/ directory.
 * 
 * The original file, ReferenceContext.java, was released 
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

#ifndef OpenGE_ReferenceContext_h
#define OpenGE_ReferenceContext_h

#include <string>

#include "GenomeLoc.h"
#include "GenomeLocParser.h"

class ReferenceContext {
public:
    const static bool UPPERCASE_REFERENCE;
    
    class ForwardingProvider {
        std::string bases;
        
    public:
        ForwardingProvider( char base ) 
        : bases(std::string(1,base))
        { }
        
        ForwardingProvider( std::string bases )
        : bases(bases)
        { }
        
        std::string getBases() { return bases; }
    };
    
    /**
     * Facilitates creation of new GenomeLocs.
     */
private:
    GenomeLocParser genomeLocParser;
    
    /**
     * The locus.
     */
    GenomeLoc locus;
    
    /**
     * The window of reference information around the current locus.
     */
    GenomeLoc * window;
    
    /**
     * The bases in the window around the current locus.  If null, then bases haven't been fetched yet
     */
    std::string * basesCache;
    
    /**
     * Lazy loader to fetch reference bases
     */
    ForwardingProvider basesProvider;
    
    /**
     * Interface to create byte[] contexts for lazy loading of the reference
     */
public:
    
    ReferenceContext( GenomeLocParser genomeLocParser, GenomeLoc locus, char base ) ;
    ReferenceContext( GenomeLocParser genomeLocParser, GenomeLoc locus, GenomeLoc window, std::string bases );
    ReferenceContext( GenomeLocParser genomeLocParser, GenomeLoc locus, GenomeLoc window, ForwardingProvider basesProvider );
    /**
     * Utility function to load bases from the provider to the cache, if necessary
     */
private:
    void fetchBasesFromProvider();
    
public:
    GenomeLocParser getGenomeLocParser();
    GenomeLoc getLocus();
    GenomeLoc getWindow();    
    char getBase();
    std::string getBases();
};

#endif
