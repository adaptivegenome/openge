//
//  ReferenceContext.h
//  OpenGE
//
//  Created by Lee Baker on 5/9/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

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
