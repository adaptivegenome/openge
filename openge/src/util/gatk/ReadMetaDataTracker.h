//
//  ReadMetaDataTracker.h
//  OpenGE
//
//  Created by Lee Baker on 5/8/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#ifndef OpenGE_ReadMetaDataTracker_h
#define OpenGE_ReadMetaDataTracker_h

#include "api/BamAlignment.h"
#include "GenomeLoc.h"
#include "GenomeLocParser.h"
#include "RODMetaDataContainer.h"

#include <string>
#include <map>


/**
 * @author aaron
 *         <p/>
 *         Class ReadMetaDataTracker
 *         <p/>
 *         a read-based meta data tracker
 */
class ReadMetaDataTracker {
private:
    /**
     * The parser, used to create new GenomeLocs.
     */
    GenomeLocParser * genomeLocParser;
    
    BamTools::BamAlignment * record;
    
    // the buffer of positions and RODs we've stored
    std::map<int, RODMetaDataContainer> mapping;
    
public:
    /**
     * create a read meta data tracker, given the read and a queue of RODatum positions
     *
     * @param record  the read to create offset from
     * @param mapping the mapping of reference ordered datum
     */
    ReadMetaDataTracker(GenomeLocParser * genomeLocParser, BamTools::BamAlignment * record, std::map<int, RODMetaDataContainer> mapping) ;
    
    /**
     * create an alignment of read position to reference ordered datum
     *
     * @param record the SAMRecord
     * @param queue  the queue (as a tree set)
     * @param cl     the class name, null if not filtered by classname
     * @param name   the datum track name, null if not filtered by name
     *
     * @return a mapping from the position in the read to the reference ordered datum
     */
private:
    std::map<int, std::set<GATKFeature> > createReadAlignment(BamTools::BamAlignment * record, std::map<int, RODMetaDataContainer> queue, std::string name);
    
    /**
     * create an alignment of read position to reference ordered datum
     *
     * @return a mapping from the position in the read to the reference ordered datum
     */
    std::map<int, std::set<GATKFeature> > createGenomeLocAlignment(const BamTools::BamAlignment & record, std::map<int, RODMetaDataContainer> mapping, std::string name);
    
    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of read offset to ROD(s)
     */
public:
    std::map<int, std::set<GATKFeature> > getReadOffsetMapping() ;    
    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of genome loc position to ROD(s)
     */
    std::map<int, std::set<GATKFeature> > getContigOffsetMapping() ;
    
    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of read offset to ROD(s)
     */
    std::map<int, std::set<GATKFeature> > getReadOffsetMapping(std::string name) ;    
    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of genome loc position to ROD(s)
     */
    std::map<int, std::set<GATKFeature> > getContigOffsetMapping(std::string name) ;       
    /**
     * get the list of all the RODS overlapping this read, without any information about their position
     * @return a Collection (no order guaranteed), of all the RODs covering this read
     */
    std::vector<GATKFeature> getAllCoveringRods() ;
};

#endif
