/*********************************************************************
 *
 * ReadMetaDataTracker.h: Port of GATK's ReadMetaDataTracker.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 8 May 2012
 *
 *********************************************************************
 *
 * This file has been ported from GATK's implementation in Java, and
 * is released under the Virginia Tech Non-Commercial Purpose License.
 * A copy of this license has been provided in  the openge/ directory.
 * 
 * The original file, ReadMetaDataTracker.java, was released 
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
    std::map<int, std::set<GATKFeature> > createGenomeLocAlignment(const BamTools::BamAlignment & record, std::map<int, RODMetaDataContainer> mapping, std::string * name);
    
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
