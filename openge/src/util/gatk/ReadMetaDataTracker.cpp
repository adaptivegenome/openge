//
//  ReadMetaDataTracker.cpp
//  OpenGE
//
//  Ported to C++ by Lee Baker on 5/8/12.
//

#include <iostream>
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
package org.broadinstitute.sting.gatk.refdata;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.datasources.providers.RODMetaDataContainer;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.*;

*/
#include "ReadMetaDataTracker.h"
#include "RODMetaDataContainer.h"
#include "GATKFeature.h"


#include <string>
#include <map>

using namespace std;
using namespace BamTools;


ReadMetaDataTracker::ReadMetaDataTracker(GenomeLocParser * genomeLocParser, BamAlignment * record, map<int, RODMetaDataContainer> mapping) 
: genomeLocParser(genomeLocParser)
, record(record)
, mapping(mapping)
{
}

map<int, set<GATKFeature> > ReadMetaDataTracker::createReadAlignment(BamAlignment * record, map<int, RODMetaDataContainer> q, string name) {
    map<int, set<GATKFeature> > ret;
    GenomeLoc location = genomeLocParser->createGenomeLoc(*record);
    int length = record->Length;
    for(map<int, RODMetaDataContainer>::iterator loc_iter = q.begin(); loc_iter != q.end(); loc_iter++) {
        const int loc = loc_iter->first;
        int position = loc - location.getStart();
        if (position >= 0 && position < length) {
            set<GATKFeature> set = q[loc].getSet(name);
            if (set.size() > 0)
                ret[position] = set;
        }
    }
    return ret;
    
}

map<int, set<GATKFeature> > ReadMetaDataTracker::createGenomeLocAlignment(const BamAlignment & record, map<int, RODMetaDataContainer> mapping, string name) {
    map<int, set<GATKFeature> > ret;
    int start = record.Position;
    int stop = record.Position + record.Length;
    for (map<int, RODMetaDataContainer>::iterator it = mapping.begin(); it != mapping.end(); it++) {
        int location = it->first;
        if (location >= start && location <= stop)
            ret[location] = mapping[location].getSet(name);
    }
    return ret;
}

map<int, set<GATKFeature> > ReadMetaDataTracker::getReadOffsetMapping() {
    return createReadAlignment(record, mapping, NULL);
}

map<int, set<GATKFeature> > ReadMetaDataTracker::getContigOffsetMapping() {
    return createGenomeLocAlignment(*record, mapping, NULL);
}

map<int, set<GATKFeature> > ReadMetaDataTracker::getReadOffsetMapping(string name) {
    return createReadAlignment(record, mapping, name);
}

map<int, set<GATKFeature> > ReadMetaDataTracker::getContigOffsetMapping(string name) {
    return createGenomeLocAlignment(*record, mapping, name);
}
#if 0
map<int, vector<GATKFeature> > ReadMetaDataTracker::getReadOffsetMapping(Class cl) {
    return createReadAlignment(record, mapping, cl, NULL);
}

map<int, vector<GATKFeature> > getContigOffsetMapping(Class cl) {
    return createGenomeLocAlignment(record, mapping, cl, null);
}
#endif

vector<GATKFeature> ReadMetaDataTracker::getAllCoveringRods() {
    /*
     
     List<GATKFeature> ret = new ArrayList<GATKFeature>();
     for (Map.Entry<Integer, RODMetaDataContainer> entry : mapping.entrySet())
     ret.addAll(entry.getValue().getSet());
     return ret;
     */
    
    vector<GATKFeature> ret;
    for (map<int, RODMetaDataContainer>::iterator it = mapping.begin(); it != mapping.end(); it++) {
        set<GATKFeature> is = it->second.getSet();
        ret.insert(ret.begin(), is.begin(), is.end());
        //ret.addAll(entry.getValue().getSet());
    }
    //for (map.Entry<int, RODMetaDataContainer> entry : mapping.entrySet())
    return ret;
}
