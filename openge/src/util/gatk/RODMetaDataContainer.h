#ifndef OpenGE_RODMetaDataContainer_h
#define OpenGE_RODMetaDataContainer_h

/*********************************************************************
 *
 * RODMetaDataContainer.cpp: Port of GATK's RODMetaDataContainer.
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
 * The original file, RODMetaDataContainer.java, was released 
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
package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.*;
*/

/**
 * 
 * @author aaron 
 * 
 * Class RODMetaDataContainer
 *
 * stores both the name and the class for each ROD.  This class assumes that:
 *
 * -Names must be unique
 * -Classes are allowed to have duplicates
 *
 * This class encapsulates the ref data associations, and provides lookup by name and by
 * class type.
 *
 */

#include <string>
#include <map>
#include <set>
#include <utility>

#include "GATKFeature.h"

class RODMetaDataContainer {
protected:
    
    // we only allow non-duplicate ROD names, a HashMap is fine
    std::map<std::string, GATKFeature> nameMap;

public:
    void addEntry(GATKFeature & data) {
        nameMap.insert(std::pair<std::string, GATKFeature>(data.getName(), data));
        //nameMap[data.getName()] = data;
    }
    
    std::set<GATKFeature> getSet(std::string name) {
        if (name.size() == 0) return getSet();
        std::set<GATKFeature> set;
        if (nameMap.count(name)) 
            set.insert(nameMap[name]);
        return set;
    }
    
    /**
     * get the feature contents of this container; the unfiltered set without their name association
     * @return
     */
    std::set<GATKFeature> getSet() {
        std::set<GATKFeature> ret;
        for(std::map<std::string, GATKFeature>::iterator nm_it = nameMap.begin(); nm_it != nameMap.end(); nm_it++)
            ret.insert(nm_it->second);
        return ret;
    }
#if 0
    // the brute force (n) search ended up being faster than sorting and binary search in all but the most extreme cases (thousands of RODs at a location).
    std::set<GATKFeature> getSet(Class cls) {
        std::set<GATKFeature> ret;
        for (pair<Class, GATKFeature> pair: classMap)
            if (pair.first.equals(cls)) ret.add(pair.second);
        return ret;
    }
#endif
};

#endif
