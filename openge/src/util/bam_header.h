#ifndef OGE_BAM_HEADER_H
#define OGE_BAM_HEADER_H
/*********************************************************************
 *
 * bam_header.h: Store information related to a BAM/SAM header
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 17 Dec 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial
 * Purpose License. A copy of this license has been provided in
 * the openge/ directory.
 *
 *********************************************************************/

#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <assert.h>

class BamSequenceRecord {
public:
	std::string assemblyId;	//AS:
	std::string checksum;	//M5:
	size_t length;	//LN:
	std::string name;	//SN:
	std::string species;	//SP:
	std::string uri;	//UR:

    BamSequenceRecord() : length(-1) {}
    BamSequenceRecord(const std::string & s);
	BamSequenceRecord(const BamSequenceRecord & r)
	: assemblyId(r.assemblyId)
	, checksum(r.checksum)
	, length(r.length)
	, name(r.name)
	, species(r.species)
	, uri(r.uri)
	{}
    
    bool operator==(const BamSequenceRecord & r) const {
        return assemblyId == r.assemblyId &&
        checksum == r.checksum &&
        length == r.length &&
        name == r.name &&
        species == r.species &&
        uri == r.uri;
    }

	const std::string & getAssemblyId() const { return assemblyId; }
	const std::string & getChecksum() const { return checksum; }
	const size_t getLength() const { return length; }
	const std::string & getName() const { return name; }
	const std::string & getSpecies() const { return species; }
	const std::string & getURI() const { return uri; }

	const std::string toString() const;	//TODO
};

class BamProgramRecord {
public:
	std::string commandLine;
	std::string id;
	std::string name;
	std::string version;
	std::string previousProgramId;
	std::string nextProgramId;

    BamProgramRecord() {}

	BamProgramRecord(const BamProgramRecord & r)
	: commandLine(r.commandLine)
	, id(r.id)
	, name(r.name)
	, version(r.version)
	, previousProgramId(r.previousProgramId)
	, nextProgramId(r.nextProgramId)
	{ }

	BamProgramRecord(const std::string & s);	//TODO implement

	const std::string & getCommandLine() const { return commandLine; }
	const std::string & getID() const { return id; }
	const std::string & getName() const { return name; }
	const std::string & getVersion() const { return version; }
	const std::string & getPreviousProgramId() const { return previousProgramId; }
	const std::string & getNextProgramId() const { return nextProgramId; }

	const std::string toString() const;	//TODO
};

class BamReadGroupRecord{
	std::string description;
	std::string flowOrder;
	std::string id;
	std::string keySequence;
	std::string library;
	std::string platformUnit;
	std::string predictedInsertionSize;
	std::string productionDate;
	std::string program;
	std::string sample;
	std::string sequencingCenter;
	std::string sequencingTechnology;
public:
	BamReadGroupRecord(const std::string & s);	//TODO
	BamReadGroupRecord(const BamReadGroupRecord & r)
	: description(r.description)
	, flowOrder(r.flowOrder)
	, id(r.id)
	, keySequence(r.keySequence)
	, library(r.library)
	, platformUnit(r.platformUnit)
	, predictedInsertionSize(r.predictedInsertionSize)
	, productionDate(r.productionDate)
	, program(r.program)
	, sample(r.sample)
	, sequencingCenter(r.sequencingCenter)
	, sequencingTechnology(r.sequencingTechnology)
	{ }

	const std::string & getDescription () const { return description; }
	const std::string & getFlowOrder() const { return flowOrder; }
	const std::string & getId() const { return id; }
	const std::string & getKeySequence() const { return keySequence; }
	const std::string & getLibrary() const { return library; }
	const std::string & getPlatformUnit() const { return platformUnit; }
	const std::string & getProductionDate() const { return productionDate; }
	const std::string & getProgram() const { return program; }
	const std::string & getSample() const { return sample; }
	const std::string & getSequencingCenter() const { return sequencingCenter; }
	const std::string & getSequencingTechnology() const { return sequencingTechnology; }

	const std::string toString() const;	//TODO
};

class BamSequenceRecords : public std::vector<BamSequenceRecord> {
	//sequence records order should be maintained, need to figure this out
    
public:
    void add(const BamSequenceRecord & r) {
        push_back(r);
    }
    
    bool contains(const std::string & n) const {
        for(const_iterator i = begin(); i != end(); i++)
            if(i->getName() == n)
                return true;
        return false;
    }
    
    BamSequenceRecord & operator[] ( const int n) {
        assert(n >= 0 && n < size());
        return *(begin() + n);
    }
    
    const BamSequenceRecord & operator[] ( const int n) const {
        assert(n >= 0 && n < size());
        return *(begin() + n);
    }

    BamSequenceRecord operator[] ( const std::string & n) {
        assert(contains(n));
        for(const_iterator i = begin(); i != end(); i++)
            if(i->getName() == n)
                return *i;
        std::cerr << "Invalid sequence " << n << ". Quitting." << std::endl;
        exit(-1);
        return *begin(); //Execution can never get here
    }
    
    const BamSequenceRecord & operator[] ( const std::string & n) const {
        assert(contains(n));
        for(const_iterator i = begin(); i != end(); i++)
            if(i->getName() == n)
                return *i;
        std::cerr << "Invalid sequence " << n << ". Quitting." << std::endl;
        exit(-1);
        return *begin(); //Execution can never get here
    }

    int indexOfString(const std::string s) const {
        assert(contains(s));
        
        for(int i = 0; i < size(); i++)
            if(at(i).getName() == s)
                return i;
        return -1;
    }

	const std::string toString() const;
};

class BamProgramRecords : public std::vector<BamProgramRecord> {
public:
    const bool contains(const std::string name) {
        for(const_iterator i = begin(); i != end(); i++)
            if(name == i->getName())
                return true;
        return false;
    }
    
    void add(const BamProgramRecord & r) {
        push_back(r);
    }

	const std::string toString();
};

class BamReadGroupRecords : public std::vector<BamReadGroupRecord> {
public:
    const bool contains(const std::string id) const {
        for(const_iterator i = begin(); i != end(); i++)
            if(id == i->getId())
                return true;
        return false;
    }

    BamReadGroupRecord operator[] ( const std::string & id) {
        assert(contains(id));
        for(const_iterator i = begin(); i != end(); i++)
            if(i->getId() == id)
                return *i;
        std::cerr << "Invalid sequence " << id << ". Quitting." << std::endl;
        exit(-1);
        return *begin(); //Execution can never get here
    }
    
    const BamReadGroupRecord & operator[] ( const std::string & id) const {
        assert(contains(id));
        for(const_iterator i = begin(); i != end(); i++)
            if(i->getId() == id)
                return *i;
        std::cerr << "Invalid sequence " << id << ". Quitting." << std::endl;
        exit(-1);
        return *begin(); //Execution can never get here
    }

	const std::string toString();
};

class BamHeader {
	BamSequenceRecords sq;
	BamProgramRecords pg;
	BamReadGroupRecords rg;
	std::vector<std::string> co;

public:
	typedef enum { SORT_UNKNOWN, SORT_UNSORTED, SORT_QUERYNAME, SORT_COORDINATE 
	} sort_order_t;

protected:
	std::string format_version;
	sort_order_t sort_order;

public:
	BamHeader()
	: sort_order(SORT_UNKNOWN)
	{}

	BamHeader(const BamHeader & h)
	: sq(h.sq)
	, pg(h.pg)
	, rg(h.rg)
	, co(h.co)
	, format_version(h.format_version)
	, sort_order(h.sort_order)
	{}

	BamHeader (const std::string & text);
    
	const BamSequenceRecords & getSequences() const { return sq; }
    BamSequenceRecords & getSequences() { return sq; }
	const BamProgramRecords & getPrograms() const { return pg; }
    BamProgramRecords & getPrograms() { return pg; }
	const BamReadGroupRecords & getReadGroups() const { return rg; }
    BamReadGroupRecords & getReadGroups() { return rg; }
	const std::vector<std::string> & getComments() const { return co; }
    std::vector<std::string> & getComments() { return co; }

	const sort_order_t getSortOrder()  const { return sort_order; }
	void setSortOrder(const sort_order_t s) { sort_order = s; }

	const std::string toString() const;	//TODO
};

#endif