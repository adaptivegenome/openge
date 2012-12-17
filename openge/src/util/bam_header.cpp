/*********************************************************************
 *
 * bam_header.cpp: Store information related to a BAM/SAM header
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

#include "bam_header.h"

#include <algorithm>
#include <sstream>
using namespace std;

//////////////////////////////////////////
// Parsing
//take a header line (minus the leading @XX\t) and convert to a list of (tag,data) pairs
vector<string> headerLineSplit(const string & line) {
	vector<string> ret;
	size_t ct = 1 + count(line.begin(), line.end(), '\t');
	ret.reserve(ct);

	string::const_iterator i = line.begin();
	while(ct--) {
		string::const_iterator next = find(i, line.end(), '\t');
		ret.push_back(string(i, next));
        i = next+1;
	}

	return ret;
}

BamSequenceRecord::BamSequenceRecord(const string & s) {
	const vector<string> segments = headerLineSplit(s);

	for(vector<string>::const_iterator seg = segments.begin(); seg != segments.end(); seg++) {
		const string tag = seg->substr(0,2);
		const string data = seg->substr(3);

		if(tag == "SN") name = data;
		else if(tag == "LN") length = atoi(data.c_str());
		else if(tag == "AS") assemblyId = data;
		else if(tag == "M5") checksum = data;
		else if(tag == "SP") species = data;
		else if(tag == "UR") uri = data;
	}

	if(name.empty() || length == -1) {
		cerr << "Mandatory field missing in header sequence line." << endl;
		exit(-1);
	}
}

BamProgramRecord::BamProgramRecord(const string & s) {
	const vector<string> segments = headerLineSplit(s);

	for(vector<string>::const_iterator seg = segments.begin(); seg != segments.end(); seg++) {
		const string tag = seg->substr(0,2);
		const string data = seg->substr(3);

		if(tag == "ID") id = data;
		else if(tag == "PN") name = data;
		else if(tag == "CL") commandLine = data;
		else if(tag == "PP") previousProgramId = data;
		else if(tag == "VN") version = data;
	}

	if(id.empty()) {
		cerr << "Mandatory field missing in header program record line." << endl;
		exit(-1);
	}
}
BamReadGroupRecord::BamReadGroupRecord(const string & s) {
	const vector<string> segments = headerLineSplit(s);

	for(vector<string>::const_iterator seg = segments.begin(); seg != segments.end(); seg++) {
		const string tag = seg->substr(0,2);
		const string data = seg->substr(3);

		if(tag == "ID") id = data;
		else if(tag == "CN") sequencingCenter = data;
		else if(tag == "DS") description = data;
		else if(tag == "DT") productionDate = data;
		else if(tag == "FO") flowOrder = data;
		else if(tag == "KS") keySequence = data;
		else if(tag == "LB") library = data;
		else if(tag == "PG") program = data;
		else if(tag == "PI") predictedInsertionSize = data;
		else if(tag == "PL") sequencingTechnology = data;
		else if(tag == "PU") platformUnit = data;
		else if(tag == "SM") sample = data;
	}

	if(id.empty()) {
		cerr << "Mandatory field missing in header read group line." << endl;
		exit(-1);
	}
}

BamHeader::BamHeader (const string & text){
	string line;
	stringstream in(text);
	while(true) {
		string line;
        getline(in, line);

		if(!in.good())
			break;
		if(line[0] != '@') {
			cerr << "Sam header format problem: line doesn't begin with a '@'. Quitting." << endl;
			cerr << "This line: " << line << endl;
			exit(-1);
		}

		if(line[3] != '\t') {
			cerr << "Sam header format problem: line doesn't have a tab after the tag. Quitting." << endl;
			cerr << "This line: " << line << endl;
			exit(-1);
		}

		string tag_in = line.substr(1, 2);
		string data_in = line.substr(4);

		if(tag_in == "CO") {
			co.push_back(data_in);
		} else if (tag_in == "RG") {
			rg.push_back(BamReadGroupRecord(data_in));
		} else if (tag_in == "SQ") {
			sq.push_back(BamSequenceRecord(data_in));
		} else if (tag_in == "PG") {
			pg.push_back(BamProgramRecord(data_in));
		} else if (tag_in == "HD") {
            const vector<string> segments = headerLineSplit(data_in);
            string sort_str;
            
            for(vector<string>::const_iterator seg = segments.begin(); seg != segments.end(); seg++) {
                const string tag = seg->substr(0,2);
                const string data = seg->substr(3);
                
                if(tag == "VN") format_version = data;
                else if(tag == "SO") sort_str = data;
            }
            
            if(sort_str.empty() || format_version.empty()) {
                cerr << "Mandatory field missing in header HD line." << endl;
                exit(-1);
            }
            if(sort_str == "unsorted")
                sort_order = BamHeader::SORT_UNSORTED;
            else if(sort_str == "coordinate")
                sort_order = BamHeader::SORT_COORDINATE;
            else if(sort_str == "queryname")
                sort_order = BamHeader::SORT_QUERYNAME;
            else if(sort_str == "unknown")
                sort_order = BamHeader::SORT_UNKNOWN;
            else {
                cerr << "Unknown sort order '" << sort_str << "'." << endl;
                exit(-1);
            }
		} else {
			cerr << "Sam header format problem: tag '" << tag_in << "' wasn't CO RG SQ PG or HD. Quitting." << endl;
			cerr << "This line: " << line << endl;
			exit(-1);
		}
	}
};

//////////////////////////////////////////
// Printing

const string BamHeader::toString() const {
	stringstream s;

	//print HD
	s << "@HD\tVN:" << format_version << "\tSO:";
	switch(sort_order) {
		case SORT_UNKNOWN: s << "unknown"; break;
		case SORT_UNSORTED: s << "unsorted"; break;
		case SORT_QUERYNAME: s << "queryname"; break;
		case SORT_COORDINATE: s << "coordinate"; break;
	}
	s << "\n";

	//print SQ
	for(BamSequenceRecords::const_iterator i = sq.begin(); i != sq.end(); i++)
		s << i->toString() << "\n";

	//print RG
	for(BamReadGroupRecords::const_iterator i = rg.begin(); i != rg.end(); i++)
		s << i->toString() << "\n";

	//print PG
	for(BamProgramRecords::const_iterator i = pg.begin(); i != pg.end(); i++)
		s << i->toString() << "\n";

	//print CO
	for(vector<string>::const_iterator i = co.begin(); i != co.end(); i++)
		s << "@CO:\t" << *i << "\n";

	return s.str();
}

const string BamSequenceRecord::toString() const {
    stringstream s;
	s << "@SQ\tSN:" << name;
	s << "\tLN:" << length;
	if(!assemblyId.empty())
		s << "\tAS:" << assemblyId;
	if(!checksum.empty())
		s << "\tM5:" << checksum;
	if(!species.empty())
		s << "\tSP:" << species;
	if(!uri.empty())
		s << "\tUR:" << uri;

	return s.str();
}

const string BamReadGroupRecord::toString() const {
	stringstream s;
	s << "@RG\tID:" << id;
	if(!sequencingCenter.empty())
		s << "\tCN:" << sequencingCenter;
	if(!description.empty())
		s << "\tDS:" << description;
	if(!productionDate.empty())
		s << "\tDT:" << productionDate;
	if(!flowOrder.empty())
		s << "\tFO:" << flowOrder;
	if(!keySequence.empty())
		s << "\tKS:" << keySequence;
	if(!keySequence.empty())
		s << "\tKS:" << keySequence;
	if(!library.empty())
		s << "\tLB:" << library;
	if(!program.empty())
		s << "\tPG:" << program;
	if(!predictedInsertionSize.empty())
		s << "\tPI:" << predictedInsertionSize;
	if(!sequencingTechnology.empty())
		s << "\tPL:" << sequencingTechnology;
	if(!platformUnit.empty())
		s << "\tPU:" << platformUnit;
	if(!sample.empty())
		s << "\tSM:" << sample;

	return s.str();
}

const string BamProgramRecord::toString() const {
	stringstream s;

	s << "@PG\tID:" << id;

	if(!name.empty())
		s << "\tPN:" << name;
	if(!commandLine.empty())
		s << "\tCL:" << commandLine;
	if(!previousProgramId.empty())
		s << "\tPP:" << previousProgramId;
	if(!version.empty())
		s << "\tVN:" << version;

	return s.str();
}
