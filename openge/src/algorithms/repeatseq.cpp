/*********************************************************************
 *
 * repeatseq.cpp: 
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 20 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial
 * Purpose License. A copy of this license has been provided in
 * the openge/ directory.
 *
 *********************************************************************/
 
#include "repeatseq.h"
#include "../util/fasta_reader.h"
#include "../util/thread_pool.h"
#include "../util/read_stream_reader.h"

#include <numeric>
#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;

//////////////////////
// from structures.cpp

Repeatseq::STRING_GT::STRING_GT(string a, Sequences b, int c, bool d, int e, int f, bool g, double h){
	GT = c;
	reads = Sequences(b);
	print = a;
	paired = d;
	MapQ = e;
	minFlank = f;
	reverse = g;
	avgBQ = h;
}

Repeatseq::STRING_GT::STRING_GT(){
	print = "";
	GT = 0;
	paired = 0;
	MapQ = 0;
	minFlank = 0;
	reverse = 0;
}
//counter struct is used in array for table files:
struct counter {
    int numGT;              //number of repeats that have a GT
    int numRepeats;         //number of total repeats
    double numRepeats2;     //number of repeats with 2 or more reads present
    double tallyC;          //a running sum of C:'s to be later divided by numRepeats
    
    counter();
};

//overloaded for sorting
bool Repeatseq::STRING_GT::operator<(const STRING_GT &other) const {
	if (minFlank >= other.minFlank) return false;
	else return true;
}

Repeatseq::Sequences::Sequences(string a, string b, string c, bool d){
	preSeq = a;
	alignedSeq = b;
	postSeq = c;
	insertions = d;
}

Repeatseq::Sequences::Sequences(){
	preSeq = "";
	alignedSeq = "";
	postSeq = "";
	insertions = 0;
}

struct VCF_INFO {
    std::string chr;
    int start;
    std::string unit;
    int length;
    int purity;
    int depth;
    double confidence;
};

Repeatseq::GT::GT(int rl, int oc, int rev, int minF, double avgbq){
	avgMinFlank = minF;
	readlength = rl;
	occurrences = oc;
	reverse = rev;
	avgBQ = avgbq;
};

counter::counter(){
	numGT = 0;
	numRepeats = 0;
	tallyC = 0;
	numRepeats2 = 0;
}

Repeatseq::Region::Region(string& region) {
	startPos = -1;
	stopPos = -1;
	size_t foundFirstColon = region.find(":");
	// we only have a single string, use the whole sequence as the target
	if (foundFirstColon == string::npos) {
		startSeq = region;
	}
	else {
		startSeq = region.substr(0, foundFirstColon);
		size_t foundRangeDots = region.find("-", foundFirstColon);
		if (foundRangeDots == string::npos) {
			startPos = atoi(region.substr(foundFirstColon + 1).c_str());
			stopPos = startPos; // just print one base if we don't give an end
		} else {
			startPos = atoi(region.substr(foundFirstColon + 1, foundRangeDots - foundRangeDots - 1).c_str());
			stopPos = atoi(region.substr(foundRangeDots + 1).c_str()); // to the start of this chromosome
		}
	}
}

string Repeatseq::Region::toString() const {
    stringstream ret;
    ret << startSeq << ":" << startPos;
    if(startPos != stopPos)
        ret << "-" << stopPos;
    return ret.str();
}

int Repeatseq::Region::length(void) const {
	if (stopPos > 0) {
		return stopPos - startPos + 1;
	} else {
		return 1;
	}
}

/////////////////////
// from repeatseq.cpp

/*
 repeatseq.cpp - main source code file for RepeatSeq
 
 See "repeatseq.h" for function & custom data structure declarations
 
 This .cpp contains functions:
 (1) runInternal() - Iterate line by line through the TRF file (calling
 print_output() on each region).
 
 (2) print_output() - This function is called for each repeat in the repeat file, and
 handles the calling of other functions to determine genotype and print
 data to files.
 
 (3) parseCigar() - Uses CIGAR sequence to align read with reference sequence.
 
 (4) printGenoPerc() - perform statistical analysis to determine most likely genotype and its
 likelihood.
 
 (5) getVCF() - print variant record to VCF file.
 */

#include "repeatseq.h"
#include <algorithm>
#include <unistd.h>

int PHI_TABLE[5][5][5][2] = { { { { 177 , 7 } , { 234 , 15 } , { 991 , 76 } , { 425 , 39 } , { 234 , 28 } } ,
    { { 324 , 23 } , { 587 , 59 } , { 250 , 28 } , { 123 , 17 } , { 764 , 110 } } ,
    { { 268 , 21 } , { 375 , 40 } , { 139 , 20 } , { 66 , 6 } , { 21 , 2 } } ,
    { { 142 , 7 } , { 25 , 1 } , { 6 , 0 } , { 7 , 0 } , { 2 , 0 } } ,
    { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } }
} ,
    { { { 948 , 35 } , { 120 , 5 } , { 458 , 29 } , { 174 , 11 } , { 75 , 7 } } ,
        { { 226 , 10 } , { 405 , 22 } , { 159 , 9 } , { 622 , 49 } , { 296 , 18 } } ,
        { { 512 , 33 } , { 111 , 8 } , { 410 , 32 } , { 170 , 17 } , { 61 , 6 } } ,
        { { 636 , 37 } , { 141 , 7 } , { 48 , 3 } , { 12 , 1 } , { 1 , 0 } } ,
        { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } }
    } ,
    { { { 364 , 11 } , { 471 , 40 } , { 176 , 28 } , { 73 , 7 } , { 33 , 8 } } ,
        { { 351 , 12 } , { 637 , 28 } , { 295 , 14 } , { 116 , 16 } , { 62 , 3 } } ,
        { { 460 , 18 } , { 80 , 0 } , { 39 , 5 } , { 13 , 0 } , { 4 , 0 } } ,
        { { 70 , 1 } , { 5 , 1 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } } ,
        { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } }
    } ,
    { { { 337 , 13 } , { 287 , 34 } , { 102 , 18 } , { 21 , 9 } , { 13 , 4 } } ,
        { { 266 , 14 } , { 277 , 37 } , { 90 , 11 } , { 28 , 5 } , { 17 , 3 } } ,
        { { 886 , 68 } , { 106 , 12 } , { 40 , 1 } , { 17 , 1 } , { 7 , 1 } } ,
        { { 170 , 7 } , { 26 , 1 } , { 13 , 0 } , { 2 , 0 } , { 0 , 0 } } ,
        { { 52 , 2 } , { 8 , 0 } , { 1 , 1 } , { 1 , 0 } , { 0 , 0 } }
    } ,
    { { { 304 , 6 } , { 302 , 24 } , { 105 , 19 } , { 38 , 7 } , { 7 , 7 } } ,
        { { 192 , 10 } , { 172 , 23 } , { 50 , 11 } , { 23 , 3 } , { 9 , 1 } } ,
        { { 300 , 11 } , { 54 , 13 } , { 15 , 4 } , { 1 , 5 } , { 1 , 0 } } ,
        { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } } ,
        { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } }
    }
};

//function declarations:
double fact(int);
double retSumFactOverIndFact(int, int, int);
std::string getVCF(std::string, std::string, std::string, int, char, bool, VCF_INFO);
double PhredToFloat(char);
std::string setToCD (std::string);
void printHeader(std::ofstream&);
void printArguments();

//function to convert phred score to probability score
inline double PhredToFloat(char chr){
	// p_right-base = 1 - 10^(-Q/10)
	double temp = chr - 33;
	return (1 - pow(10,temp/-10));
}

double log_factorial[10000] = {};
string VERSION = "0.6.4";

void Repeatseq::RepeatseqJob::runJob() {
    repeatseq->print_output(region, vcfFile, oFile, callsFile, reads, somatic_reads);
    
    for(vector<OGERead *>::const_iterator i = reads.begin(); i != reads.end(); i++) {
        OGERead::deallocate(*i);
    }
    complete = true;
    repeatseq->flushWrites();
}

bool readOverlapsRegion(const OGERead & read, const BamRegion & region) {
    if(read.getRefID() != region.LeftRefID)
        return false;
    if(read.GetEndPosition() >= region.LeftPosition && read.getPosition() < region.RightPosition)
        return true;
    return false;
}

void Repeatseq::flushWrites() {
    static bool busy = false;
    if(busy) return;    //unsafe, but should cut down on spinlock contention. Not every thread needs to call this.
    busy = true;
    write_flush_lock.lock();

    for(int i = last_deleted_job; i < jobs.size(); i++) {
        if(!jobs[i])
            continue;
        
        if(jobs[i]->complete) {
            //flush
            if (makeRepeatseqFile){ oFile << jobs[i]->oFile.str(); }
            if (makeCallsFile){ callsFile << jobs[i]->callsFile.str(); }
            if (makeVcfFile) vcfFile << jobs[i]->vcfFile.str();
            
            //delete
            delete jobs[i];
            jobs[i]= NULL;
            last_deleted_job = i;
        } else {
            break;  //stop flushing once we find a job that isn't complete
        }
            
    }
    
    write_flush_lock.unlock();
    busy = false;
}

Repeatseq::RegionStringComparator::RegionStringComparator(const BamTools::SamSequenceDictionary & d)
: d(d)
{}

bool Repeatseq::RegionStringComparator::operator() (string a, string b) const {
    const Repeatseq::Region r1(a);
    const Repeatseq::Region r2(b);
    int rid1 = d.IndexOfString(r1.startSeq);
    int rid2 = d.IndexOfString(r2.startSeq);
    if(rid1 != rid2)
        return rid1 < rid2;
    return r1.startPos < r2.startPos;
}

int Repeatseq::runInternal() {

    srand( time(NULL) );
    
    //load log_factorial vector
    //for (int i=1,val=0 ; i < 100000; ++i){ <--- this change breaks the code, val must be floating point
    double val = 0;
    for (int i=1; i < 10000; ++i){ // scaled down to 10k for speed
        val += log(i);
        log_factorial[i] = val;
    }

    //create index filepaths & output filepaths (ensuring output is to current directory):
    string output_filename = setToCD(output_filename_base + paramString + ".repeatseq");
    string calls_filename = setToCD(output_filename_base + paramString + ".calls");
    string vcf_filename = setToCD(output_filename_base + paramString + ".vcf");
    
    //open input & output filestreams:
    if (makeRepeatseqFile){ oFile.open(output_filename.c_str()); }
    if (makeCallsFile){ callsFile.open(calls_filename.c_str()); }
    if (makeVcfFile) vcfFile.open(vcf_filename.c_str());
    ifstream range_file(intervals_filename.c_str());
    if (!range_file.is_open()) { throw "Unable to open input range file."; }
    
    //print VCF header information:
    printHeader(vcfFile);

    fasta_reader.open(fasta_filename);

    //read in the region file
    vector<string> regions;
    string line;
    while(getline(range_file,line))
        regions.push_back(line);
    
    BamTools::SamSequenceDictionary sequence_dictionary = getHeader().Sequences;
    ogeSortMt(regions.begin(), regions.end(), RegionStringComparator(sequence_dictionary));

    int num_jobs = regions.size();
    jobs.reserve(num_jobs);
    
    MultiReader somatic_reader;
    if(!somatic_input.empty())
        somatic_reader.open(somatic_input);
    
    deque<OGERead *> read_buffer, read_buffer_somatic;
    
    //set up threads to actually print the output
    for(int job_id = 0; job_id != num_jobs; job_id++) {
        const Region r(regions[job_id]);
        if(!sequence_dictionary.Contains(r.startSeq)) {
            cerr << "Sequence " << r.startSeq << " from regions file is not in BAM header. Quitting." << endl;
            exit(-1);
        }
        int RID = sequence_dictionary.IndexOfString(r.startSeq);
        const BamRegion btregion(RID, r.startPos, RID, r.stopPos);

        // add reads to a buffer until we are sure we have all reads that might be in this region
        while(true) {
            OGERead * read = getInputAlignment();
            
            if(!read)   //have read the whole file, we are done.
                break;
            if(read->getRefID() < RID || (read->getRefID() == RID && read->GetEndPosition()+1 < r.startPos) ) {
                putOutputAlignment(read);
                continue;
            }
            read_buffer.push_front(read);
            
            if(read->getRefID() > RID || (read->getRefID() == RID && read->getPosition() > r.stopPos+1))
                break;
        }
        
        //remove reads from the left side that shouldn't be in this region
        while(!read_buffer.empty()) {
            OGERead * back = read_buffer.back();
            
            if(back->getRefID() < RID || (back->getRefID() == RID && back->GetEndPosition()+1 < r.startPos) ) {
                read_buffer.pop_back();
                putOutputAlignment(back);
            } else
                break;
        }
        
        if(!somatic_input.empty()) {
            
            // add reads to a buffer until we are sure we have all reads that might be in this region
            while(true) {
                OGERead * read = somatic_reader.read();
                
                if(!read)   //have read the whole file, we are done.
                    break;
                if(read->getRefID() < RID || (read->getRefID() == RID && read->GetEndPosition()+1 < r.startPos) ) {
                    OGERead::deallocate(read);
                    continue;
                }
                read_buffer_somatic.push_front(read);
                
                if(read->getRefID() > RID || (read->getRefID() == RID && read->getPosition() > r.stopPos+1))
                    break;
            }
            
            //remove reads from the left side that shouldn't be in this region
            while(!read_buffer_somatic.empty()) {
                OGERead * back = read_buffer_somatic.back();
                
                if(back->getRefID() < RID || (back->getRefID() == RID && back->GetEndPosition()+1 < r.startPos) ) {
                    read_buffer_somatic.pop_back();
                    OGERead::deallocate(back);
                } else
                    break;
            }
        }
        
        RepeatseqJob * job = new RepeatseqJob(this, job_id, regions[job_id]);
        
        //add all reads that overlap to the buffer
        for(deque<OGERead *>::const_reverse_iterator i = read_buffer.rbegin(); i != read_buffer.rend(); i++) {
            if((*i)->GetEndPosition() >= btregion.LeftPosition && (*i)->getPosition() < btregion.RightPosition) {
            //if( readOverlapsRegion(**i, btregion)) {
                OGERead * read = OGERead::allocate();
                *read = **i;
                job->reads.push_back(read);
            }
        }
        
        if(!somatic_input.empty()) {
            //add all reads that overlap to the buffer
            for(deque<OGERead *>::const_reverse_iterator i = read_buffer_somatic.rbegin(); i != read_buffer_somatic.rend(); i++) {
                if((*i)->GetEndPosition() >= btregion.LeftPosition && (*i)->getPosition() < btregion.RightPosition) {
                    //if( readOverlapsRegion(**i, btregion)) {
                    OGERead * read = OGERead::allocate();
                    *read = **i;
                    job->somatic_reads.push_back(read);
                }
            }
        }

        //start the job
        ThreadPool::sharedPool()->addJob(job);
        jobs.push_back(job);
    }

    if(!somatic_input.empty())
        somatic_reader.close();

    //wait for all workers to finish
    ThreadPool::sharedPool()->waitForJobCompletion();
    
    flushWrites();

    return 0;
}

inline string parseCigar(stringstream &cigarSeq, string &alignedSeq, const string &QS, vector<string> & insertions, const int alignStart, const int refStart, const int LR_CHARS_TO_PRINT, double &avgBQ){
	int reserveSize = alignedSeq.length() + 500;
    
	//reserve sufficient space (so iterators remain valid)
	alignedSeq.reserve(reserveSize);
	string tempInsertions = "";
	tempInsertions.reserve(reserveSize);
	
	//iterators & other variables
	string::iterator it=alignedSeq.begin();
	string::iterator START;
	char cigChar;
	int cigLength;
	bool STARTset = 0;
	int posLeft = refStart - alignStart;
	int posLeftINS = refStart - alignStart - LR_CHARS_TO_PRINT;
	bool firstRun = true;
	
	//determine average base quality:
	avgBQ = 0;
	for (int i=0; i<QS.length(); ++i){ avgBQ += PhredToFloat(QS[i]); }
	avgBQ /= QS.length();
    
	cigarSeq.clear();
	while (!cigarSeq.eof()) {
		//parse:
		cigLength = -1;
		cigarSeq >> cigLength;
		cigarSeq >> cigChar;
		
		//Perform operations on aligned seq:
		switch(cigChar) {
			case 'M':                   //MATCH to the reference
				for (int i = cigLength; i>0; i--) {
					if (posLeft>0) {
						posLeft--;
						posLeftINS--;
					}
					else if (!STARTset) {
						START = it;
						STARTset = 1;
					}
					it++;
				}
				
				break;
				
			case 'I':                  //INSERTION to the reference
				tempInsertions = "";
				*(it-1) += 32;	//convert previous letter to lower case (to mark the following insertion)
				
				for (int i = cigLength; i>0; i--) {
					tempInsertions += *it + 1;
					*it = 'd';         //convert to d's to remove later
					++it;
				}
				if (posLeftINS <= 0) { insertions.push_back(tempInsertions); }
				
				break;
				
			case 'D':                       //DELETION from the reference
				for (int i = cigLength; i>0; i--) {
					alignedSeq.insert(it, 1, '-');	//inserts - into the aligned seq
					posLeft--;
					posLeftINS--;
					if (posLeft < 0 && !STARTset) {
						START = it;
						STARTset = 1;
					}
					++it;
				}
				break;
				
			case 'N':       //SKIPPED region from the reference
				return "";	//fail the read (return null string)
				
			case 'S':                       //SOFT CLIP on the read (clipped sequence present in <seq>)
				if (firstRun && !STARTset) posLeft+=cigLength;
				else if (firstRun) START -= cigLength;
				
				for (int i = cigLength; i>0; i--) {
					if (posLeft>0) {
						posLeft--;
						posLeftINS--;
					}
					else if (!STARTset) {
						START = it;
						STARTset = 1;
					}
					
					*it = 'S';				//mark as soft-clipped
					++it;
				}
				break;
				
			case 'H':   //HARD CLIP on the read (clipped sequence NOT present in <seq>)
				break;
				
			case 'P':   //PADDING (silent deletion from the padded reference sequence)
				if (posLeft>0) {
					posLeft--;
					posLeftINS--;
				}
				else if (!STARTset) {
					START = it;
					STARTset = 1;
				}
				
				it += cigLength;
				break;
		}       //end case
		firstRun = 0;
	}
	
	int offset = alignStart - refStart;
	if (!STARTset) {
		START = alignedSeq.begin();
		while ((offset--) > 0) { alignedSeq.insert(alignedSeq.begin(),1,'x'); }
	}
	int numD = 0;
	for (string::iterator ii = START; ii > START - LR_CHARS_TO_PRINT && ii >= alignedSeq.begin(); --ii) {
		if (*ii == 'd') numD++;
	}
	
	string temp = "";
	temp.reserve(500);
	
	string::iterator ii = START;
	for (int i = 0; i < numD + LR_CHARS_TO_PRINT; ++i) {
		if (ii > alignedSeq.begin()) {
			--ii;
			temp.insert(0,1,*ii);
		}
		else {
			temp.insert(0,1,'x');
		}
	}
	
	ii = START;
	for (int i = 0; i < alignStart - refStart; ++i){ temp += 'x'; }
	while (ii < alignedSeq.end()) { temp += *(ii++); }
	
	return temp; //return modified string
}

std::string Repeatseq::getReference(const Region & target) const {
    
	//ensure target doesn't overrun end of chromosome
	if (target.startPos+target.length() > fasta_reader.getSequenceLength(target.startSeq)+1) throw "Target range is outside of chromosome.\n exiting..";
	
	string sequence;                // holds reference sequence
	//if asked to print entire sequence:
	if (target.startPos == -1) sequence = fasta_reader.readSequence(target.startSeq, 0, fasta_reader.getSequenceLength(target.startSeq));
	
	//when start position is exactly at the beginning of the chromosome:
	else if (target.startPos == 1)
		sequence = " "
		+ fasta_reader.readSequence(target.startSeq, target.startPos - 1, target.length())
		+ " "
		+ fasta_reader.readSequence(target.startSeq, target.startPos - 1 + target.length(), LR_CHARS_TO_PRINT);
	
	//when start position is within 20 of beginning of chromosome:
	else if (target.startPos < 1 + LR_CHARS_TO_PRINT)
		sequence = fasta_reader.readSequence(target.startSeq, 0, target.startPos - 1)
		+ " "
		+ fasta_reader.readSequence(target.startSeq, target.startPos - 1, target.length())
		+ " "
		+ fasta_reader.readSequence(target.startSeq, target.startPos - 1 + target.length(), LR_CHARS_TO_PRINT);
	
	
	//when end position is exactly at the end of chromosome:
	else if (target.startPos+target.length() ==  fasta_reader.getSequenceLength(target.startSeq)+1)
		sequence = fasta_reader.readSequence(target.startSeq, target.startPos - 1 - LR_CHARS_TO_PRINT, LR_CHARS_TO_PRINT)
		+ " "
		+ fasta_reader.readSequence(target.startSeq, target.startPos - 1, target.length())
		+ " ";
	
	//when end position is within 20 of end of chromosome:
	else if (target.startPos+target.length()+LR_CHARS_TO_PRINT > fasta_reader.getSequenceLength(target.startSeq)+1)
		sequence = fasta_reader.readSequence(target.startSeq, target.startPos - 1 - LR_CHARS_TO_PRINT, LR_CHARS_TO_PRINT)
		+ " "
		+ fasta_reader.readSequence(target.startSeq, target.startPos - 1, target.length())
		+ " "
		+ fasta_reader.readSequence(target.startSeq, target.startPos - 1 + target.length(), fasta_reader.getSequenceLength(target.startSeq)-target.startPos-target.length()+1);
	
	//all other cases:
	else sequence = fasta_reader.readSequence(target.startSeq, target.startPos - 1 - LR_CHARS_TO_PRINT, LR_CHARS_TO_PRINT)
		+ " "
		+ fasta_reader.readSequence(target.startSeq, target.startPos - 1, target.length())
		+ " "
		+ fasta_reader.readSequence(target.startSeq, target.startPos - 1 + target.length(), LR_CHARS_TO_PRINT);
    
	std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
    return sequence;
}

std::vector<Repeatseq::STRING_GT> Repeatseq::makeToPrint(const std::vector<OGERead *> reads, const Region & target, const std::string & leftReference, const std::string & centerReference, const std::string & rightReference, int & numStars, int & depth) const {
    
	stringstream ssPrint;                   //where data to print will be stored
	string PreAlignedPost = "";             //contains all 3 strings to be printed
	vector<STRING_GT> toPrint;
	toPrint.reserve(100);
	
	//cout << "trying " << target.startSeq << ":" << target.startPos - 1 << "-" << target.stopPos - 1 << endl;
	// iterate through alignments in this region,
	for(vector<OGERead *>::const_iterator i = reads.begin(); i != reads.end(); i++) {
        
        const OGERead & al = **i;
		//cout << " found\n";
        vector<string> insertions;
		ssPrint.str("");
		stringstream cigarSeq;
		int gtBonus = 0;
		
        const vector<CigarOp> & cigar_data = al.getCigarData();
		if (cigar_data.empty()) {
			numStars++;
			continue;
			//if CIGAR is not there, it's * case..
			//so increment numStars and get next alignment
		}
		
		//load cigarSeq
		for ( vector<CigarOp>::const_iterator it=cigar_data.begin(); it < cigar_data.end(); it++ ) {
			cigarSeq << it->Length;
			cigarSeq << it->Type;
		}
		
		//run parseCigar:
		double avgBQ;
        string query_bases = al.getQueryBases();
        string qualities = al.getQualities();
		string PreAlignedPost = parseCigar(cigarSeq, query_bases, al.getQualities(), insertions, al.getPosition() + 1, target.startPos, LR_CHARS_TO_PRINT, avgBQ);
		if (PreAlignedPost == ""){
			//If an 'N' or other problem was found
			cout << "N found-- Possible Error!\n";
			continue;
		}
		
		//adjust for d's
		for (int a = PreAlignedPost.find('d',0); a!=-1; a=PreAlignedPost.find('d',0)) {
			if ( (a + 1) > LR_CHARS_TO_PRINT && (a + 1) < LR_CHARS_TO_PRINT + target.length()) gtBonus+=1;
			PreAlignedPost.erase(a,1);
		}
		
		//set strings to print based off of value input
		string PreSeq, AlignedSeq, PostSeq;
		
		//if there's not enough characters to make it through PreSeq, skip read
		if (PreAlignedPost.length() < LR_CHARS_TO_PRINT+1) continue;
		
		//Split PreAlignedPost into 3 substrings
		PreSeq = PreAlignedPost.substr(0,LR_CHARS_TO_PRINT);
		AlignedSeq = PreAlignedPost.substr(LR_CHARS_TO_PRINT, target.length());
		if (AlignedSeq.length() < target.length()) AlignedSeq.resize(target.length(),'x');
		else PostSeq = PreAlignedPost.substr(LR_CHARS_TO_PRINT + target.length(), LR_CHARS_TO_PRINT);
		PostSeq.resize(LR_CHARS_TO_PRINT,'x');
		
		if (AlignedSeq[target.length()/2] != 'x') ++depth;      //increment depth (if middle character is NOT an x)
		int numMatchesL = 0, numMatchesR = 0;
		int minflank = 0;
        
		// if first and last characters of sequence range are present in read, print it's information:
		if (AlignedSeq[0] != ' ' && AlignedSeq[0] != 'x' && AlignedSeq[0] != 'X' && AlignedSeq[0] != 'S') {
			if (AlignedSeq[AlignedSeq.length()-1]!= 'x' && AlignedSeq[AlignedSeq.length()-1]!= ' ' && AlignedSeq[AlignedSeq.length()-1]!='X' && AlignedSeq[AlignedSeq.length()-1] != 'S') {
				string toprintPre = string(PreSeq);
				string toprintAligned = string(AlignedSeq);
				string toprintPost = string(PostSeq);
				
				bool hasinsertions = (! insertions.empty());
				if (hasinsertions){
					//PROCESS SEQUENCE:
					//put insertions back in pre-sequence (as lower case) here
					for (int i = 0; i < toprintPre.length();){
						if (toprintPre[i] > 96 && toprintPre[i] != 'x'){	//is lowercase
							toprintPre[i++] -= 32;							//convert to uppercase
							if (i == toprintPre.length()) toprintAligned = insertions.front() + toprintAligned;
							else toprintPre.insert(i,insertions.front());
							insertions.erase(insertions.begin());
						}
						else ++i;
					}
					//put insertions back in Aligned-sequence (as lower case) here
					for (int i = 0; i < toprintAligned.length();){
						if (toprintAligned[i] > 96 && toprintAligned[i] != 'x'){	//is lowercase
							toprintAligned[i++] -= 32;							//convert to uppercase
							if (i == toprintAligned.length()) toprintPost = insertions.front() + toprintPost;
							else toprintAligned.insert(i,insertions.front());
							insertions.erase(insertions.begin());
						}
						else ++i;
					}
					//put insertions back in Post-sequence (as lower case) here
					for (int i = 0; i < toprintPost.length();){
						if (toprintPost[i] > 96 && toprintPost[i] != 'x'){	//is lowercase
							toprintPost[i++] -= 32;							//convert to uppercase
							if (i == toprintPost.length()) toprintPost += insertions.front();
							else toprintPost.insert(i,insertions.front());
							insertions.erase(insertions.begin());
						}
						else ++i;
					}
				}
				
				ssPrint << " " << (al.getPosition() + 1) << " ";   //start position
				
				//Determine & print read size information:
				int readSize = 0;
                const vector<CigarOp> & cigar_data = al.getCigarData();
				for (vector<CigarOp>::const_iterator it=cigar_data.begin(); it < cigar_data.end(); it++){
					if (it->Type == 'M' || it->Type == 'I' || it->Type == 'S' || it->Type == '=' || it->Type == 'X'){
						readSize += it->Length;         //increment readsize by the length
					}
				}
				ssPrint << readSize << " ";      //read size
				
				//FILTER based on min/max read length restrictions:
				if (readLengthMin && readSize < readLengthMin){ continue; }
				if (readLengthMax && readSize > readLengthMax){ continue; }
                
				//Determine consecutive matching flanking bases (LEFT):
				string::const_iterator i = PreSeq.end()-1;
				string::const_iterator i2 = leftReference.end()-1;
				bool consStreak = 1;
				numMatchesL = 0;
				for (int ctr = 0; ctr < PreSeq.length(); ++ctr ) {      //-1 compensates for matching null character @ end of all strings
					if ((*i != *i2) && (*i != *i2 + 32)) {
						consStreak = 0;
						if (ctr < 3){
							if (*i == 'x' || *i == 'S' || (*i2 != '-' && *i == '-') || (*i2 == '-' && *i != '-' )){
								continue; //fail the read
							}
						}
					}
					else if (consStreak){ ++numMatchesL;}
					--i; --i2;
				}
				
				//Determine consecutive matching flanking bases (RIGHT):
				i = PostSeq.begin();
				i2 = rightReference.begin();
				consStreak = 1;
				numMatchesR = 0;
				for (int ctr = 0; ctr < PostSeq.length(); ctr++) {
					if ((*i != *i2) && (*i != *i2 + 32)){
						consStreak = 0;
						if (ctr < 3){
							if (*i == 'x' || *i == 'S' || (*i2 != '-' && *i == '-') || (*i2 == '-' && *i != '-' )){
								continue; //fail the read
							}
						}
					}
					else{
						if (consStreak) ++numMatchesR;
					}
					++i; ++i2;
				}
				
				// Set minflank & print matching # of consecutive bases to the left/right of repeat
				if (numMatchesR < minflank) minflank = numMatchesR;
				else { minflank = numMatchesL; }
				ssPrint << numMatchesL << " " << numMatchesR << " ";
				
				//FILTER based on consecutive flank bases
				if (numMatchesL < consLeftFlank) continue;
				if (numMatchesR < consRightFlank) continue;
				
				//Print avgBQ:
				ssPrint << "B:" << double(int(10000*avgBQ))/10000 << " ";
                
				//FILTER based on MapQ, then print MapQ
				if (al.getMapQuality() < MapQuality) continue;  //MapQuality Filter
				ssPrint << "M:" << al.getMapQuality() << " ";
				
				//PRINT FLAG STRING:
				ssPrint << "F:";
				if (al.IsPaired()) ssPrint << 'p';
				if (al.IsProperPair()) ssPrint << 'P';
				if (!al.IsMapped()) ssPrint << 'u';
				if (!al.IsMateMapped()) ssPrint << 'U';
				if (al.IsReverseStrand()) ssPrint << 'r';
				if (al.IsMateReverseStrand()) ssPrint << 'R';
				if (al.IsFirstMate()) ssPrint << '1';
				if (al.IsSecondMate()) ssPrint << '2';
				if (!al.IsPrimaryAlignment()) ssPrint << 's';
				if (al.IsFailedQC()) ssPrint << 'f';
				if (al.IsDuplicate()) ssPrint << 'd';
				
				//print CIGAR string:
				ssPrint << " C:";
				for (vector<CigarOp>::const_iterator it=cigar_data.begin(); it < cigar_data.end(); it++) {
					ssPrint << it->Length;
					ssPrint << it->Type;
				}
				
				//-MULTI filter (check for XT:A:R tag):
				string stringXT;
				al.GetTag("XT",stringXT);
				if (multi && stringXT.find('R',0) != -1) continue;  //if stringXT contains R, ignore read
				
				//-PP filter (check if read is properly paired):
				if (properlyPaired && !al.IsProperPair()){ continue; }
				
				ssPrint << " ID:" << al.getName() << endl;
				
				toPrint.push_back( STRING_GT(ssPrint.str(), Sequences(toprintPre, toprintAligned, toprintPost, hasinsertions), AlignedSeq.length() + gtBonus, al.IsProperPair(), al.getMapQuality(), minflank, al.IsReverseStrand(), avgBQ) );
			}
		}        //end if statements
		
	} //end while loop
    
    
	
	//push reference sequences into vectors for expansion & printing:
	toPrint.insert( toPrint.begin(), STRING_GT("\n", Sequences(leftReference, centerReference, rightReference, 0), 0, 0, 0, 0, 0, 0.0) );
	
	// If any of the reads have insertions, expand the reads without inserted bases so all reads are fully printed:
	if (true){
		//Handle PRE-SEQ:
		for (int index = 0, limit = LR_CHARS_TO_PRINT + 1; index < limit; ++index){
			for (vector<STRING_GT>::iterator jt=toPrint.begin(); jt < toPrint.end(); jt++){
				if (index >= jt->reads.preSeq.length()) continue;
				if (jt->reads.preSeq[index] == 'B' || jt->reads.preSeq[index] == 'U' || jt->reads.preSeq[index] == 'D' || jt->reads.preSeq[index] == 'H' || jt->reads.preSeq[index] == 'O'){
					limit++;
					
					vector<STRING_GT>::iterator pt;
					if (jt + 1 == toPrint.end()) pt = toPrint.begin();
					else pt = jt+1;
					
					while(pt != jt){
						if (pt->reads.preSeq[index] == 'B' || pt->reads.preSeq[index] == 'U' || pt->reads.preSeq[index] == 'D' || pt->reads.preSeq[index] == 'H' || pt->reads.preSeq[index] == 'O') pt->reads.preSeq[index] -=1;
						else pt->reads.preSeq.insert(index,"-");
						
						if (pt + 1 == toPrint.end()) pt = toPrint.begin();
						else pt++;
					}
					
					jt->reads.preSeq[index] -= 1;
					
				}
			}
		}
		//Handle ALIGNED-SEQ:
		for (int index = 0, limit = target.length() + 1; index < limit; ++index){
			for (vector<STRING_GT>::iterator jt=toPrint.begin(); jt < toPrint.end(); jt++){
				if (index >= jt->reads.alignedSeq.length()) continue;
				if (jt->reads.alignedSeq[index] == 'B' || jt->reads.alignedSeq[index] == 'U' || jt->reads.alignedSeq[index] == 'D' || jt->reads.alignedSeq[index] == 'H' || jt->reads.alignedSeq[index] == 'O'){
					limit++;
					
					vector<STRING_GT>::iterator pt;
					if (jt + 1 == toPrint.end()) pt = toPrint.begin();
					else pt = jt+1;
					
					while(pt != jt){
						if (pt->reads.alignedSeq[index] == 'B' || pt->reads.alignedSeq[index] == 'U' || pt->reads.alignedSeq[index] == 'D' || pt->reads.alignedSeq[index] == 'H' || pt->reads.alignedSeq[index] == 'O') pt->reads.alignedSeq[index] -=1;
						else pt->reads.alignedSeq.insert(index,"-");
						
						if (pt + 1 == toPrint.end()) pt = toPrint.begin();
						else pt++;
					}
					
					jt->reads.alignedSeq[index] -= 1;
					
				}
			}
		}
		//Handle POST-SEQ:
		for (int index = 0, limit = LR_CHARS_TO_PRINT + 1; index < limit; ++index){
			for (vector<STRING_GT>::iterator jt=toPrint.begin(); jt < toPrint.end(); jt++){
				if (index >= jt->reads.postSeq.length()) continue;
				if (jt->reads.postSeq[index] == 'B' || jt->reads.postSeq[index] == 'U' || jt->reads.postSeq[index] == 'D' || jt->reads.postSeq[index] == 'H' || jt->reads.postSeq[index] == 'O'){
					limit++;
					
					vector<STRING_GT>::iterator pt;
					if (jt + 1 == toPrint.end()) pt = toPrint.begin();
					else pt = jt+1;
					
					while(pt != jt){
						if (pt->reads.postSeq[index] == 'B' || pt->reads.postSeq[index] == 'U' || pt->reads.postSeq[index] == 'D' || pt->reads.postSeq[index] == 'H' || pt->reads.postSeq[index] == 'O') pt->reads.postSeq[index] -=1;
						else pt->reads.postSeq.insert(index,"-");
						
						if (pt + 1 == toPrint.end()) pt = toPrint.begin();
						else pt++;
					}
					
					jt->reads.postSeq[index] -= 1;
					
				}
			}
		}
		
		// fix for insertions/deletions immediately following repeat:
		int index = 0;
		while(toPrint.begin()->reads.postSeq[index] == '-'){ ++index; }
		for (vector<STRING_GT>::iterator jt=toPrint.begin(); jt < toPrint.end(); jt++){
			jt->reads.alignedSeq += jt->reads.postSeq.substr(0, index);
			jt->reads.postSeq.erase(0,index);
			
			if (jt->GT){ //if it's not the reference..
				string repeat = jt->reads.alignedSeq;
				repeat.erase(std::remove (repeat.begin(), repeat.end(), '-'), repeat.end());
				
				jt->GT = repeat.length();
			}
		}
	}
    
    return toPrint;
}

vector<Repeatseq::GT> Repeatseq::makeGenotypeInformation(vector<STRING_GT> & toPrint) const {
    vector<GT> vectorGT;
	vectorGT.reserve(100);
	// Build VectorGT from toPrint:
	for (vector<STRING_GT>::iterator tP=toPrint.begin(); tP < toPrint.end(); ++tP) {
		if (tP->GT == 0) continue; //ignore reference
		if (!vectorGT.empty()) {
			for (vector<GT>::iterator it = vectorGT.begin(); 1; ) {
				if (it->readlength == tP->GT) {
					it->occurrences += 1;
					it->avgBQ += tP->avgBQ;
					it->avgMinFlank += tP->minFlank;
					if (tP->reverse) it->reverse += 1;
					break;
				}
				else {
					++it;
					if (it == vectorGT.end()) {
						GT a(tP->GT, 1, tP->reverse, tP->minFlank, tP->avgBQ);
						vectorGT.insert(vectorGT.end(), a);
						break;
					}
				}
			}
		}
		else {
			GT a(tP->GT, 1, tP->reverse, tP->minFlank, tP->avgBQ);
			vectorGT.insert(vectorGT.end(), a);
		}
	}
    
	//sort vectorGT by occurrences..
	sort(vectorGT.begin(), vectorGT.end());
    
	//average out BQs & flanks
	for (vector<GT>::iterator it = vectorGT.begin(); it < vectorGT.end(); ++it) {
		it->avgBQ /= it->occurrences;
		it->avgMinFlank /= it->occurrences;
	}
    
    return vectorGT;
}

inline void Repeatseq::print_output(const string & region_line, stringstream &vcf_buffer,  stringstream &o_buffer, stringstream &calls_buffer, const vector<OGERead *> & reads, const vector<OGERead *> & reads_somatic) const {

	// parse region argument:
	string secondColumn = region_line.substr(region_line.find('\t',0)+1,-1); // text string to the right of tab
	if (secondColumn == "") cout << "missing information after the tab in region file for " << region_line << ".\ncontinuing..." << endl;
	string region = region_line.substr(0,region_line.find('\t',0));          //erases all of region string after tab
	
	// parse secondColumn:
	if (int(secondColumn.find('_',0)) == -1) {
		cout << "improper second column found for " << region << ".\ncontinuing with next region..." << endl;
		return;
	}
	int unitLength = atoi(secondColumn.substr(0,secondColumn.find('_',0)).c_str());
	string UnitSeq = secondColumn.substr(secondColumn.rfind('_')+1);
	
	int pos = 0;
	for (int i = 0; i < 3; ++i) pos = secondColumn.find('_',pos + 1);
	++pos; //increment past fourth '_'
	double purity = atof(secondColumn.substr(pos,secondColumn.find('_',pos)).c_str());
	
	Region target(region);
	if (target.startPos > target.stopPos) throw "Invalid input file...";
	
    string sequence = getReference(target);
	
	int firstSpace = sequence.find(' ',0);
	int secondSpace = sequence.find(' ',firstSpace+1);
	
	string leftReference, centerReference, rightReference;
	if (firstSpace != 0) leftReference = sequence.substr(0,firstSpace);
	else leftReference = "";
	centerReference = sequence.substr(firstSpace+1,secondSpace-firstSpace-1);
	if (secondSpace != -1) rightReference = sequence.substr(secondSpace+1,-1);
	else rightReference = "";

	// define our region of interest:
	// debug-cout << "region: " << target.startSeq << ":" << target.startPos-1 << "-" << target.stopPos-1 << endl;
	//BamRegion bamRegion(reader.GetReferenceID(target.startSeq), target.startPos - 1,reader.GetReferenceID(target.startSeq), target.stopPos - 1);
	//reader.SetRegion(bamRegion);
	
	// prep for getting alignment info
	int depth = 0;
	int numStars = 0;
	int depth_somatic = 0;
	int numStars_somatic = 0;

	vector<STRING_GT> toPrint = makeToPrint(reads, target, leftReference, centerReference, rightReference, numStars, depth);
	vector<GT> vectorGT = makeGenotypeInformation(toPrint);
    
	vector<STRING_GT> toPrint_somatic;
	vector<GT> vectorGT_somatic;
    
    if(error_model == ERROR_SOMATIC_A || error_model == ERROR_SOMATIC_B) {
        toPrint_somatic = makeToPrint(reads_somatic, target, leftReference, centerReference, rightReference, numStars_somatic, depth_somatic);
        vectorGT_somatic = makeGenotypeInformation(toPrint);
    }

    {
        int numReads = toPrint.size() - 1;  // subtract 1 for the reference.

        double avgMapQ;
        {
            int occ = 0;
            int tallyMapQ = 0;
            for (vector<STRING_GT>::iterator it=toPrint.begin(); it < toPrint.end(); ++it) {
                tallyMapQ = tallyMapQ + it->MapQ;
                ++occ;
            }
            if (occ) avgMapQ = double(tallyMapQ)/double(occ);
            else avgMapQ = -1;
        }

        int majGT = 0;
        int occurMajGT = 0;
        int totalOccurrences = 0;
        double concordance = 0;
        
        //output header line
        o_buffer << "~" << region << " ";
        o_buffer << secondColumn;
        o_buffer << " REF:" << target.length();
        o_buffer << " A:";
        if (!vectorGT.size()) {
            o_buffer << "NA ";
            concordance = -1;
            majGT = -1;
        }
        else {
            if (vectorGT.size() == 1) {
                if (numReads == 1) {
                    concordance = -1;
                    majGT = vectorGT.begin()->readlength;
                    o_buffer << "NA ";
                }
                else {
                    concordance = 1;
                    majGT = vectorGT.begin()->readlength;
                    o_buffer << vectorGT.begin()->readlength << " ";
                    //oFile << "(" << vectorGT.begin()->avgBQ << ")"; //temp
                }
            }
            else {
                for (vector<GT>::iterator it=vectorGT.begin(); it < vectorGT.end(); it++) {
                    o_buffer << it->readlength << "[" << it->occurrences << "] " ;
                    //oFile << "(" << it->avgBQ << ") ";
                    if (it->occurrences >= occurMajGT) {
                        occurMajGT = it->occurrences;
                        if (it->readlength > majGT) majGT = it->readlength;
                    }
                    totalOccurrences += it->occurrences;
                }
                concordance = double(double(occurMajGT)-1.00) / double(double(totalOccurrences)-1.00);
            }
        }
        
        //concordance = # of reads that support the majority GT / total number of reads
        if (concordance < 0) o_buffer << "C:NA";
        else o_buffer << "C:" << concordance;
        
        o_buffer << " D:" << depth << " R:" << numReads << " S:" << numStars;
        if (avgMapQ >= 0) o_buffer << " M:" << double(int(100*avgMapQ))/100;
        else o_buffer << " M:NA";
        
        o_buffer << " GT:";
        calls_buffer << region << "\t" << secondColumn << "\t";
        vector<int> vGT;
        double conf = 0;
        if (vectorGT.size()  == 0) {
            o_buffer << "NA L:NA" << endl;
            calls_buffer << "NA\tNA\n";
        }
        else if (vectorGT.size() > 9){          // if more than 9 GTs are present
            o_buffer << "NA L:NA" << endl;
            calls_buffer << "NA\tNA\n";
        }
        else if (concordance >= 0.99){          //no need to compute confidence if all the reads agree
            o_buffer << majGT << " L:50" << endl;
            calls_buffer << majGT << "L:50" << endl;
        }
        else {
            if(error_model == ERROR_SOMATIC_B)
                vGT = somaticConfidence(vectorGT, vectorGT_somatic, target, error_model, conf);
            else
                vGT = printGenoPerc(vectorGT, target.length(), unitLength, conf, mode);

            if (numReads <= 1){ conf = 0; }
            //write genotypes to calls & repeats file
            if (vGT.size() == 0) { throw "vGT.size() == 0.. ERROR!\n"; }
            else if (vGT.size() == 1 && conf > 3.02) { o_buffer << vGT[0] << " L:" << conf << "\n"; calls_buffer << vGT[0] << '\t' << conf << '\n'; }
            else if (vGT.size() == 2 && conf > 3.02) { o_buffer << vGT[0] << "h" << vGT[1] << " L:" << conf << "\n"; calls_buffer << vGT[0] << "h" << vGT[1] << '\t' << conf << '\n'; }
            else{ o_buffer << "NA L:" << conf << endl; calls_buffer << "NA\tNA\n"; }
        }
        
        // Set info for printing VCF file
        //structure to be passed to VCF-writing function:
        VCF_INFO INFO;
        INFO.chr = target.startSeq;
        INFO.start = target.startPos + 1;
        INFO.unit = UnitSeq;
        INFO.length = target.length();
        INFO.purity = purity;
        INFO.depth = numReads;
        INFO.confidence = conf;
        
        // GO THROUGH VECTOR AND PRINT ALL REMAINING
        if (toPrint.size()>1){ //if there are reads present..
            string REF = toPrint[0].reads.alignedSeq;
            bool homo = false;
            if (vGT.size() == 1) homo = true;
            
            for (vector<STRING_GT>::iterator it=toPrint.begin(); it < toPrint.end(); it++) {
                // print .repeats file:
                o_buffer << it->reads.preSeq << " " << it->reads.alignedSeq << " " << it->reads.postSeq << it->print;
                
                // finished printing to .repeats file.
                if (vGT.size() != 0 && conf > 3.02){
                    if (vGT.size() > 1 || vGT[0] != target.length() /*there's been a mutation*/){
                        // print .vcf file:
                        vector<int>::iterator tempgt = std::find(vGT.begin(), vGT.end(), it->GT);
                        if (tempgt != vGT.end() && it->GT != target.length()){
                            //debug vcf << "VCF record for " << REF << " --> " << it->reads.alignedSeq << "..\n";
                            
                            // the read represents one of our genotypes..
                            string vcfRecord = getVCF(it->reads.alignedSeq, REF, target.startSeq, target.startPos, *(leftReference.end()-1), homo, INFO);
                            vcf_buffer << vcfRecord;
                            
                            //remove the genotype from the genotype list..
                            vGT.erase( tempgt );
                        }
                        // finished printing to .vcf file.
                    }
                }
                
                // continue iterating through each read..
            }
        }
    }
    
    if(error_model == ERROR_SOMATIC_A || error_model == ERROR_SOMATIC_B)
    {
        int numReads = toPrint.size() - 1;  // subtract 1 for the reference.
        
        double avgMapQ;
        {
            int occ = 0;
            int tallyMapQ = 0;
            for (vector<STRING_GT>::iterator it=toPrint.begin(); it < toPrint.end(); ++it) {
                tallyMapQ = tallyMapQ + it->MapQ;
                ++occ;
            }
            if (occ) avgMapQ = double(tallyMapQ)/double(occ);
            else avgMapQ = -1;
        }
        
        int majGT = 0;
        int occurMajGT = 0;
        int totalOccurrences = 0;
        double concordance = 0;
        
        //output header line
        o_buffer << "~~" << region << " ";
        o_buffer << secondColumn;
        o_buffer << " REF:" << target.length();
        o_buffer << " A:";
        if (!vectorGT.size()) {
            o_buffer << "NA ";
            concordance = -1;
            majGT = -1;
        }
        else {
            if (vectorGT.size() == 1) {
                if (numReads == 1) {
                    concordance = -1;
                    majGT = vectorGT.begin()->readlength;
                    o_buffer << "NA ";
                }
                else {
                    concordance = 1;
                    majGT = vectorGT.begin()->readlength;
                    o_buffer << vectorGT.begin()->readlength << " ";
                    //oFile << "(" << vectorGT.begin()->avgBQ << ")"; //temp
                }
            }
            else {
                for (vector<GT>::iterator it=vectorGT.begin(); it < vectorGT.end(); it++) {
                    o_buffer << it->readlength << "[" << it->occurrences << "] " ;
                    //oFile << "(" << it->avgBQ << ") ";
                    if (it->occurrences >= occurMajGT) {
                        occurMajGT = it->occurrences;
                        if (it->readlength > majGT) majGT = it->readlength;
                    }
                    totalOccurrences += it->occurrences;
                }
                concordance = double(double(occurMajGT)-1.00) / double(double(totalOccurrences)-1.00);
            }
        }
        
        //concordance = # of reads that support the majority GT / total number of reads
        if (concordance < 0) o_buffer << "C:NA";
        else o_buffer << "C:" << concordance;
        
        o_buffer << " D:" << depth << " R:" << numReads << " S:" << numStars;
        if (avgMapQ >= 0) o_buffer << " M:" << double(int(100*avgMapQ))/100;
        else o_buffer << " M:NA";
        
        o_buffer << " GT:";
        calls_buffer << region << "\t" << secondColumn << "\t";
        vector<int> vGT;
        double conf = 0;
        if (vectorGT.size()  == 0) {
            o_buffer << "NA L:NA" << endl;
            calls_buffer << "NA\tNA\n";
        }
        else if (vectorGT.size() > 9){          // if more than 9 GTs are present
            o_buffer << "NA L:NA" << endl;
            calls_buffer << "NA\tNA\n";
        }
        else if (concordance >= 0.99){          //no need to compute confidence if all the reads agree
            o_buffer << majGT << " L:50" << endl;
            calls_buffer << majGT << "L:50" << endl;
        }
        else {
            vGT = somaticConfidence(vectorGT_somatic, vectorGT, target, error_model, conf);
            if (numReads <= 1){ conf = 0; }
            //write genotypes to calls & repeats file
            if (vGT.size() == 0) { throw "vGT.size() == 0.. ERROR!\n"; }
            else if (vGT.size() == 1 && conf > 3.02) { o_buffer << vGT[0] << " L:" << conf << "\n"; calls_buffer << vGT[0] << '\t' << conf << '\n'; }
            else if (vGT.size() == 2 && conf > 3.02) { o_buffer << vGT[0] << "h" << vGT[1] << " L:" << conf << "\n"; calls_buffer << vGT[0] << "h" << vGT[1] << '\t' << conf << '\n'; }
            else{ o_buffer << "NA L:" << conf << endl; calls_buffer << "NA\tNA\n"; }
        }
        
        // Set info for printing VCF file
        //structure to be passed to VCF-writing function:
        VCF_INFO INFO;
        INFO.chr = target.startSeq;
        INFO.start = target.startPos + 1;
        INFO.unit = UnitSeq;
        INFO.length = target.length();
        INFO.purity = purity;
        INFO.depth = numReads;
        INFO.confidence = conf;
        
        // GO THROUGH VECTOR AND PRINT ALL REMAINING
        if (toPrint.size()>1){ //if there are reads present..
            string REF = toPrint[0].reads.alignedSeq;
            bool homo = false;
            if (vGT.size() == 1) homo = true;
            
            for (vector<STRING_GT>::iterator it=toPrint.begin(); it < toPrint.end(); it++) {
                // print .repeats file:
                o_buffer << it->reads.preSeq << " " << it->reads.alignedSeq << " " << it->reads.postSeq << it->print;
                
                // finished printing to .repeats file.
                if (vGT.size() != 0 && conf > 3.02){
                    if (vGT.size() > 1 || vGT[0] != target.length() /*there's been a mutation*/){
                        // print .vcf file:
                        vector<int>::iterator tempgt = std::find(vGT.begin(), vGT.end(), it->GT);
                        if (tempgt != vGT.end() && it->GT != target.length()){
                            //debug vcf << "VCF record for " << REF << " --> " << it->reads.alignedSeq << "..\n";
                            
                            // the read represents one of our genotypes..
                            string vcfRecord = getVCF(it->reads.alignedSeq, REF, target.startSeq, target.startPos, *(leftReference.end()-1), homo, INFO);
                            vcf_buffer << vcfRecord;
                            
                            //remove the genotype from the genotype list..
                            vGT.erase( tempgt );
                        }
                        // finished printing to .vcf file.
                    }
                }
                
                // continue iterating through each read..
            }
        }
    }
	
	return;
}

inline double fact ( int n ){
    double fact = 1;
    while (n > 1) fact *= n--;
    return fact;
}

inline int nCr (int n, int r){
    return fact(n)/fact(r)/fact(n-r);
}

class tagAndRead{
public:
    string m_name;
    double m_pX;
    tagAndRead(string a, double b){
        m_name = a;
        m_pX = b;
    }
};

inline bool compareTAR(tagAndRead a, tagAndRead b){
    return (a.m_pX > b.m_pX);
}

inline double retBetaMult(int* vector, int alleles){
	double value = 1, sum = 0;
    // alleles + 1 --> 2 if homozygous, 3 if hetero
    for (int i = 0; i < alleles + 1; ++i) {
		value += log_factorial[vector[i]-1];
		sum += vector[i];
	}
    value -= log_factorial[int(sum) - 1];
	return value;
	
}

inline double multinomial_beta(const vector<double> & alpha) {
    double numerator = 1.;
    for(vector<double>::const_iterator i = alpha.begin(); i != alpha.end(); i++)
        numerator *= tgamma(*i);
    
    double alpha_sum = accumulate(alpha.begin(),alpha.end(),0);
    
    return numerator / tgamma(alpha_sum);
}

inline double dirichlet(const vector<double> & alpha, const vector<double> & x) {
    assert(alpha.size() == x.size());
    
    double prod = 1.;
    
    for(int k = 0; k < alpha.size(); k++)
        prod *= pow(x[k], alpha[k] - 1);
    
    return prod / multinomial_beta(alpha);
}

struct somatic_caller_data {
    string name;
    double p_x, p_x_gi, p_gi_x, p_gi;
    vector<double> alpha;
    vector<double> x;
    somatic_caller_data(string name)
    : name(name)
    , p_x(0), p_x_gi(0), p_gi_x(0), p_gi(0)
    {}
    
    bool operator<(const somatic_caller_data & d) const {
        return this->p_gi_x > d.p_gi_x;
    }
};

inline vector<int> Repeatseq::somaticConfidence(vector<GT> & vectorGT, const vector<GT> & vectorGT_reference, const Region & target, error_model_t model, double &confidence) const {
    
    int total_x = 0, total_alpha = 0;
    map<int,int> x_counts, alpha_counts;    //map of count# to occurrences
    
    for(vector<GT>::const_iterator i = vectorGT.begin(); i != vectorGT.end(); i++) {
        total_x += i->occurrences;
        x_counts[i->readlength] = i->occurrences;
    }
    for(vector<GT>::const_iterator i = vectorGT_reference.begin(); i != vectorGT_reference.end(); i++) {
        total_alpha += i->occurrences;
        alpha_counts[i->readlength] = i->occurrences;
    }

    vector<GT> vectorGT_union;
    set_union(vectorGT.begin(), vectorGT.end(), vectorGT_reference.begin(), vectorGT_reference.end(), back_inserter(vectorGT_union));
    
    set<GT> allGTs(vectorGT_union.begin(), vectorGT_union.end());
    
    int n = total_x;

	vector<somatic_caller_data> pXarray;
	allGTs.insert(Repeatseq::GT(0,0,0,0,0.0)); //allows locus to be considered homozygous

    
    // STEP 2
    for (set<GT>::const_iterator it = allGTs.begin(); it !=  allGTs.end(); ++it){
        set<GT>::const_iterator itnext = it;
        itnext++;
        for (set<GT>::const_iterator jt = itnext; jt !=  allGTs.end(); ++jt){
            set<GT>::const_iterator jtnext = jt;
            if(jt->occurrences != 0 || jt->readlength != 0)
                jtnext++;
            for (set<GT>::const_iterator kt = jtnext; kt !=  allGTs.end(); ++kt){
                int alleles = 1;
                int alpha_error = total_alpha;
                int x_error = total_x;

                vector<double> alpha, x;
                
                stringstream tempss;
                tempss << it->readlength;
                if (jt->occurrences != 0) {
                    tempss << "h" << jt->readlength;
                    alleles++;
                }
                if (kt->occurrences != 0) {
                    tempss << "h" << kt->readlength;
                    alleles++;
                }
                
                cerr << tempss.str() << ":";
                
                ///////////////////////////
                // build alpha and X

                if(it->readlength != 0 || it->occurrences != 0) {
                    alpha.push_back(alpha_counts.count(it->readlength) ? it->occurrences : 0);
                    x.push_back(x_counts.count(it->readlength) ? x_counts[it->readlength] : 0);
                }
                
                if(jt->readlength != 0 || jt->occurrences != 0) {
                    alpha.push_back(alpha_counts.count(jt->readlength) ? jt->occurrences : 0);
                    x.push_back(x_counts.count(jt->readlength) ? x_counts[jt->readlength] : 0);
                }
                
                if(kt->readlength != 0 || kt->occurrences != 0) {
                    alpha.push_back(alpha_counts.count(kt->readlength) ? kt->occurrences : 0);
                    x.push_back(x_counts.count(kt->readlength) ? x_counts[kt->readlength] : 0);
                }
                
                for(vector<double>::const_iterator i = x.begin(); i != x.end(); i++)
                    x_error -= *i;
                for(vector<double>::const_iterator i = alpha.begin(); i != alpha.end(); i++)
                    alpha_error -= *i;

                x.push_back(x_error);
                alpha.push_back(alpha_error);

                assert(x.size() == alpha.size());
                
                pXarray.push_back(somatic_caller_data(tempss.str()));
                pXarray.back().alpha = alpha;
                pXarray.back().x = x;
                
                cerr << "alpha = [";
                for(vector<double>::const_iterator i = alpha.begin(); i != alpha.end(); i++)
                    cerr << " " << *i;
                cerr << " ] ";
                cerr << "x = [";
                for(vector<double>::const_iterator i = x.begin(); i != x.end(); i++)
                    cerr << " " << *i;
                cerr << " ] ";
                
                cerr << endl;
            }
        }
    }

    //STEP 3
    for(vector<somatic_caller_data>::iterator i = pXarray.begin(); i != pXarray.end(); i++) {
        int k = i->alpha.size();

        const vector<double> & alpha = i->alpha;
        const vector<double> & x = i->x;
        
        //calculate p(gi), described under eq 3
        i->p_gi = 1. / (double)k;
        
        //equation 2
        vector<double> alpha_x_1(k);    //alpha + x + 1
        for( int j = 0; j != k; j++)
            alpha_x_1[j] = alpha[j] + x[j] + 1;
        
        vector<double> alpha_1(k);  //alpha + 1
        for( int j = 0; j != k; j++)
            alpha_1[j] = alpha[j] + 1;
        
        double prod_x_fact = 1; // product from j = 1..k of (xj)!
        for( int j = 0; j != k; j++)
            prod_x_fact *= fact(x[j]);
        
        i->p_x_gi = (multinomial_beta(alpha_x_1) / multinomial_beta(alpha_1)) * (fact(n) / prod_x_fact);
    }

    // sum from j = 1..k of pi(x|gj)*pi(gj)
    double sum_pxgj_pgj = 0;
    for(vector<somatic_caller_data>::const_iterator j = pXarray.begin(); j != pXarray.end(); j++)
        sum_pxgj_pgj += j->p_x_gi * j->p_gi;

    for(vector<somatic_caller_data>::iterator i = pXarray.begin(); i != pXarray.end(); i++) {
        //equation 3
        i->p_gi_x = (i->p_x_gi * i->p_gi) / sum_pxgj_pgj;
    }

    sort(pXarray.begin(), pXarray.end());
    for(vector<somatic_caller_data>::const_iterator i = pXarray.begin(); i != pXarray.end(); i++)
        cerr << "P(" << i->name << "): " << i->p_gi_x << endl;

    //build return vector
    vector<int> gts;
    gts.push_back( atoi(pXarray.begin()->name.c_str()) );

    while(true) {
        int hpos = pXarray.begin()->name.find('h');
        if (hpos != -1){
            gts.push_back( atoi(pXarray.begin()->name.c_str()) );
        }
        else
            return gts;
    }
}

inline vector<int> Repeatseq::printGenoPerc(vector<GT> vectorGT, int ref_length, int unit_size, double &confidence, int mode) const {
	if (ref_length > 70) ref_length = 70;
	if (unit_size > 5) unit_size = 5;
	else if (unit_size < 1) unit_size = 1;
	for (vector<GT>::iterator it = vectorGT.begin(); it < vectorGT.end(); ++it){
		it->avgBQ = -30*log10(it->avgBQ);
		if (it->avgBQ < 0){ it->avgBQ = 0; }
		else if (it->avgBQ > 4){ it->avgBQ = 4; }
	}
    
	vector<tagAndRead> pXarray;
	vector<int> gts;
	stringstream toReturn;
    
	vectorGT.push_back(Repeatseq::GT(0,0,0,0,0.0)); //allows locus to be considered homozygous
    double pXtotal = 0;
    string name;
    
	// Calculate LOCAL_PHI
    int mostCommon = 0, secondCommon = 0; double totalSum = 0;
	for (vector<Repeatseq::GT>::iterator it = vectorGT.begin(); it < vectorGT.end(); ++it){
        int tempOccur = it->occurrences;
		if (tempOccur > mostCommon) {secondCommon = mostCommon; mostCommon = tempOccur;}
        else if (tempOccur > secondCommon) {secondCommon = tempOccur;}
        totalSum += tempOccur;
    }
	double LOCAL_PHI;
	if (mode == 1){ LOCAL_PHI = double(totalSum-mostCommon)/totalSum; }
    else if (mode == 2){ LOCAL_PHI = double(totalSum-(mostCommon+secondCommon))/totalSum; }
	else{ LOCAL_PHI = 0; }
	
    for (vector<GT>::iterator it = vectorGT.begin(); it < vectorGT.end(); ++it){
        for (vector<GT>::iterator jt = it+1; jt < vectorGT.end(); ++jt){
            int alleles = 1, errorOccurrences = 0;
            for (vector<GT>::iterator errt = vectorGT.begin(); errt < vectorGT.end(); ++errt){
                if (errt != jt && errt != it) { errorOccurrences += errt->occurrences; }
            }
            
            stringstream tempss;
            if (jt->occurrences != 0) {
                tempss << it->readlength << "h" << jt->readlength;
                alleles = 2;
                name = tempss.str();
            }
            else {
                tempss << it->readlength;
                name = tempss.str();
            }
            
            //if haploid mode is enabled, ensure only one allele is present:
            if (mode == 1 && alleles == 2) continue;
			
            //determine likelihood:
            int* ERROR_TABLE_1 = PHI_TABLE[unit_size-1][ref_length/15][int(it->avgBQ)];
            int* ERROR_TABLE_2 = PHI_TABLE[unit_size-1][ref_length/15][int(jt->avgBQ)];
            
            int ERROR_1[2];
            int ERROR_2[2];
            
            // Set temporary error rate arrays
            if (it->occurrences == 0){ ERROR_1[0] = 0; ERROR_1[1] = 0;}
            else { ERROR_1[0] = ERROR_TABLE_1[1]; ERROR_1[1] = ERROR_TABLE_1[0];}
            if (jt->occurrences == 0){ ERROR_2[0] = 0; ERROR_2[1] = 0;}
            else { ERROR_2[0] = ERROR_TABLE_2[1]; ERROR_2[1] = ERROR_TABLE_2[0];}
            
            // Test if LOCAL error rate is greater than error rate from table
            //  use the greater of the two
            //if (!manualErrorRate){
            //	if (LOCAL_PHI > float(ERROR_1[0])/(ERROR_1[0] + ERROR_1[1])) {
            //		ERROR_1[0] = LOCAL_PHI*int(totalSum);
            //		ERROR_1[1] = int(totalSum) - ERROR_1[0];
            //	}
            //	if (LOCAL_PHI > float(ERROR_2[0])/(ERROR_2[0] + ERROR_2[1])) {
            //		ERROR_2[0] = LOCAL_PHI*int(totalSum);
            //		ERROR_2[1] = int(totalSum) - ERROR_2[0];
            //	}
            //}
            
            int v_numerator[3];
            v_numerator[0] = 1 + ERROR_1[1] + it->occurrences;
            v_numerator[1] = 1 + ERROR_2[1] + jt->occurrences;
            if (alleles == 2){
                v_numerator[2] = 1 + ERROR_1[0] + ERROR_2[0] + errorOccurrences;
            }
            else {
                v_numerator[1] = 1 + ERROR_1[0] + ERROR_2[0] + errorOccurrences;
                v_numerator[2] = -1;
            }
            
            int v_denom[3];
            v_denom[0] = 1 + ERROR_1[1];
            v_denom[1] = 1 + ERROR_2[1];
            if (alleles == 2){
                v_denom[2] = 1 + ERROR_1[0] + ERROR_2[0];
            }
            else {
                v_denom[1] = 1 + ERROR_1[0] + ERROR_2[0];
                v_denom[2] = -1;
            }
            
            
            // Calculate NUMERATOR & DENOMINATOR from arrays
	    	double NUMERATOR = retBetaMult(v_numerator, alleles);
            double DENOM = retBetaMult(v_denom, alleles);
            
            //add genotype & likelihood to pXarray:
            tagAndRead temp = tagAndRead(name, exp(log(retSumFactOverIndFact(it->occurrences,jt->occurrences,errorOccurrences))+NUMERATOR-DENOM));
            pXarray.push_back(temp);
            pXtotal += temp.m_pX;
        }
    }
	
    for (vector<tagAndRead>::iterator it = pXarray.begin(); it < pXarray.end(); ++it){
		it->m_pX /= pXtotal;
	}
    
	// sort, based on likelihood
	sort(pXarray.begin(), pXarray.end(), compareTAR);
	
	// set gts, based on sorted pXarray
	int hpos = pXarray.begin()->m_name.find('h');
	if (hpos == -1){
		//homozygous..
		gts.push_back( atoi(pXarray.begin()->m_name.c_str()) );
	}
	else {
		//heterozygous...
		gts.push_back( atoi(pXarray.begin()->m_name.substr(0, hpos).c_str()) );
		gts.push_back( atoi(pXarray.begin()->m_name.substr(hpos+1, -1).c_str()) );
	}
    
    cerr << "pX: ";
    
    for (vector<tagAndRead>::iterator it = pXarray.begin(); it < pXarray.end(); ++it) {
        cerr << it->m_name << "(" << it->m_pX << ") ";
    }
    cerr << endl;
	
	// set confidence value
	confidence = -10*log10(1-pXarray.begin()->m_pX);
	//cout << confidence << endl << endl;
	if (confidence > 50) confidence = 50; //impose our upper bound to confidence..
	
	//check for NaN --> set to 0
	if (confidence != confidence) {	confidence = 0;	}
    
	return gts;
}

// retSumFactOverIndFact(a,b,c) returns fact(a+b+c))/(fact(a)*fact(b)*fact(c))
// (avoids overflow better by avoiding directly computing fact(a+b+c)..)
double retSumFactOverIndFact(int a, int b, int c){
	double val = 1;
	
	// set max/min values
	int max = a, min1 = b, min2 = c;
	if (b > max && b > c){
		max = b;
		min1 = a;
		min2 = c;
	}
	else if (c > max) {
		max = c;
		min1 = a;
		min2 = b;
	}
	
	// determine value to return..
	for (int i = 1; i <= min1; ++i){
		++max;
		val *= max;
		val /= i;
	}
	for (int i = 1; i <= min2; ++i){
		++max;
		val *= max;
		val /= i;
	}
	
	return val;
}


string getVCF(string alignment, string reference, string chr, int start, char precBase, bool homozygous, VCF_INFO info){
	stringstream vcf;
	int begin, end = -1, bothInsOffset = 0;
	
	//assumes alignment & reference are the same length..
	for(int index = alignment.length()-1; index >= 0; --index){
		if ((alignment[index] != reference[index])){
			if (alignment[index] == '-' || reference[index] == '-'){
				end = index;
				while (((alignment[index] != reference[index]) || alignment[index] == '-') && index >= 0){ index -= 1; }
				begin = index;
				while (index >= 0){
					if (reference[index] == '-' && alignment[index] != '-'){ begin -= 1; }
					else if (reference[index] == '-' && alignment[index] == '-'){ bothInsOffset += 1; }
					index -= 1;
				}
				break;
			}
		}
	}
	
	if (end == -1) { return ""; } // no difference was found..
	
	start += begin;
	if (begin == -1) {
		reference = precBase + reference.substr(0, end+1);
		alignment = precBase + alignment.substr(0, end+1);
	}
	else{
		//call getVCF to recursively get VCF for any indels earlier in the sequences...
		vcf << getVCF(alignment.substr(0,begin), reference.substr(0,begin), chr, start-begin, precBase, homozygous, info);
		
		reference = reference.substr(begin, end-begin+1);
		alignment = alignment.substr(begin, end-begin+1);
	}
	
	//remove -'s
	reference.erase( std::remove(reference.begin(), reference.end(), '-'), reference.end() );
	alignment.erase( std::remove(alignment.begin(), alignment.end(), '-'), alignment.end() );
	
	vcf << chr << '\t';
	vcf << start-bothInsOffset << '\t';
	vcf << "." << '\t'; //ID
	vcf << reference << '\t';
	vcf << alignment << '\t';
	vcf << info.confidence << '\t'; //qual -> -10log_10 prob(call in ALT is wrong)
	if (info.confidence > 0.8) vcf << "PASS\t"; //filter
	else vcf << ".\t";
	vcf << "RU=" << info.unit << ";DP=" << info.depth << ";RL=" << info.length << "\t"; //info
	vcf << "GT\t"; //format
	if (homozygous){ vcf << "1/1\n"; }
	else{ vcf << "1/0\n"; }
	
	return vcf.str();
}

//function to ensure filepath is in the current directory
string setToCD (string filepath){
	if (filepath.rfind('/') != -1){ filepath = filepath.substr( filepath.rfind('/') + 1, -1); }
	return filepath;
}

void printHeader(ofstream &vcf){
	vcf << "##fileformat=VCFv4.0" << endl;
	vcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	vcf << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << endl;
	vcf << "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeated Unit\">" << endl;
	vcf << "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Reference Length of Microsatellite\">" << endl;
	vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" << endl;
}
