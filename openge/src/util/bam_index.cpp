/*********************************************************************
 *
 * bam_index.cpp: Generate an index for a BAM file.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 28 Dec 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial
 * Purpose License. A copy of this license has been provided in
 * the openge/ directory.
 *
 *********************************************************************/

#include "bam_index.h"
#include "bgzf_output_stream.h"
#include "bgzf_input_stream.h"

using namespace std;

// Samtools merges chunks that are at least this close together: see bam_index.c:41 in samtools
const int BAM_MIN_CHUNK_GAP = 32768;

BamIndex::BamIndexSequence::BamIndexBin::BamIndexBin(uint64_t unmapped_reads, uint64_t mapped_reads, uint64_t data_start, uint64_t data_stop) {
    chunks.push_back(pair<uint64_t, uint64_t>(data_start, data_stop));
    chunks.push_back(pair<uint64_t, uint64_t>(mapped_reads, unmapped_reads));
}

void BamIndex::BamIndexSequence::BamIndexBin::addRead(int start_pos, uint64_t file_start, uint64_t file_stop) {
	// there are four cases for addRead-
	// * the read could be directly after a previous chunk, so we should append the chunk
	// * the read could be directly before another chunk, so we should append the chunk
	// * the read could bridge the gap between two chunks
	// * the read could be adjacent to nothing, and we should create a new chunk
	//
	// in practice, since we know we will see reads in sort order, and more importantly, in
	// the order they are in the file, we know that only cases #1 and #4 above are possible
    //
    // samtools merges chunks that are less than BAM_MIN_CHUNK_GAP away from each other; we do the same
    //
    
	// case #4
    if(!chunks.empty() && file_start - chunks.back().second < BAM_MIN_CHUNK_GAP)
        chunks.back().second = file_stop;
    else // case #1
		chunks.push_back(pair<uint64_t, uint64_t>(file_start, file_stop));
}
void BamIndex::BamIndexSequence::BamIndexBin::read(std::ifstream & stream) {
    int32_t size;
    stream.read((char *) &size, sizeof(size));
    for(int i = 0; i < size; i++) {
        uint64_t first, second;
        stream.read((char *)&first, sizeof(first));
        stream.read((char *)&second, sizeof(second));
        chunks.push_back(pair<uint64_t,uint64_t>(first, second));
    }
}

void BamIndex::BamIndexSequence::BamIndexBin::write(ofstream & stream) const {
	int32_t size = chunks.size();
	stream.write((const char *)&size, sizeof(int32_t));
    
	for(vector<pair<uint64_t, uint64_t> >::const_iterator i = chunks.begin(); i != chunks.end(); i++) {
		stream.write((const char *)&i->first, sizeof(uint64_t));
		stream.write((const char *)&i->second, sizeof(uint64_t));
	}
}

void BamIndex::BamIndexSequence::BamIndexBin::remap(BgzfOutputStream * remapper_stream) {
    for(vector<pair<uint64_t, uint64_t> >::iterator i = chunks.begin(); i != chunks.end(); i++)
        *i = pair<uint64_t, uint64_t>(remapper_stream->mapWriteLocationToBgzfPosition(i->first),remapper_stream->mapWriteLocationToBgzfPosition(i->second));
    
    // coalesce chunks, again to match samtools behaviour
    /*bool changed = false;
    do {
        changed = false;
        for(vector<pair<uint64_t, uint64_t> >::iterator i = chunks.begin(); i != chunks.end(); i++) {
            if((i + 1) == chunks.end())
                break;
            
            if(((i+1)->first >> 16) != (i->second >> 16) && ( ((i+1)->first >> 16) - (i->second >> 16)) <= BAM_MIN_CHUNK_GAP) {
                i->second = (i+1)->second;
                cerr << "removing chunk " << ((i+1)->first >> 16) << " " << (i->second >> 16) << endl;
                chunks.erase(i+1);
                changed = true;
                break;
            }
        }
    } while(changed);*/
}

void BamIndex::BamIndexSequence::fillMissing() {
    
    // fill in leading zeros
    vector<uint64_t>::iterator first = lower_bound(linear_index.begin(), linear_index.end(), 0);
    //fill(linear_index.begin(), first, *first);
    
    //fill in blocks of zeros where necessary
    for(vector<uint64_t>::iterator i = first; i != (linear_index.end()-1); i++) {
        vector<uint64_t>::iterator j = i + 1;
        if(*j == 0)
            *j = *i;
    }
}

BamIndex::BamIndexSequence::BamIndexSequence(const BamSequenceRecord & record) {
}

void BamIndex::BamIndexSequence::setMetadataFrame(uint64_t unmapped_reads, uint64_t mapped_reads, uint64_t data_start, uint64_t data_stop) {
    bins[37450] = new BamIndexBin(unmapped_reads, mapped_reads, data_start, data_stop);
}

void BamIndex::BamIndexSequence::read(std::ifstream & stream) {
    int32_t num_bins;
    stream.read((char *) & num_bins, sizeof(num_bins));
    
    for(int i = 0; i < num_bins; i++) {
        uint32_t k;
        BamIndexBin * v = new BamIndexBin();
        stream.read((char *) &k, sizeof(k));
        bins[k] = v;
        v->read(stream);
    }
    
    int32_t intervals;
    stream.read((char *) &intervals, sizeof(intervals));
    linear_index.resize(intervals);
    stream.read((char *) &linear_index[0], sizeof(linear_index[0]) * intervals);
}

void BamIndex::BamIndexSequence::write(ofstream & stream) const {
	int32_t num_bins = bins.size();
    stream.write((const char *) & num_bins, sizeof(num_bins));
    
	for(map<uint32_t,BamIndexBin *>::const_iterator i = bins.begin(); i != bins.end(); i++) {
		uint32_t bin = i->first;
		stream.write((const char *)&bin, sizeof(bin));
		i->second->write(stream);
	}
    
    vector<uint64_t>::const_iterator beginning_of_zeros = linear_index.end();
    while(beginning_of_zeros != linear_index.begin() && *(beginning_of_zeros - 1) == 0)
        beginning_of_zeros--;

	int32_t intervals = distance(linear_index.begin(), beginning_of_zeros);
	stream.write((const char *)&intervals, sizeof(intervals));

	for(vector<uint64_t>::const_iterator i = linear_index.begin(); i != beginning_of_zeros; i++)
		stream.write((const char *) &*i, sizeof(*i));
}

void BamIndex::BamIndexSequence::addRead(const OGERead * read, int end_pos, int bin, uint64_t file_start, uint64_t file_stop) {
	//linear index
    if(read->IsMapped()) {
        int ix_start = read->getPosition() / (1 << 14);
        int ix_end = (end_pos-1) / (1 << 14);
        
        if(linear_index.size() <= ix_end) {
            linear_index.resize(ix_end +1, 0);
        }
        for(int ix = ix_start; ix <= ix_end; ++ix) {
            if(linear_index[ix] == 0)
                linear_index[ix] = file_start;
            else
                linear_index[ix] = min(linear_index[ix], file_start);
        }
        fillMissing();
        
        assert(ix_start < linear_index.size());
        assert(ix_end < linear_index.size());
    }
    
	//normal index
	if(bins.end() == bins.find(bin)) {
		bins[bin] = new BamIndexBin();
	}
    
	bins[bin]->addRead(read->getPosition(), file_start, file_stop);
}

void BamIndex::BamIndexSequence::remap(BgzfOutputStream * remapper_stream) {
    for(vector<uint64_t>::iterator i = linear_index.begin(); i != linear_index.end(); i++)
        *i = remapper_stream->mapWriteLocationToBgzfPosition(*i);
    
	for(map<uint32_t,BamIndexBin *>::iterator i = bins.begin(); i != bins.end(); i++)
        i->second->remap(remapper_stream);
}

BamIndex::BamIndex(const BamHeader & h)
: metadata(h.getSequences().size())
, num_coordless_reads(0)
{
    const BamSequenceRecords sequence_records = h.getSequences();
	for(BamSequenceRecords::const_iterator i = sequence_records.begin(); i != sequence_records.end(); i++)
		sequences.push_back(new BamIndexSequence(*i));
}

BamIndex::~BamIndex() {
	for(vector<BamIndexSequence *>::const_iterator i = sequences.begin(); i != sequences.end(); i++)
		delete *i;
}

void BamIndex::addRead(const OGERead * read, int end_pos, int bin, uint64_t file_start, uint64_t file_stop) {
    
    if(read->getRefID() != -1) {
        metadata_t & m = metadata[read->getRefID()];
        read->IsMapped() ? m.num_mapped_reads++ : m.num_unmapped_reads++;
        
        m.read_start_position = min(m.read_start_position, file_start);
        m.read_stop_position = max(m.read_stop_position, file_stop);
    }
    
    if(read->getPosition() == -1)
        num_coordless_reads++;

	assert(read->getRefID() < sequences.size() || read->getRefID() == -1);
    if(read->getRefID() != -1 && read->getPosition() != -1)
        sequences[read->getRefID()]->addRead(read, end_pos, bin, file_start, file_stop);
}

void BamIndex::readFile(const std::string & filename, BgzfInputStream * remapper_stream) {
    assert(remapper_stream == NULL);    //remapping not supported
    ifstream f;
    f.open(filename.c_str());
    
    if(f.fail())
		cerr << "Warning: failed to open BAM index file " << filename << "." << endl;
    
    char magic[5] = {0};
    f.read(magic, 4);
    uint32_t seq_ct;
    f.read((char *) &seq_ct, sizeof(seq_ct));
    
    for(int i = 0; i < seq_ct; i++) {
        sequences[i]->read(f);
    }
    
    f.read((char *)&num_coordless_reads, sizeof(num_coordless_reads));
    
    f.close();
}

void BamIndex::writeFile(const string & filename, BgzfOutputStream * remapper_stream) const {
    return;
	ofstream f;
	f.open(filename.c_str());
    
	if(f.fail())
		cerr << "Warning: failed to open BAM index file " << filename << "." << endl;
    
	const char * magic = "BAI\1";
	f.write(magic, 4);
    
	uint32_t sequence_ct = sequences.size();
	f.write((const char *)&sequence_ct, sizeof(sequence_ct));
    
    if(remapper_stream)
        for(vector<BamIndexSequence *>::const_iterator i = sequences.begin(); i != sequences.end(); i++)
            (*i)->remap(remapper_stream);

    // some metadata fields can't be remapped as they arent file offsets, so we have to do this manually
    for(int i = 0; i < sequences.size(); i++) {
        const metadata_t & m = metadata[i];
        if(m.num_mapped_reads + m.num_unmapped_reads != 0) {
            if(remapper_stream)
                sequences[i]->setMetadataFrame(m.num_unmapped_reads, m.num_mapped_reads, remapper_stream->mapWriteLocationToBgzfPosition(m.read_start_position), remapper_stream->mapWriteLocationToBgzfPosition(m.read_stop_position));
            else
                sequences[i]->setMetadataFrame(m.num_unmapped_reads, m.num_mapped_reads, m.read_start_position, m.read_stop_position);
        }
    }

    for(vector<BamIndexSequence *>::const_iterator i = sequences.begin(); i != sequences.end(); i++)
        (*i)->write(f);
    
    f.write((const char *)&num_coordless_reads, sizeof(num_coordless_reads));
    
    f.close();
}
