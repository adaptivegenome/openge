#ifndef OGE_BAM_SERIALIZER_H
#define OGE_BAM_SERIALIZER_H

/*********************************************************************
 *
 * bam_serializer.h: Write data to an uncompressed bam stream.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 28 August 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial
 * Purpose License. A copy of this license has been provided in
 * the openge/ directory.
 *
 *********************************************************************/

#include "read_stream_writer.h"
#include "bgzf_output_stream.h"
#include "bam_index.h"

template <class output_stream_t>
class BamSerializer : public ReadStreamWriter {
public:
    BamSerializer(bool generate_index = false) : generate_index(generate_index), index(NULL) {}
    virtual bool open(const std::string & filename, const BamHeader & header);
    virtual void close();
    virtual bool is_open() const { return output_stream.is_open(); }

    virtual bool write(const OGERead & alignment);
    
    // access the real output stream object, in case we need to change some setting,
    // like compression level for bgzfstream.
    output_stream_t & getOutputStream() { return output_stream; }
protected:
    bool generate_index;
    std::string filename;
    output_stream_t output_stream;
    BamIndex * index;
    size_t write_offset;
};

template <class output_stream_t>
bool BamSerializer<output_stream_t>::open(const std::string & filename, const BamHeader & header) {
    this->filename = filename;
    output_stream.open(filename.c_str());
    
    if(output_stream.fail()) return false;
    
    //magic header
    output_stream.write("BAM\1", 4);
    
    //header text
    std::string header_txt = header.toString();
    int header_size = header_txt.size();
    char zero = 0;
    output_stream.write((char *)&header_size, 4);
    output_stream.write(header_txt.c_str(), header_size);
    
    //references
    int seq_size = header.getSequences().size();
    output_stream.write((char *) &seq_size, sizeof(seq_size));
    write_offset = 4 + 4 + header_size + sizeof(seq_size);
    for(BamSequenceRecords::const_iterator i = header.getSequences().begin(); i != header.getSequences().end(); i++) {
        int name_size = i->getName().size() + 1;
        int length = i->getLength();
        output_stream.write((char *)&name_size, 4);
        output_stream.write(i->getName().c_str(), i->getName().size());
        output_stream.write(&zero,1);
        output_stream.write((char *)&length,4);
        write_offset += 4 + 1 + i->getName().size() + 4;
    }
    
    if(generate_index && header.getSortOrder() == BamHeader::SORT_COORDINATE && NULL != dynamic_cast<BgzfOutputStream *>(&output_stream)) {
        index = new BamIndex(header);
    }

    return !output_stream.fail() ;
}

template <class output_stream_t>
void BamSerializer<output_stream_t>::close() {
    output_stream.close();
    if(index)
        index->writeFile(filename + ".bai", dynamic_cast<BgzfOutputStream *>(&output_stream));
}

// calculates minimum bin for a BAM alignment interval [begin, end)
// Taken from BAM specification.
uint32_t inline CalculateMinimumBin(const int beg, int end) {
    
    --end;
    if ( (beg >> 14) == (end >> 14) ) return 4681 + (beg >> 14);
    if ( (beg >> 17) == (end >> 17) ) return  585 + (beg >> 17);
    if ( (beg >> 20) == (end >> 20) ) return   73 + (beg >> 20);
    if ( (beg >> 23) == (end >> 23) ) return    9 + (beg >> 23);
    if ( (beg >> 26) == (end >> 26) ) return    1 + (beg >> 26);
    return 0;
}

// The implementation of this function is originally from BamWriter_p.cpp from bamtools.
// Bamtools is released under the BSD license.
template <class output_stream_t>
bool BamSerializer<output_stream_t>::write(const OGERead & al) {
    // calculate char lengths
    const unsigned int nameLength         = al.getNameLength();
    const unsigned int numCigarOperations = al.getNumCigarOps();
    const unsigned int queryLength        = al.getLength();
    
    const unsigned int end_pos = al.GetEndPosition();

    // no way to tell if alignment's bin is already defined (there is no default, invalid value)
    // so we'll go ahead calculate its bin ID before storing
    const uint32_t alignmentBin = CalculateMinimumBin(al.getPosition(), end_pos);

    // write the block size
    unsigned int blockSize = al.getSupportData().getBlockLength();
    output_stream.write((char*)&blockSize, 4);
    
    // assign the BAM core data
    uint32_t buffer[8];
    buffer[0] = al.getRefID();
    buffer[1] = al.getPosition();
    buffer[2] = (alignmentBin << 16) | (al.getMapQuality() << 8) | nameLength;
    buffer[3] = (al.getAlignmentFlag() << 16) | numCigarOperations;
    buffer[4] = queryLength;
    buffer[5] = al.getMateRefID();
    buffer[6] = al.getMatePosition();
    buffer[7] = al.getInsertSize();
    
    // write the BAM core
    output_stream.write((char*)&buffer, 32);
    
    const std::string & char_data = al.getSupportData().getAllCharData();

    output_stream.write(&char_data[0] , char_data.size());
    
    if(index) {
        size_t write_len = char_data.size() + 32 + 4;
        index->addRead(&al, end_pos, alignmentBin, write_offset, write_offset + write_len);
        write_offset += write_len;
    }

    return true;
}

typedef BamSerializer<std::ofstream> RawBamWriter;

#endif