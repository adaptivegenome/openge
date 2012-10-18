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

template <class output_stream_t>
class BamSerializer : public ReadStreamWriter {
public:
    virtual bool open(const std::string & filename, const BamTools::SamHeader & header);
    virtual void close();
    virtual bool is_open() const { return output_stream.is_open(); }

    virtual bool write(const OGERead & alignment);
    
    // access the real output stream object, in case we need to change some setting,
    // like compression level for bgzfstream.
    output_stream_t & getOutputStream() { return output_stream; }
protected:
    output_stream_t output_stream;
};

template <class output_stream_t>
bool BamSerializer<output_stream_t>::open(const std::string & filename, const BamTools::SamHeader & header) {
    
    output_stream.open(filename.c_str());
    
    if(output_stream.fail()) return false;
    
    //magic header
    output_stream.write("BAM\1", 4);
    
    //header text
    std::string header_txt = header.ToString();
    int header_size = header_txt.size();
    char zero = 0;
    output_stream.write((char *)&header_size, 4);
    output_stream.write(header_txt.c_str(), header_size);
    
    //references
    int seq_size = header.Sequences.Size();
    output_stream.write((char *) &seq_size, sizeof(seq_size));
    for(BamTools::SamSequenceConstIterator i = header.Sequences.Begin(); i != header.Sequences.End(); i++) {
        int name_size = i->Name.size() + 1;
        int length = atoi(i->Length.c_str());
        output_stream.write((char *)&name_size, 4);
        output_stream.write(i->Name.c_str(), i->Name.size());
        output_stream.write(&zero,1);
        output_stream.write((char *)&length,4);
    }

    return !output_stream.fail() ;
}

template <class output_stream_t>
void BamSerializer<output_stream_t>::close() {
    output_stream.close();
}

// calculates minimum bin for a BAM alignment interval [begin, end)
// Taken from BAM specification.
uint32_t inline CalculateMinimumBin(const int begin, int end) {
    --end;
    if ( (begin >> 14) == (end >> 14) ) return 4681 + (begin >> 14);
    if ( (begin >> 17) == (end >> 17) ) return  585 + (begin >> 17);
    if ( (begin >> 20) == (end >> 20) ) return   73 + (begin >> 20);
    if ( (begin >> 23) == (end >> 23) ) return    9 + (begin >> 23);
    if ( (begin >> 26) == (end >> 26) ) return    1 + (begin >> 26);
    return 0;
}

// The implementation of this function is originally from BamWriter_p.cpp from bamtools.
// Bamtools is released under the BSD license.
template <class output_stream_t>
bool BamSerializer<output_stream_t>::write(const OGERead & al) {
    // calculate char lengths
    const unsigned int nameLength         = al.getName().size() + 1;
    const unsigned int numCigarOperations = al.getCigarData().size();
    const unsigned int queryLength        = al.getQueryBases().size();
    const unsigned int tagDataLength      = al.getTagData().size();

    // no way to tell if alignment's bin is already defined (there is no default, invalid value)
    // so we'll go ahead calculate its bin ID before storing
    const uint32_t alignmentBin = CalculateMinimumBin(al.getPosition(), al.GetEndPosition());

    // create our packed cigar string
    const unsigned int packedCigarLength = al.getCigarDataLength();
    
    // encode the query
    const unsigned int encodedQueryLength = al.getQueryBasesLength();
    
    // write the block size
    const unsigned int dataBlockSize = nameLength +
    packedCigarLength +
    encodedQueryLength +
    queryLength +
    tagDataLength;
    unsigned int blockSize = 32 + dataBlockSize;
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

    return true;
}

typedef BamSerializer<std::ofstream> RawBamWriter;

#endif