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

#include <api/SamHeader.h>
#include <api/BamAlignment.h>

template <class output_stream_t>
class BamSerializer {
public:
    bool open(const std::string & filename, const BamTools::SamHeader & header);
    void close();
    bool is_open() { return output_stream.is_open(); }

    bool write(const BamTools::BamAlignment & alignment);
    
    // access the real output stream object, in case we need to change some setting,
    // like compression level for bgzfstream.
    output_stream_t & getOutputStream() { return output_stream; }
protected:
    output_stream_t output_stream;
    
    
    void CreatePackedCigar(const std::vector<BamTools::CigarOp>& cigarOperations, std::string & packedCigar);
    void EncodeQuerySequence(const std::string& query, std::string& encodedQuery);

};

template <class output_stream_t>
bool BamSerializer<output_stream_t>::open(const std::string & filename, const BamTools::SamHeader & header) {
    
    output_stream.open(filename.c_str());
    
    if(output_stream.fail()) return false;
    
    //magic header
    output_stream.write("BAM\1", 4);
    
    //header text
    std::string header_txt = header.ToString();
    int header_size = header_txt.size() + 1;
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
void BamSerializer<output_stream_t>::CreatePackedCigar(const std::vector<BamTools::CigarOp>& cigarOperations, std::string& packedCigar) {
    
    // initialize
    const size_t numCigarOperations = cigarOperations.size();
    packedCigar.resize(numCigarOperations * 4);
    
    // pack the cigar data into the string
    unsigned int* pPackedCigar = (unsigned int*)packedCigar.data();
    
    // iterate over cigar operations
    std::vector<BamTools::CigarOp>::const_iterator coIter = cigarOperations.begin();
    std::vector<BamTools::CigarOp>::const_iterator coEnd  = cigarOperations.end();
    for ( ; coIter != coEnd; ++coIter ) {
        
        // store op in packedCigar ("MIDNSHP=X")
        uint8_t cigarOp;
        switch ( coIter->Type ) {
            case 'M': cigarOp = 0; break;
            case 'I': cigarOp = 1; break;
            case 'D': cigarOp = 2; break;
            case 'N': cigarOp = 3; break;
            case 'S': cigarOp = 4; break;
            case 'H': cigarOp = 5; break;
            case 'P': cigarOp = 6; break;
            case '=': cigarOp = 7; break;
            case 'X': cigarOp = 8; break;
            default:
                std::cerr << std::string("BamSerializer: invalid CIGAR operation type ") << coIter->Type << " . Aborting." << std::endl;
                exit(-1);
        }
        
        *pPackedCigar = coIter->Length << 4 | cigarOp;
        pPackedCigar++;
    }
}

// The implementation of this function is originally from BamWriter_p.cpp from bamtools.
// Bamtools is released under the BSD license.
template <class output_stream_t>
void BamSerializer<output_stream_t>::EncodeQuerySequence(const std::string& query, std::string& encodedQuery) {
    
    // prepare the encoded query string
    const size_t queryLength = query.size();
    const size_t encodedQueryLength = static_cast<size_t>((queryLength+1)/2);
    encodedQuery.resize(encodedQueryLength);
    char* pEncodedQuery = (char*)encodedQuery.data();
    const char* pQuery = (const char*)query.data();
    
    // walk through original query sequence, encoding its bases (letter to position) "=ACMGRSVTWYHKDBN"
    unsigned char nucleotideCode;
    bool useHighWord = true;
    while ( *pQuery ) {
        switch ( *pQuery ) {
            case '=': nucleotideCode =  0; break;
            case 'A': nucleotideCode =  1; break;
            case 'C': nucleotideCode =  2; break;
            case 'M': nucleotideCode =  3; break;
            case 'G': nucleotideCode =  4; break;
            case 'R': nucleotideCode =  5; break;
            case 'S': nucleotideCode =  6; break;
            case 'V': nucleotideCode =  7; break;
            case 'T': nucleotideCode =  8; break;
            case 'W': nucleotideCode =  9; break;
            case 'Y': nucleotideCode = 10; break;
            case 'H': nucleotideCode = 11; break;
            case 'K': nucleotideCode = 12; break;
            case 'D': nucleotideCode = 13; break;
            case 'B': nucleotideCode = 14; break;
            case 'N': nucleotideCode = 15; break;
            default:
                std::cerr << std::string("BamSerializer: invalid sequence base: ") << *pQuery << ". Aborting." << std::endl;
                exit(-1);
        }
        
        // pack the nucleotide code
        if ( useHighWord ) {
            *pEncodedQuery = nucleotideCode << 4;
            useHighWord = false;
        } else {
            *pEncodedQuery |= nucleotideCode;
            ++pEncodedQuery;
            useHighWord = true;
        }
        
        // increment the query position
        ++pQuery;
    }
}

// The implementation of this function is originally from BamWriter_p.cpp from bamtools.
// Bamtools is released under the BSD license.
template <class output_stream_t>
bool BamSerializer<output_stream_t>::write(const BamTools::BamAlignment & al) {
    // calculate char lengths
    const unsigned int nameLength         = al.getName().size() + 1;
    const unsigned int numCigarOperations = al.getCigarData().size();
    const unsigned int queryLength        = al.getQueryBases().size();
    const unsigned int tagDataLength      = al.getTagData().size();

    // no way to tell if alignment's bin is already defined (there is no default, invalid value)
    // so we'll go ahead calculate its bin ID before storing
    const uint32_t alignmentBin = CalculateMinimumBin(al.getPosition(), al.GetEndPosition());
    
    // create our packed cigar string
    std::string packedCigar;
    CreatePackedCigar(al.getCigarData(), packedCigar);
    const unsigned int packedCigarLength = packedCigar.size();
    
    // encode the query
    std::string encodedQuery;
    EncodeQuerySequence(al.getQueryBases(), encodedQuery);
    const unsigned int encodedQueryLength = encodedQuery.size();
    
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
    
    // write the query name
    output_stream.write(al.getName().c_str(), nameLength);
    
    output_stream.write(packedCigar.data(), packedCigarLength);
    
    // write the encoded query sequence
    output_stream.write(encodedQuery.data(), encodedQueryLength);
    
    // write the base qualities
    char* pBaseQualities = (char*)al.getQualities().data();
    for ( size_t i = 0; i < queryLength; ++i )
        pBaseQualities[i] -= 33; // FASTQ conversion
    output_stream.write(pBaseQualities, queryLength);
    
    output_stream.write(al.getTagData().data(), tagDataLength);

    return true;
}

typedef BamSerializer<std::ofstream> RawBamWriter;

#endif