#ifndef OGE_BAM_DESERIALIZER_H
#define OGE_BAM_DESERIALIZER_H

/*********************************************************************
 *
 * bam_deserializer.h: Read data from an uncompressed bam stream.
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

#include "read_stream_reader.h"

template <class input_stream_t>
class BamDeserializer : public ReadStreamReader {
public:
    virtual bool open(const std::string & filename);
    virtual const BamTools::SamHeader & getHeader() const { return header; };
    virtual void close();
    virtual BamTools::BamAlignment * read();
    virtual bool is_open() { return input_stream.is_open(); }
protected:
    input_stream_t input_stream;
    BamTools::SamHeader header;
};

template <class input_stream_t>
bool BamDeserializer<input_stream_t>::open(const std::string & filename) {
    input_stream.open(filename.c_str());
    
    if(input_stream.fail()) {
        std::cerr << "Open failed." << std::endl;
        return false;
    }
    
    //check magic values
    char magic[4];
    input_stream.read((char *)&magic, 4);
    
    if(input_stream.fail() || magic[0] != 'B' || magic[1] != 'A' || magic[2] != 'M' || magic[3] != 1) {
        std::cerr << "Error reading BAM stream header magic bytes. Aborting." << std::endl;
        exit(-1);
    }
    
    //read header text
    int text_len;
    input_stream.read((char *)&text_len, 4);
    
    if(input_stream.fail()) {
        std::cerr << "Error reading BAM stream header text length. Aborting." << std::endl;
        exit(-1);
    }

    std::string text(text_len+1, (char)0);
    input_stream.read(&text[0], text_len);
    text[text_len] = 0;

    if(input_stream.fail()) {
        std::cerr << "Error reading BAM stream header text. Aborting." << std::endl;
        exit(-1);
    }
    
    if(input_stream.fail() || magic[0] != 'B' || magic[1] != 'A' || magic[2] != 'M' || magic[3] != 1) {
        std::cerr << "Error reading BAM stream header magic bytes. Aborting." << std::endl;
        exit(-1);
    }
    
    header.SetHeaderText(text);
    
    //read in references, checking values as we go
    
    int reference_ct;
    input_stream.read((char *) &reference_ct, 4);
    
    if(input_stream.fail()) {
        std::cerr << "Error reading BAM stream header text length. Aborting." << std::endl;
        exit(-1);
    }

    if(header.Sequences.Size() != reference_ct) {
        std::cerr << "WARNING: BAM header text sequence data count doesn't match reference sequence list. Is this file corrupted?" << std::endl;
    }

    
    for(int i = 0; i < reference_ct; i++) {
        //read name length
        int name_length;
        input_stream.read((char *)&name_length, 4);
        if(input_stream.fail()) {
            std::cerr << "Error reading BAM stream reference sequence name length. Aborting." << std::endl;
            exit(-1);
        }
        
        //read name
        std::string name(name_length-1, (char)0);
        input_stream.read(&name[0], name_length);
        name[name_length-1] = (char)0;
        if(input_stream.fail()) {
            std::cerr << "Error reading BAM stream reference sequence. Aborting." << std::endl;
            exit(-1);
        }
        
        //read length
        int length;
        input_stream.read((char *)&length, 4);
        if(input_stream.fail()) {
            std::cerr << "Error reading BAM stream header reference sequence length. Aborting." << std::endl;
            exit(-1);
        }

        BamTools::SamSequence & seq = header.Sequences[i];
        
        int seq_length = atoi(seq.Length.c_str());
        if(seq.Name.compare(name) != 0 || seq_length != length) {
            std::cerr << "WARNING: BAM header text doesn't match sequence information. Is this file corrupted?" << std::endl;
            exit(-1);
        }
    }
    
    return true;
}

template <class input_stream_t>
void BamDeserializer<input_stream_t>::close() {
    input_stream.close();
}

template <class input_stream_t>
BamTools::BamAlignment * BamDeserializer<input_stream_t>::read() {
    BamTools::BamAlignment * al = BamTools::BamAlignment::allocate();

    // read in the 'block length' value, make sure it's not zero
    char buffer[sizeof(uint32_t)];

    input_stream.read(buffer, sizeof(uint32_t));
    if(input_stream.eof())
        return NULL;

    if ( input_stream.fail() ) {
        std::cerr << "Expected more bytes reading BAM core. Is this file truncated or corrupted?" << std::endl;
        delete al;
        return NULL;
    }

    BamTools::BamAlignment::BamAlignmentSupportData & SupportData = al->getSupportData();
    SupportData.BlockLength = BamTools::UnpackUnsignedInt(buffer);

    if ( SupportData.BlockLength == 0 )
        if ( SupportData.BlockLength == 0 )
            return NULL;
    
    // read in core alignment data, make sure the right size of data was read
    char x[32];
    input_stream.read(x, 32);
    if ( input_stream.fail() ) {
        std::cerr << "Expected more bytes reading BAM core. Is this file truncated or corrupted?" << std::endl;
        delete al;
        return NULL;
    }

    // set BamAlignment 'core' and 'support' data
    al->setRefID(BamTools::UnpackSignedInt(&x[0]));
    al->setPosition(BamTools::UnpackSignedInt(&x[4]));
    
    unsigned int tempValue = BamTools::UnpackUnsignedInt(&x[8]);
    al->setBin(tempValue >> 16);
    al->setMapQuality(tempValue >> 8 & 0xff);
    SupportData.QueryNameLength = tempValue & 0xff;
    
    tempValue = BamTools::UnpackUnsignedInt(&x[12]);
    al->setAlignmentFlag(tempValue >> 16);
    SupportData.NumCigarOperations = tempValue & 0xffff;
    
    SupportData.QuerySequenceLength = BamTools::UnpackUnsignedInt(&x[16]);
    al->setMateRefID(BamTools::UnpackSignedInt(&x[20]));
    al->setMatePosition(BamTools::UnpackSignedInt(&x[24]));
    al->setInsertSize(BamTools::UnpackSignedInt(&x[28]));
    
    // read in character data - make sure proper data size was read
    bool readCharDataOK = false;
    const unsigned int dataLength = SupportData.BlockLength - 32;
    SupportData.AllCharData.resize(dataLength);
    input_stream.read(&SupportData.AllCharData[0], dataLength);
    if ( !input_stream.fail() ) {
        // set success flag
        readCharDataOK = true;
    } else {
        std::cerr << "Expected more bytes reading BAM char data. Is this file truncated or corrupted?" << std::endl;
        delete al;
        return NULL;
    }

    // return success/failure
    return al;
}


#endif