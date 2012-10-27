// ***************************************************************************
// BamAlignment.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 4 April 2012 (DB)
// ---------------------------------------------------------------------------
// Provides the BamAlignment data structure
// ***************************************************************************

#include "api/BamAlignment.h"
#include "api/BamConstants.h"

#include <sstream>

using namespace BamTools;
using namespace std;

/*! \class BamTools::BamAlignment
    \brief The main BAM alignment data structure.

    Provides methods to query/modify BAM alignment data fields.
*/
/*! \var BamAlignment::Name
    \brief read name
*/
/*! \var BamAlignment::Length
    \brief length of query sequence
*/
/*! \var BamAlignment::QueryBases
    \brief 'original' sequence (as reported from sequencing machine)

    \note Setting this field to "*" indicates that the sequence is not to be stored on output.
    In this case, the contents of the Qualities field should be invalidated as well (cleared or marked as "*").
*/
/*! \var BamAlignment::Qualities
    \brief FASTQ qualities (ASCII characters, not numeric values)

    \note Setting this field to "*" indicates to BamWriter that the quality scores are not to be stored,
    but instead will be output as a sequence of '0xFF'. Otherwise, QueryBases must not be a "*" and
    the length of this field should equal the length of QueryBases.
*/
/*! \var BamAlignment::TagData
    \brief tag data (use the provided methods to query/modify)
*/
/*! \var BamAlignment::RefID
    \brief ID number for reference sequence
*/
/*! \var BamAlignment::Position
    \brief position (0-based) where alignment starts
*/
/*! \var BamAlignment::Bin
    \brief BAM (standard) index bin number for this alignment
*/
/*! \var BamAlignment::MapQuality
    \brief mapping quality score
*/
/*! \var BamAlignment::AlignmentFlag
    \brief alignment bit-flag (use the provided methods to query/modify)
*/
/*! \var BamAlignment::CigarData
    \brief CIGAR operations for this alignment
*/
/*! \var BamAlignment::MateRefID
    \brief ID number for reference sequence where alignment's mate was aligned
*/
/*! \var BamAlignment::MatePosition
    \brief position (0-based) where alignment's mate starts
*/
/*! \var BamAlignment::InsertSize
    \brief mate-pair insert size
*/
/*! \var BamAlignment::Filename
    \brief name of BAM file which this alignment comes from
*/

/*! \fn BamAlignment::BamAlignment(void)
    \brief constructor
*/
BamAlignment::BamAlignment(void)
    : RefID(-1)
    , Position(-1)
    , Bin(0)
    , MapQuality(0)
    , AlignmentFlag(0)
    , MateRefID(-1)
    , MatePosition(-1)
    , InsertSize(0)
{ }

/*! \fn BamAlignment::BamAlignment(const BamAlignment& other)
    \brief copy constructor
*/
BamAlignment::BamAlignment(const BamAlignment& other)
    : Name(other.Name)
    , QueryBases(other.QueryBases)
    , Qualities(other.Qualities)
    , TagData(other.TagData)
    , RefID(other.RefID)
    , Position(other.Position)
    , Bin(other.Bin)
    , MapQuality(other.MapQuality)
    , AlignmentFlag(other.AlignmentFlag)
    , CigarData(other.CigarData)
    , MateRefID(other.MateRefID)
    , MatePosition(other.MatePosition)
    , InsertSize(other.InsertSize)
    , SupportData(other.SupportData)
{ }

/*! \fn BamAlignment::~BamAlignment(void)
    \brief destructor
*/
BamAlignment::~BamAlignment(void) { }

void BamAlignment::clear() {
    RefID = -1;
    Position = -1;
    Bin = 0;
    MapQuality = 0;
    AlignmentFlag = 0;
    MateRefID = -1;
    MatePosition = -1;
    InsertSize = 0;

    CigarData.clear();
    Name.clear();
    Qualities.clear();
    QueryBases.clear();
    TagData.clear();
    SupportData.clear();
}

void DecodeSequenceData(const string & encoded_sequence, string & decoded_sequence, size_t size) {
    
    decoded_sequence.resize(size);
    for ( size_t i = 0; i < size; ++i ) {
        const char singleBase = Constants::BAM_DNA_LOOKUP[ ( (encoded_sequence[(i/2)] >> (4*(1-(i%2)))) & 0xf ) ];
        decoded_sequence[i] = singleBase;
    }
}

// The implementation of this function is originally from BamWriter_p.cpp from bamtools.
// Bamtools is released under the BSD license.
void CreatePackedCigar(const std::vector<BamTools::CigarOp>& cigarOperations, std::string& packedCigar) {
    
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
void EncodeQuerySequence(const std::string& query, std::string& encodedQuery) {
    
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

bool BamAlignment::BuildName(void) const {
    lazy_load_lock.lock();
    
    if(!Name.empty()) {
        lazy_load_lock.unlock();
        return false;
    }
    
    Name = SupportData.getName();

    lazy_load_lock.unlock();
    
    return true;
}

bool BamAlignment::BuildQualitiesData(void) const {
    lazy_load_lock.lock();
    
    if(!Qualities.empty()) {
        lazy_load_lock.unlock();
        return false;
    }
    
    Qualities = SupportData.getQual();
    
    lazy_load_lock.unlock();
    
    return true;
}

bool BamAlignment::BuildQueryBasesData(void) const {
    lazy_load_lock.lock();
    
    if(!QueryBases.empty()) {
        lazy_load_lock.unlock();
        return false;
    }
    
    QueryBases = SupportData.getSeq();
    
    lazy_load_lock.unlock();
    
    return true;
}

bool BamAlignment::BuildTagData(void) const {
    lazy_load_lock.lock();
    
    if(!TagData.empty()) {
        lazy_load_lock.unlock();
        return false;
    }
    
    TagData = SupportData.getTagData();
    
    lazy_load_lock.unlock();
    
    return true;
}

bool BamAlignment::BuildCigarData(void) const {
    lazy_load_lock.lock();
    
    if(!CigarData.empty()) {
        lazy_load_lock.unlock();
        return false;
    }
    
    CigarData = SupportData.getCigar();
    
    lazy_load_lock.unlock();
    
    return true;
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

/*! \fn bool BamAlignment::FindTag(const std::string& tag, char*& pTagData, const unsigned int& tagDataLength, unsigned int& numBytesParsed) const
    \internal

    Searches for requested tag in BAM tag data.

    \param[in]     tag            requested 2-character tag name
    \param[in,out] pTagData       pointer to current position in BamAlignment::TagData
    \param[in]     tagDataLength  length of BamAlignment::TagData
    \param[in,out] numBytesParsed number of bytes parsed so far

    \return \c true if found

    \post If \a tag is found, \a pTagData will point to the byte where the tag data begins.
          \a numBytesParsed will correspond to the position in the full TagData string.

*/
bool BamAlignment::FindTag(const std::string& tag,
                           char*& pTagData,
                           const unsigned int& tagDataLength,
                           unsigned int& numBytesParsed) const
{
    
    BuildTagData();

    while ( numBytesParsed < tagDataLength ) {

        const char* pTagType        = pTagData;
        const char* pTagStorageType = pTagData + 2;
        pTagData       += 3;
        numBytesParsed += 3;

        // check the current tag, return true on match
        if ( strncmp(pTagType, tag.c_str(), 2) == 0 )
            return true;

        // get the storage class and find the next tag
        if ( *pTagStorageType == '\0' ) return false;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return false;
        if ( *pTagData == '\0' ) return false;
    }

    // checked all tags, none match
    return false;
}

/*! \fn int BamAlignment::GetEndPosition(bool usePadded = false, bool closedInterval = false) const
    \brief Calculates alignment end position, based on its starting position and CIGAR data.

    \warning The position returned now represents a zero-based, HALF-OPEN interval.
    In previous versions of BamTools (0.x & 1.x) all intervals were treated
    as zero-based, CLOSED.

    \param[in] usePadded      Allow inserted bases to affect the reported position. Default is
                              false, so that reported position stays synced with reference
                              coordinates.
    \param[in] closedInterval Setting this to true will return a 0-based end coordinate. Default is
                              false, so that his value represents a standard, half-open interval.

    \return alignment end position
*/
int BamAlignment::GetEndPosition(bool usePadded, bool closedInterval) const {

    // initialize alignment end to starting position
    int alignEnd = getPosition();

    // iterate over cigar operations
    vector<CigarOp>::const_iterator cigarIter = getCigarData().begin();
    vector<CigarOp>::const_iterator cigarEnd  = getCigarData().end();
    for ( ; cigarIter != cigarEnd; ++cigarIter) {
        const CigarOp& op = (*cigarIter);

        switch ( op.Type ) {

            // increase end position on CIGAR chars [DMXN=]
            case Constants::BAM_CIGAR_DEL_CHAR      :
            case Constants::BAM_CIGAR_MATCH_CHAR    :
            case Constants::BAM_CIGAR_MISMATCH_CHAR :
            case Constants::BAM_CIGAR_REFSKIP_CHAR  :
            case Constants::BAM_CIGAR_SEQMATCH_CHAR :
                alignEnd += op.Length;
                break;

            // increase end position on insertion, only if @usePadded is true
            case Constants::BAM_CIGAR_INS_CHAR :
                if ( usePadded )
                    alignEnd += op.Length;
                break;

            // all other CIGAR chars do not affect end position
            default :
                break;
        }
    }

    // adjust for closedInterval, if requested
    if ( closedInterval )
        alignEnd -= 1;

    // return result
    return alignEnd;
}

/*! \fn std::string BamAlignment::GetErrorString(void) const
    \brief Returns a human-readable description of the last error that occurred

    This method allows elimination of STDERR pollution. Developers of client code
    may choose how the messages are displayed to the user, if at all.

    \return error description
*/
std::string BamAlignment::GetErrorString(void) const {
    return ErrorString;
}

/*! \fn bool BamAlignment::GetTagType(const std::string& tag, char& type) const
    \brief Retrieves the BAM tag type-code associated with requested tag name.

    \param[in]  tag  2-character tag name
    \param[out] type retrieved (1-character) type-code

    \return \c true if found
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::GetTagType(const std::string& tag, char& type) const {
  
    BuildTagData();

    // skip if no tags present
    if ( TagData.empty() ) {
        // TODO: set error string?
        return false;
    }

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag not found, return failure
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ){
        // TODO: set error string?
        return false;
    }

    // otherwise, retrieve & validate tag type code
    type = *(pTagData - 1);
    switch ( type ) {
        case (Constants::BAM_TAG_TYPE_ASCII)  :
        case (Constants::BAM_TAG_TYPE_INT8)   :
        case (Constants::BAM_TAG_TYPE_UINT8)  :
        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
        case (Constants::BAM_TAG_TYPE_FLOAT)  :
        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
        case (Constants::BAM_TAG_TYPE_ARRAY)  :
            return true;

        // unknown tag type
        default:
            const string message = string("invalid tag type: ") + type;
            SetErrorString("BamAlignment::GetTagType", message);
            return false;
    }
}

/*! \fn bool BamAlignment::HasTag(const std::string& tag) const
    \brief Returns true if alignment has a record for requested tag.

    \param[in] tag 2-character tag name
    \return \c true if alignment has a record for tag
*/
bool BamAlignment::HasTag(const std::string& tag) const {

    BuildTagData();

    // return false if no tag data present
    if ( TagData.empty() )
        return false;

    // localize the tag data for lookup
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if result of tag lookup
    return FindTag(tag, pTagData, tagDataLength, numBytesParsed);
}

/*! \fn bool BamAlignment::IsDuplicate(void) const
    \return \c true if this read is a PCR duplicate
*/
bool BamAlignment::IsDuplicate(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_DUPLICATE) != 0 );
}

/*! \fn bool BamAlignment::IsFailedQC(void) const
    \return \c true if this read failed quality control
*/
bool BamAlignment::IsFailedQC(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_QC_FAILED) != 0 );
}

/*! \fn bool BamAlignment::IsFirstMate(void) const
    \return \c true if alignment is first mate on paired-end read
*/
bool BamAlignment::IsFirstMate(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_READ_1) != 0 );
}

/*! \fn bool BamAlignment::IsMapped(void) const
    \return \c true if alignment is mapped
*/
bool BamAlignment::IsMapped(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_UNMAPPED) == 0 );
}

/*! \fn bool BamAlignment::IsMateMapped(void) const
    \return \c true if alignment's mate is mapped
*/
bool BamAlignment::IsMateMapped(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_MATE_UNMAPPED) == 0 );
}

/*! \fn bool BamAlignment::IsMateReverseStrand(void) const
    \return \c true if alignment's mate mapped to reverse strand
*/
bool BamAlignment::IsMateReverseStrand(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_MATE_REVERSE_STRAND) != 0 );
}

/*! \fn bool BamAlignment::IsPaired(void) const
    \return \c true if alignment part of paired-end read
*/
bool BamAlignment::IsPaired(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_PAIRED) != 0 );
}

/*! \fn bool BamAlignment::IsPrimaryAlignment(void) const
    \return \c true if reported position is primary alignment
*/
bool BamAlignment::IsPrimaryAlignment(void) const  {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_SECONDARY) == 0 );
}

/*! \fn bool BamAlignment::IsProperPair(void) const
    \return \c true if alignment is part of read that satisfied paired-end resolution
*/
bool BamAlignment::IsProperPair(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_PROPER_PAIR) != 0 );
}

/*! \fn bool BamAlignment::IsReverseStrand(void) const
    \return \c true if alignment mapped to reverse strand
*/
bool BamAlignment::IsReverseStrand(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_REVERSE_STRAND) != 0 );
}

/*! \fn bool BamAlignment::IsSecondMate(void) const
    \return \c true if alignment is second mate on read
*/
bool BamAlignment::IsSecondMate(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_READ_2) != 0 );
}

/*! \fn bool BamAlignment::IsValidSize(const std::string& tag, const std::string& type) const
    \internal

    Checks that tag name & type strings are expected sizes.

    \param tag[in]  BAM tag name
    \param type[in] BAM tag type-code
    \return \c true if both input strings are valid sizes
*/
bool BamAlignment::IsValidSize(const std::string& tag, const std::string& type) const {
    return (tag.size()  == Constants::BAM_TAG_TAGSIZE) &&
           (type.size() == Constants::BAM_TAG_TYPESIZE);
}

/*! \fn void BamAlignment::RemoveTag(const std::string& tag)
    \brief Removes field from BAM tags.

    \param[in] tag 2-character name of field to remove
*/
void BamAlignment::RemoveTag(const std::string& tag) {
  
    BuildTagData();

    // skip if no tags available
    if ( TagData.empty() )
        return;
  
    // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;

    // skip if tag not found
    if  ( !FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) )
        return;

    // otherwise, remove it
    char *  newTagData = (char *) alloca(originalTagDataLength);

    // copy original tag data up til desired tag
    pTagData       -= 3;
    numBytesParsed -= 3;
    const unsigned int beginningTagDataLength = numBytesParsed;
    newTagDataLength += beginningTagDataLength;
    memcpy(newTagData, pOriginalTagData, numBytesParsed);

    // attemp to skip to next tag
    const char* pTagStorageType = pTagData + 2;
    pTagData       += 3;
    numBytesParsed += 3;
    if ( SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) {

        // squeeze remaining tag data
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagDataLength = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + beginningTagDataLength, pTagData, endTagDataLength );

        // save modified tag data in alignment
        TagData.assign(newTagData, beginningTagDataLength + endTagDataLength);
    }

    SupportData.setTagData(TagData);
}

/*! \fn void BamAlignment::SetErrorString(const std::string& where, const std::string& what) const
    \internal

    Sets a formatted error string for this alignment.

    \param[in] where class/method where error occurred
    \param[in] what  description of error
*/
void BamAlignment::SetErrorString(const std::string& where, const std::string& what) const {
    static const string SEPARATOR = ": ";
    ErrorString = where + SEPARATOR + what;
}

/*! \fn void BamAlignment::SetIsDuplicate(bool ok)
    \brief Sets value of "PCR duplicate" flag to \a ok.
*/
void BamAlignment::SetIsDuplicate(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_DUPLICATE;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_DUPLICATE;
}

/*! \fn void BamAlignment::SetIsFailedQC(bool ok)
    \brief Sets "failed quality control" flag to \a ok.
*/
void BamAlignment::SetIsFailedQC(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_QC_FAILED;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_QC_FAILED;
}

/*! \fn void BamAlignment::SetIsFirstMate(bool ok)
    \brief Sets "alignment is first mate" flag to \a ok.
*/
void BamAlignment::SetIsFirstMate(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_READ_1;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_READ_1;
}

/*! \fn void BamAlignment::SetIsMapped(bool ok)
    \brief Sets "alignment is mapped" flag to \a ok.
*/
void BamAlignment::SetIsMapped(bool ok) {
    if (ok) AlignmentFlag &= ~Constants::BAM_ALIGNMENT_UNMAPPED;
    else    AlignmentFlag |=  Constants::BAM_ALIGNMENT_UNMAPPED;
}

/*! \fn void BamAlignment::SetIsMateMapped(bool ok)
    \brief Sets "alignment's mate is mapped" flag to \a ok.
*/
void BamAlignment::SetIsMateMapped(bool ok) {
    if (ok) AlignmentFlag &= ~Constants::BAM_ALIGNMENT_MATE_UNMAPPED;
    else    AlignmentFlag |=  Constants::BAM_ALIGNMENT_MATE_UNMAPPED;
}

/*! \fn void BamAlignment::SetIsMateReverseStrand(bool ok)
    \brief Sets "alignment's mate mapped to reverse strand" flag to \a ok.
*/
void BamAlignment::SetIsMateReverseStrand(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_MATE_REVERSE_STRAND;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_MATE_REVERSE_STRAND;
}

/*! \fn void BamAlignment::SetIsPaired(bool ok)
    \brief Sets "alignment part of paired-end read" flag to \a ok.
*/
void BamAlignment::SetIsPaired(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_PAIRED;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_PAIRED;
}

/*! \fn void BamAlignment::SetIsPrimaryAlignment(bool ok)
    \brief Sets "position is primary alignment" flag to \a ok.
*/
void BamAlignment::SetIsPrimaryAlignment(bool ok) {
    if (ok) AlignmentFlag &= ~Constants::BAM_ALIGNMENT_SECONDARY;
    else    AlignmentFlag |=  Constants::BAM_ALIGNMENT_SECONDARY;
}

/*! \fn void BamAlignment::SetIsProperPair(bool ok)
    \brief Sets "alignment is part of read that satisfied paired-end resolution" flag to \a ok.
*/
void BamAlignment::SetIsProperPair(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_PROPER_PAIR;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_PROPER_PAIR;
}

/*! \fn void BamAlignment::SetIsReverseStrand(bool ok)
    \brief Sets "alignment mapped to reverse strand" flag to \a ok.
*/
void BamAlignment::SetIsReverseStrand(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_REVERSE_STRAND;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_REVERSE_STRAND;
}

/*! \fn void BamAlignment::SetIsSecondMate(bool ok)
    \brief Sets "alignment is second mate on read" flag to \a ok.
*/
void BamAlignment::SetIsSecondMate(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_READ_2;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_READ_2;
}

/*! \fn bool BamAlignment::SkipToNextTag(const char storageType, char*& pTagData, unsigned int& numBytesParsed) const
    \internal

    Moves to next available tag in tag data string

    \param[in]     storageType    BAM tag type-code that determines how far to move cursor
    \param[in,out] pTagData       pointer to current position (cursor) in tag string
    \param[in,out] numBytesParsed report of how many bytes were parsed (cumulatively)

    \return \c if storageType was a recognized BAM tag type

    \post \a pTagData       will point to the byte where the next tag data begins.
          \a numBytesParsed will correspond to the cursor's position in the full TagData string.
*/
bool BamAlignment::SkipToNextTag(const char storageType,
                                 char*& pTagData,
                                 unsigned int& numBytesParsed) const
{
    switch (storageType) {

        case (Constants::BAM_TAG_TYPE_ASCII) :
        case (Constants::BAM_TAG_TYPE_INT8)  :
        case (Constants::BAM_TAG_TYPE_UINT8) :
            ++numBytesParsed;
            ++pTagData;
            break;

        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
            numBytesParsed += sizeof(uint16_t);
            pTagData       += sizeof(uint16_t);
            break;

        case (Constants::BAM_TAG_TYPE_FLOAT)  :
        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
            numBytesParsed += sizeof(uint32_t);
            pTagData       += sizeof(uint32_t);
            break;

        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
            while( *pTagData ) {
                ++numBytesParsed;
                ++pTagData;
            }
            // increment for null-terminator
            ++numBytesParsed;
            ++pTagData;
            break;

        case (Constants::BAM_TAG_TYPE_ARRAY) :

        {
            // read array type
            const char arrayType = *pTagData;
            ++numBytesParsed;
            ++pTagData;

            // read number of elements
            int32_t numElements;
            memcpy(&numElements, pTagData, sizeof(uint32_t)); // already endian-swapped, if needed
            numBytesParsed += sizeof(uint32_t);
            pTagData       += sizeof(uint32_t);

            // calculate number of bytes to skip
            int bytesToSkip = 0;
            switch (arrayType) {
                case (Constants::BAM_TAG_TYPE_INT8)  :
                case (Constants::BAM_TAG_TYPE_UINT8) :
                    bytesToSkip = numElements;
                    break;
                case (Constants::BAM_TAG_TYPE_INT16)  :
                case (Constants::BAM_TAG_TYPE_UINT16) :
                    bytesToSkip = numElements*sizeof(uint16_t);
                    break;
                case (Constants::BAM_TAG_TYPE_FLOAT)  :
                case (Constants::BAM_TAG_TYPE_INT32)  :
                case (Constants::BAM_TAG_TYPE_UINT32) :
                    bytesToSkip = numElements*sizeof(uint32_t);
                    break;
                default:
                    const string message = string("invalid binary array type: ") + arrayType;
                    SetErrorString("BamAlignment::SkipToNextTag", message);
                    return false;
            }

            // skip binary array contents
            numBytesParsed += bytesToSkip;
            pTagData       += bytesToSkip;
            break;
        }

        default:
            const string message = string("invalid tag type: ") + storageType;
            SetErrorString("BamAlignment::SkipToNextTag", message);
            return false;
    }

    // if we get here, tag skipped OK - return success
    return true;
}

string cigarToString(const vector<CigarOp> cigar)
{
    stringstream ss;
    for(vector<CigarOp>::const_iterator i = cigar.begin(); i != cigar.end(); i++)
        ss << i->Length << i->Type;

    return string(ss.str());
}

//#pragma mark BamAlignmentSupportData
void BamAlignment::BamAlignmentSupportData::setCigar(const std::vector<CigarOp> & cigar) {
    string encoded_cigar;
    CreatePackedCigar(cigar, encoded_cigar);
    AllCharData.replace(beginCigar(), endCigar(), encoded_cigar);
    NumCigarOperations = cigar.size();
}

const std::vector<CigarOp> BamAlignment::BamAlignmentSupportData::getCigar() const {
    string encoded_cigar(beginCigar(), endCigar());
    uint32_t * cigarData = (uint32_t *) &encoded_cigar[0];
    vector<CigarOp> CigarData;

    CigarData.reserve(NumCigarOperations);

    for ( unsigned int i = 0; i < NumCigarOperations; ++i ) {
        
        // build CigarOp structure
        CigarOp op;
        op.Length = (cigarData[i] >> Constants::BAM_CIGAR_SHIFT);
        op.Type   = Constants::BAM_CIGAR_LOOKUP[ (cigarData[i] & Constants::BAM_CIGAR_MASK) ];
        
        // save CigarOp
        CigarData.push_back(op);
    }
    return CigarData;
}

void BamAlignment::BamAlignmentSupportData::setSeq(const std::string & seq) {
    string encoded_seq;
    EncodeQuerySequence(seq, encoded_seq);
    AllCharData.replace(beginSeq(), endSeq(), encoded_seq);
    QuerySequenceLength = seq.size();
}

const std::string BamAlignment::BamAlignmentSupportData::getSeq() const {
    string decoded;

    DecodeSequenceData(string(beginSeq(), endSeq()), decoded, QuerySequenceLength);

    return decoded;
}

void BamAlignment::BamAlignmentSupportData::setQual(const std::string & seq) {
    AllCharData.replace(beginQual(), endQual(), seq);
    for(string::iterator i = beginQual(); i != endQual(); i++)
        *i -= 33;
}

const std::string BamAlignment::BamAlignmentSupportData::getQual() const {
    string ret(beginQual(), endQual());
    for(uint32_t i = 0; i < QuerySequenceLength; i++)
        ret[i] += 33;
    return ret;
}

