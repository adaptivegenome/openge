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
#include "api/internal/utils/BamThreadPool.h"

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
/*! \var BamAlignment::AlignedBases
    \brief 'aligned' sequence (includes any indels, padding, clipping)

    This field will be completely empty after reading from BamReader/BamMultiReader when
    QueryBases is empty.
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
    : hasAlignedBasesData(false)
    , hasQualitiesData(false)
    , hasQueryBasesData(false)
    , hasTagData(false)
    , hasCigarData(false)
    , stringDataDirty(false)
    , RefID(-1)
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
    : hasAlignedBasesData(false)
    , hasQualitiesData(false)
    , hasQueryBasesData(false)
    , hasTagData(false)
    , hasCigarData(false)
    , stringDataDirty(false)
    , Name(other.Name)
    , QueryBases(other.QueryBases)
    , AlignedBases(other.AlignedBases)
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
    hasAlignedBasesData = false;
    hasQualitiesData = false;
    hasQueryBasesData = false;
    hasTagData = false;
    hasCigarData = false;
    stringDataDirty = false;
    RefID = -1;
    Position = -1;
    Bin = 0;
    MapQuality = 0;
    AlignmentFlag = 0;
    MateRefID = -1;
    MatePosition = -1;
    InsertSize = 0;
    
    AlignedBases.clear();
    CigarData.clear();
    Name.clear();
    Qualities.clear();
    QueryBases.clear();
    TagData.clear();
}

bool BamAlignment::BuildQueryBasesData() const {
    lazy_load_lock.lock();
    if(hasQueryBasesData) {
        lazy_load_lock.unlock();
        return true;
    }

    hasQueryBasesData = true;

    // calculate character lengths/offsets
    const unsigned int seqDataOffset  = SupportData.QueryNameLength + (SupportData.NumCigarOperations*4);
    const unsigned int qualDataOffset = seqDataOffset + (SupportData.QuerySequenceLength+1)/2;
    
    // check offsets to see what char data exists
    const bool hasSeqData  = ( seqDataOffset  < qualDataOffset );
    // save query sequence
    QueryBases.clear();
    if ( hasSeqData ) {
        const char* seqData = SupportData.AllCharData.data() + seqDataOffset;
        QueryBases.reserve(SupportData.QuerySequenceLength);
        for ( size_t i = 0; i < SupportData.QuerySequenceLength; ++i ) {
            const char singleBase = Constants::BAM_DNA_LOOKUP[ ( (seqData[(i/2)] >> (4*(1-(i%2)))) & 0xf ) ];
            QueryBases.append(1, singleBase);
        }
    }

    lazy_load_lock.unlock();

    return true;
}

bool BamAlignment::BuildQualitiesData() const {
    lazy_load_lock.lock();
    if(hasQualitiesData) {
        lazy_load_lock.unlock();
        return true;
    }

    hasQualitiesData = true;

    // calculate character lengths/offsets
    const unsigned int seqDataOffset  = SupportData.QueryNameLength + (SupportData.NumCigarOperations*4);
    const unsigned int qualDataOffset = seqDataOffset + (SupportData.QuerySequenceLength+1)/2;
    const unsigned int tagDataOffset  = qualDataOffset + SupportData.QuerySequenceLength;
    
    // check offsets to see what char data exists
    const bool hasQualData = ( qualDataOffset < tagDataOffset );

    Qualities.clear();
    if ( hasQualData ) {
        const char* qualData = SupportData.AllCharData.data() + qualDataOffset;

        // if marked as unstored (sequence of 0xFF) - don't do conversion, just fill with 0xFFs
        if ( qualData[0] == (char)0xFF )
            Qualities.resize(SupportData.QuerySequenceLength, (char)0xFF);

        // otherwise convert from numeric QV to 'FASTQ-style' ASCII character
        else {
            Qualities.reserve(SupportData.QuerySequenceLength);
            for ( size_t i = 0; i < SupportData.QuerySequenceLength; ++i )
                Qualities.append(1, qualData[i]+33);
        }
    }
    
    lazy_load_lock.unlock();
    
    return true;
}

bool BamAlignment::BuildAlignedBasesData() const {
    lazy_load_lock.lock();
    if(hasAlignedBasesData) {
        lazy_load_lock.unlock();
        return true;
    }

    hasAlignedBasesData = true;

    // clear previous AlignedBases
    AlignedBases.clear();

    // if QueryBases has data, build AlignedBases using CIGAR data
    // otherwise, AlignedBases will remain empty (this case IS allowed)
    if ( !QueryBases.empty() && QueryBases != "*" ) {

        // resize AlignedBases
        AlignedBases.reserve(SupportData.QuerySequenceLength);

        // iterate over CigarOps
        int k = 0;
        vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
        vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
        for ( ; cigarIter != cigarEnd; ++cigarIter ) {
            const CigarOp& op = (*cigarIter);

            switch ( op.Type ) {

                // for 'M', 'I', '=', 'X' - write bases
                case (Constants::BAM_CIGAR_MATCH_CHAR)    :
                case (Constants::BAM_CIGAR_INS_CHAR)      :
                case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
                case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
                    AlignedBases.append(QueryBases.substr(k, op.Length));
                    // fall through

                // for 'S' - soft clip, do not write bases
                // but increment placeholder 'k'
                case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
                    k += op.Length;
                    break;

                // for 'D' - write gap character
                case (Constants::BAM_CIGAR_DEL_CHAR) :
                    AlignedBases.append(op.Length, Constants::BAM_DNA_DEL);
                    break;

                // for 'P' - write padding character
                case (Constants::BAM_CIGAR_PAD_CHAR) :
                    AlignedBases.append( op.Length, Constants::BAM_DNA_PAD );
                    break;

                // for 'N' - write N's, skip bases in original query sequence
                case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
                    AlignedBases.append( op.Length, Constants::BAM_DNA_N );
                    break;

                // for 'H' - hard clip, do nothing to AlignedBases, move to next op
                case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
                    break;

                // invalid CIGAR op-code
                default:
                    const string message = string("invalid CIGAR operation type: ") + op.Type;
                    SetErrorString("BamAlignment::BuildAlignedBasesData", message);
                    return false;
            }
        }
    }

    lazy_load_lock.unlock();

    return true;
}

bool BamAlignment::BuildTagData() const {
    lazy_load_lock.lock();
    if(hasTagData) {
        lazy_load_lock.unlock();
        return true;
    }

    hasTagData = true;

    // calculate character lengths/offsets
    const unsigned int dataLength     = SupportData.BlockLength - Constants::BAM_CORE_SIZE;
    const unsigned int seqDataOffset  = SupportData.QueryNameLength + (SupportData.NumCigarOperations*4);
    const unsigned int qualDataOffset = seqDataOffset + (SupportData.QuerySequenceLength+1)/2;
    const unsigned int tagDataOffset  = qualDataOffset + SupportData.QuerySequenceLength;
    const unsigned int tagDataLength  = dataLength - tagDataOffset;
    
    // check offsets to see what char data exists
    const bool hasTagData  = ( tagDataOffset  < dataLength );
    TagData.clear();
    if ( hasTagData ) {

        char* tagData = (((char*)SupportData.AllCharData.data()) + tagDataOffset);

        // store tagData in alignment
        TagData.resize(tagDataLength);
        memcpy((char*)(TagData.data()), tagData, tagDataLength);
    }
    
    lazy_load_lock.unlock();

    return true;
}

bool BamAlignment::BuildCigarData() const {
    
    lazy_load_lock.lock();

    if(hasCigarData) {
        lazy_load_lock.unlock();
        return true;
    }

    hasCigarData = true;

    // save CIGAR ops
    // need to calculate this here so that  BamAlignment::GetEndPosition() performs correctly,
    // even when GetNextAlignment() is called
    CigarData.reserve(SupportData.NumCigarOperations);
    const unsigned int cigarDataOffset = SupportData.QueryNameLength;
    uint32_t* cigarData = (uint32_t*)(&SupportData.AllCharData[0] + cigarDataOffset);
    for ( unsigned int i = 0; i < SupportData.NumCigarOperations; ++i ) {

        // build CigarOp structure
        CigarOp op;
        op.Length = (cigarData[i] >> Constants::BAM_CIGAR_SHIFT);
        op.Type   = Constants::BAM_CIGAR_LOOKUP[ (cigarData[i] & Constants::BAM_CIGAR_MASK) ];
        
        // save CigarOp
        CigarData.push_back(op);
    }
    
    lazy_load_lock.unlock();
    
    return true;
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

void BamAlignment::FlushCharData(void) const {
    
    BuildTagData();
    BuildQueryBasesData();
    BuildCigarData();
    BuildQualitiesData();

    lazy_load_lock.lock();
    if(!stringDataDirty) {
        lazy_load_lock.unlock();
        return;
    }

    stringDataDirty = false;

    // calculate char lengths
    const unsigned int nameLength         = Name.size() + 1;
    const unsigned int numCigarOperations = CigarData.size();
    const unsigned int queryLength        = QueryBases.size();
    const unsigned int tagDataLength      = TagData.size();
    
    // no way to tell if alignment's bin is already defined (there is no default, invalid value)
    // so we'll go ahead calculate its bin ID before storing
    const uint32_t alignmentBin = CalculateMinimumBin(getPosition(), GetEndPosition());
    
    // create our packed cigar string
    std::string packedCigar;
    CreatePackedCigar(CigarData, packedCigar);
    const unsigned int packedCigarLength = packedCigar.size();
    
    stringstream output_stream;
    
    // encode the query
    std::string encodedQuery;
    EncodeQuerySequence(QueryBases, encodedQuery);
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
    buffer[0] = getRefID();
    buffer[1] = getPosition();
    buffer[2] = (alignmentBin << 16) | (getMapQuality() << 8) | nameLength;
    buffer[3] = (getAlignmentFlag() << 16) | numCigarOperations;
    buffer[4] = queryLength;
    buffer[5] = getMateRefID();
    buffer[6] = getMatePosition();
    buffer[7] = getInsertSize();
    
    // write the BAM core
    output_stream.write((char*)&buffer, 32);
    
    // write the query name
    output_stream.write(Name.c_str(), nameLength);
    
    output_stream.write(packedCigar.data(), packedCigarLength);
    
    // write the encoded query sequence
    output_stream.write(encodedQuery.data(), encodedQueryLength);
    
    // write the base qualities
    char* pBaseQualities = (char*)Qualities.data();
    char * qualities = (char *) alloca(Qualities.size());
    for ( size_t i = 0; i < queryLength; ++i )
        qualities[i] = pBaseQualities[i] - 33; // FASTQ conversion
    output_stream.write(qualities, queryLength);
    
    output_stream.write(TagData.data(), tagDataLength);
    
    SupportData.AllCharData = output_stream.str();
    lazy_load_lock.unlock();
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
    int alignEnd = Position;

    // iterate over cigar operations
    vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
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

/*! \fn bool BamAlignment::GetSoftClips(std::vector<int>& clipSizes, std::vector<int>& readPositions, std::vector<int>& genomePositions, bool usePadded = false) const
    \brief Identifies if an alignment has a soft clip. If so, identifies the
           sizes of the soft clips, as well as their positions in the read and reference.

    \param[out] clipSizes       vector of the sizes of each soft clip in the alignment
    \param[out] readPositions   vector of the 0-based read locations of each soft clip in the alignment.
                                These positions are basically indexes within the read, not genomic positions.
    \param[out] genomePositions vector of the 0-based genome locations of each soft clip in the alignment
    \param[in]  usePadded       inserted bases affect reported position. Default is false, so that
                                reported position stays 'sync-ed' with reference coordinates.

    \return \c true if any soft clips were found in the alignment
*/
bool BamAlignment::GetSoftClips(vector<int>& clipSizes,
                                vector<int>& readPositions,
                                vector<int>& genomePositions,
                                bool usePadded) const
{
    // initialize positions & flags
    int refPosition  = Position;
    int readPosition = 0;
    bool softClipFound = false;
    bool firstCigarOp  = true;

    // iterate over cigar operations
    vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter) {
        const CigarOp& op = (*cigarIter);

        switch ( op.Type ) {

            // increase both read & genome positions on CIGAR chars [DMXN=]
            case Constants::BAM_CIGAR_DEL_CHAR      :
            case Constants::BAM_CIGAR_MATCH_CHAR    :
            case Constants::BAM_CIGAR_MISMATCH_CHAR :
            case Constants::BAM_CIGAR_REFSKIP_CHAR  :
            case Constants::BAM_CIGAR_SEQMATCH_CHAR :
                refPosition  += op.Length;
                readPosition += op.Length;
                break;

            // increase read position on insertion, genome position only if @usePadded is true
            case Constants::BAM_CIGAR_INS_CHAR :
                readPosition += op.Length;
                if ( usePadded )
                    refPosition += op.Length;
                break;

            case Constants::BAM_CIGAR_SOFTCLIP_CHAR :

                softClipFound = true;

                //////////////////////////////////////////////////////////////////////////////
                // if we are dealing with the *first* CIGAR operation
                // for this alignment, we increment the read position so that
                // the read and genome position of the clip are referring to the same base.
                // For example, in the alignment below, the ref position would be 4, yet
                //              the read position would be 0. Thus, to "sync" the two,
                //              we need to increment the read position by the length of the
                //              soft clip.
                // Read:  ATCGTTTCGTCCCTGC
                // Ref:   GGGATTTCGTCCCTGC
                // Cigar: SSSSMMMMMMMMMMMM
                //
                // NOTE: This only needs to be done if the soft clip is the _first_ CIGAR op.
                //////////////////////////////////////////////////////////////////////////////
                if ( firstCigarOp )
                    readPosition += op.Length;

                // track the soft clip's size, read position, and genome position
                clipSizes.push_back(op.Length);
                readPositions.push_back(readPosition);
                genomePositions.push_back(refPosition);

            // any other CIGAR operations have no effect
            default :
                break;
        }

        // clear our "first pass" flag
        firstCigarOp = false;
    }

    // return whether any soft clips found
    return softClipFound;
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
    RaiiBuffer newTagData(originalTagDataLength);

    // copy original tag data up til desired tag
    pTagData       -= 3;
    numBytesParsed -= 3;
    const unsigned int beginningTagDataLength = numBytesParsed;
    newTagDataLength += beginningTagDataLength;
    memcpy(newTagData.Buffer, pOriginalTagData, numBytesParsed);

    // attemp to skip to next tag
    const char* pTagStorageType = pTagData + 2;
    pTagData       += 3;
    numBytesParsed += 3;
    if ( SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) {

        // squeeze remaining tag data
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagDataLength = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData.Buffer + beginningTagDataLength, pTagData, endTagDataLength );

        // save modified tag data in alignment
        TagData.assign(newTagData.Buffer, beginningTagDataLength + endTagDataLength);
    }
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

//////////////
// cache allocations for performance

// this class should be a subclass of BamAlignment, but the way the API exports are structured,
// we would have to export more classes with the API. Instead, we'll just do this temporarily
// until bamtools is replaced.
class BamAlignClearJob : public BamThreadJob {
public:
    vector<BamAlignment *> * cached_allocations, *cached_allocations_cleared;
    BamSpinlock * allocator_spinlock;
    bool * clean_thread_running;
    virtual void runJob();
};

void BamAlignClearJob::runJob() {
    vector<BamAlignment *> to_clear;    //thread private copy of to-be-cleared alignments
    to_clear.reserve(100);
    static int ct = 0;
    while (true) {
        allocator_spinlock->lock();
        while(!cached_allocations->empty() && to_clear.size() < 100) {
            BamAlignment * al = cached_allocations->back();
            cached_allocations->pop_back();
            ct++;
            to_clear.push_back(al);
        }
        
        allocator_spinlock->unlock();
        
        if(to_clear.empty()) break;

        for(vector<BamAlignment *>::iterator i = to_clear.begin(); i!= to_clear.end(); i++)
            (*i)->clear();
        
        allocator_spinlock->lock();
        for(vector<BamAlignment *>::iterator i = to_clear.begin(); i!= to_clear.end(); i++)
            cached_allocations_cleared->push_back(*i);
        allocator_spinlock->unlock();
        
        to_clear.clear();
    }
    *clean_thread_running = false;
}

BamAlignment * BamAlignment::allocate() {
    BamAlignment * ret = NULL;
    allocator_spinlock.lock();

    if(!cached_allocations_cleared.empty()) {
        ret = cached_allocations_cleared.back();
        cached_allocations_cleared.pop_back();
    }

    allocator_spinlock.unlock();
    if(!ret) ret = new BamAlignment();
    //else ret->clear();

    return ret;
}

BamThreadPool allocator_pool;

void BamAlignment::deallocate(BamAlignment * al) {
    allocator_spinlock.lock();
    cached_allocations.push_back(al);
    allocator_spinlock.unlock();
    
    if(cached_allocations.size() > 100 && !clean_thread_running) {
        clean_thread_running = true;
        BamAlignClearJob * job = new BamAlignClearJob;
        job->allocator_spinlock = &allocator_spinlock;
        job->cached_allocations = & cached_allocations;
        job->clean_thread_running = &clean_thread_running;
        job->cached_allocations_cleared = &cached_allocations_cleared;
        allocator_pool.addJob(job);
    }
}

void BamAlignment::clearCachedAllocations() {
    allocator_pool.waitForJobCompletion();
    allocator_spinlock.lock();
    while(!cached_allocations.empty()) {
        BamAlignment * al = cached_allocations.back();
        delete(al);
        cached_allocations.pop_back();
    }
    allocator_spinlock.unlock();
}

BamSpinlock BamAlignment::allocator_spinlock;
std::vector<BamAlignment *> BamAlignment::cached_allocations;
std::vector<BamAlignment *> BamAlignment::cached_allocations_cleared;
bool BamAlignment::clean_thread_running = false;
