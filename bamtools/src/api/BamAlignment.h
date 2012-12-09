// ***************************************************************************
// BamAlignment.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 16 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the BamAlignment data structure
// ***************************************************************************

#ifndef BAMALIGNMENT_H
#define BAMALIGNMENT_H

#include "api/api_global.h"
#include "api/BamAux.h"
#include "api/BamConstants.h"
#include "api/BamParallelismSettings.h"
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>

template <class> class BamDeserializer;

std::string cigarToString(const std::vector<BamTools::CigarOp> cigar);

namespace BamTools {
    
    // BamAlignment data structure
    class API_EXPORT BamAlignment {
        
        // constructors & destructor
    public:
        BamAlignment(void);
        BamAlignment(const BamAlignment& other);
        ~BamAlignment(void);
        void clear();    //go back to a freshly constructed state
        
        // queries against alignment flags
    public:
        bool IsDuplicate(void) const;         // returns true if this read is a PCR duplicate
        bool IsFailedQC(void) const;          // returns true if this read failed quality control
        bool IsFirstMate(void) const;         // returns true if alignment is first mate on read
        bool IsMapped(void) const;            // returns true if alignment is mapped
        bool IsMateMapped(void) const;        // returns true if alignment's mate is mapped
        bool IsMateReverseStrand(void) const; // returns true if alignment's mate mapped to reverse strand
        bool IsPaired(void) const;            // returns true if alignment part of paired-end read
        bool IsPrimaryAlignment(void) const;  // returns true if reported position is primary alignment
        bool IsProperPair(void) const;        // returns true if alignment is part of read that satisfied paired-end resolution
        bool IsReverseStrand(void) const;     // returns true if alignment mapped to reverse strand
        bool IsSecondMate(void) const;        // returns true if alignment is second mate on read
        
        // manipulate alignment flags
    public:
        void SetIsDuplicate(bool ok);         // sets value of "PCR duplicate" flag
        void SetIsFailedQC(bool ok);          // sets value of "failed quality control" flag
        void SetIsFirstMate(bool ok);         // sets value of "alignment is first mate" flag
        void SetIsMapped(bool ok);            // sets value of "alignment is mapped" flag
        void SetIsMateMapped(bool ok);        // sets value of "alignment's mate is mapped" flag
        void SetIsMateReverseStrand(bool ok); // sets value of "alignment's mate mapped to reverse strand" flag
        void SetIsPaired(bool ok);            // sets value of "alignment part of paired-end read" flag
        void SetIsPrimaryAlignment(bool ok);  // sets value of "position is primary alignment" flag
        void SetIsProperPair(bool ok);        // sets value of "alignment is part of read that satisfied paired-end resolution" flag
        void SetIsReverseStrand(bool ok);     // sets value of "alignment mapped to reverse strand" flag
        void SetIsSecondMate(bool ok);        // sets value of "alignment is second mate on read" flag
        
        // tag data access methods
    public:
        
        // add a new tag
        template<typename T> bool AddTag(const std::string& tag, const std::string& type, const T& value);
        template<typename T> bool AddTag(const std::string& tag, const std::vector<T>& values);
        
        // edit (or append) tag
        template<typename T> bool EditTag(const std::string& tag, const std::string& type, const T& value);
        template<typename T> bool EditTag(const std::string& tag, const std::vector<T>& values);
        
        // retrieves tag data
        template<typename T> bool GetTag(const std::string& tag, T& destination) const;
        template<typename T> bool GetTag(const std::string& tag, std::vector<T>& destination) const;
        
        // retrieves the SAM/BAM type-code for requested tag name
        bool GetTagType(const std::string& tag, char& type) const;
        
        // returns true if alignment has a record for this tag name
        bool HasTag(const std::string& tag) const;
        
        // removes a tag
        void RemoveTag(const std::string& tag);
        
    public:
        // calculates alignment end position
        int GetEndPosition(bool usePadded = false, bool closedInterval = false) const;
        
        // returns a description of the last error that occurred
        std::string GetErrorString(void) const;

        // public data fields
    protected:
        int32_t     RefID;              // ID number for reference sequence
        int32_t     Position;           // position (0-based) where alignment starts
        uint16_t    Bin;                // BAM (standard) index bin number for this alignment
        uint16_t    MapQuality;         // mapping quality score
        uint32_t    AlignmentFlag;      // alignment bit-flag (use provided methods to query/modify)
        
        int32_t     MateRefID;          // ID number for reference sequence where alignment's mate was aligned
        int32_t     MatePosition;       // position (0-based) where alignment's mate starts
        int32_t     InsertSize;         // mate-pair insert size
        
    public:
        const std::string getName() const { return SupportData.getName(); }
        uint32_t getNameLength() const { return SupportData.getQueryNameLength(); }
        int32_t getLength() const { return SupportData.getQuerySequenceLength(); }
        const std::string getQueryBases() const { return SupportData.getSeq(); }
        uint32_t getQueryBasesLength() const { return SupportData.getQuerySequenceLength(); }
        const std::string getQualities() const { return SupportData.getQual(); }
        const std::string getTagData() const { return SupportData.getTagData(); }
        int32_t getRefID() const { return RefID; }
        int32_t getPosition() const { return Position; }
        uint16_t getBin() const { return Bin; }
        uint16_t getMapQuality() const { return MapQuality; }
        uint32_t getAlignmentFlag() const { return AlignmentFlag; }
        const std::vector<CigarOp> getCigarData() const { return SupportData.getCigar(); }
        uint32_t getNumCigarOps() const { return SupportData.getNumCigarOperations(); }
        int32_t getMateRefID() const { return MateRefID; }
        int32_t getMatePosition() const { return MatePosition; }
        int32_t getInsertSize() const { return InsertSize; }
        
        void setName(const std::string & newName) { SupportData.setName(newName); };
        void setQueryBases(const std::string & newQueryBases) { SupportData.setSeq(newQueryBases); };
        void setQualities(const std::string & newQualities) { SupportData.setQual(newQualities); };
        void setTagData(const std::string & newTagData) { SupportData.setTagData(newTagData); };
        void setRefID(int32_t newRefID) { RefID = newRefID; }
        void setPosition(int32_t newPosition) { Position = newPosition; }
        void setBin(uint16_t newBin) { Bin = newBin; }
        void setMapQuality(uint16_t newMapQuality) { MapQuality = newMapQuality; }
        void setAlignmentFlag(uint32_t newAlignmentFlag) { AlignmentFlag = newAlignmentFlag; }
        void setCigarData(const std::vector<CigarOp> & newCigarData) { SupportData.setCigar(newCigarData); }
        void setMateRefID(int32_t newMateRefID) { MateRefID = newMateRefID; }
        void setMatePosition(int32_t newMatePosition) { MatePosition = newMatePosition; }
        void setInsertSize(int32_t newInsertSize) { InsertSize = newInsertSize; }
        
        
        //! \internal
        // internal utility methods
    private:
        bool FindTag(const std::string& tag,
                     char*& pTagData,
                     const unsigned int& tagDataLength,
                     unsigned int& numBytesParsed) const;
        bool IsValidSize(const std::string& tag, const std::string& type) const;
        void SetErrorString(const std::string& where, const std::string& what) const;
        bool SkipToNextTag(const char storageType,
                           char*& pTagData,
                           unsigned int& numBytesParsed) const;
        
        // internal data
    public:
        
        class TagDataView {
            const char * p_data;
            size_t length;
        public:
            TagDataView(const char * start, size_t length) : p_data(start), length(length) {}
            const char & operator[](int i) const { assert(i < length && i >= 0); return p_data[i]; }
            size_t size() const { return length; }
            bool empty() const { return length == 0; }
            const char * data() const {return p_data; }
        };
        
        class BamAlignmentSupportData {
        protected:
            // data members
            std::string AllCharData;
            uint32_t    NumCigarOperations;
            uint32_t    QueryNameLength;
            uint32_t    QuerySequenceLength;
        public:
            // constructor
            BamAlignmentSupportData(void)
            : NumCigarOperations(0)
            , QueryNameLength(0)
            , QuerySequenceLength(0)
            {
                //put a null character in the string. This way, we recognize
                AllCharData = " ";
                AllCharData[0] = 0;
            }
            uint32_t getBlockLength() const { return AllCharData.size() + 32; }
            void clear() { AllCharData.clear(); NumCigarOperations = 0; QueryNameLength = 0; QuerySequenceLength = 0; }
            
        protected:
            // convenience
            std::string::iterator beginName() { return AllCharData.begin(); }
            std::string::iterator endName() { return AllCharData.begin() + QueryNameLength; }
            std::string::iterator beginCigar() { return endName(); }
            std::string::iterator endCigar() { return beginCigar() + (NumCigarOperations * 4); }
            std::string::iterator beginSeq() { return endCigar(); }
            std::string::iterator endSeq() { return beginSeq() + ((QuerySequenceLength + 1)/2); }
            std::string::iterator beginQual() { return endSeq(); }
            std::string::iterator endQual() { return beginQual() + QuerySequenceLength; }
            std::string::iterator beginTagData() { return endQual(); }
            std::string::iterator endTagData() { return AllCharData.end(); }

            std::string::const_iterator beginName() const { return AllCharData.begin(); }
            std::string::const_iterator endName() const { return AllCharData.begin() + QueryNameLength; }
            std::string::const_iterator beginCigar() const { return endName(); }
            std::string::const_iterator endCigar() const { return beginCigar() + (NumCigarOperations * 4); }
            std::string::const_iterator beginSeq() const { return endCigar(); }
            std::string::const_iterator endSeq() const { return beginSeq() + ((QuerySequenceLength + 1)/2); }
            std::string::const_iterator beginQual() const { return endSeq(); }
            std::string::const_iterator endQual() const { return beginQual() + QuerySequenceLength; }
            std::string::const_iterator beginTagData() const { return endQual(); }
            std::string::const_iterator endTagData() const { return AllCharData.end(); }
        public:
            //accessors (with encoding)
            void setName(const std::string & name) { AllCharData.replace(beginName(), endName(), name.c_str(), name.size() + 1);  QueryNameLength = name.size()+1;}
            const std::string getName() const { return std::string(&AllCharData[0], QueryNameLength-1); }
            
            void setCigar(const std::vector<CigarOp> & cigar);
            const std::vector<CigarOp> getCigar() const;
            
            void setSeq(const std::string & seq);
            const std::string getSeq() const;
            
            void setQual(const std::string & qual);
            const std::string getQual() const;
            
            void setTagData(const std::string & data) { AllCharData.replace(beginTagData(), endTagData(), data); }
            const std::string getTagData() const { return std::string(beginTagData(), endTagData()); }
            const TagDataView getTagDataView() const { return TagDataView(&*beginTagData(), distance(beginTagData(), endTagData())); }
            
            uint32_t getQueryNameLength() const { return QueryNameLength; }
            uint32_t getQuerySequenceLength() const { return QuerySequenceLength; }
            uint32_t getNumCigarOperations() const { return NumCigarOperations; }
            
            void setData(const char * data, size_t data_len, uint32_t num_cigar, uint32_t seq_len, uint32_t name_len) { AllCharData.assign(data, data_len); NumCigarOperations = num_cigar; QuerySequenceLength = seq_len; QueryNameLength = name_len; }
            
            const std::string & getAllCharData() const { return AllCharData;}
        };
        const std::string & getBamEncodedStringData() const { return SupportData.getAllCharData(); }

        void setBamStringData(const char * data, size_t data_len, uint32_t num_cigar, uint32_t seq_len, uint32_t name_len) {
            SupportData.setData(data, data_len, num_cigar, seq_len, name_len);
        }
    protected:
        template <class> friend class BamDeserializer;
    public: //FIXME!!! this is unnecessarily public, and should only be available to friends (see above).
        const BamAlignmentSupportData & getSupportData() const { return SupportData; }
    private:
        mutable BamAlignmentSupportData SupportData;
        
        mutable std::string ErrorString; // mutable to allow updates even in logically const methods
        //! \endinternal
        
    public:
        std::string cigarString() const { return cigarToString(SupportData.getCigar()); }
    };
    
    // ---------------------------------------------------------
    // BamAlignment tag access methods
    
    /*! \fn bool AddTag(const std::string& tag, const std::string& type, const T& value)
     \brief Adds a field to the BAM tags.
     
     Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.
     
     \param[in] tag   2-character tag name
     \param[in] type  1-character tag type
     \param[in] value data to store
     \return \c true if the \b new tag was added successfully
     \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
     */
    template<typename T>
    inline bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const T& value) {
        
        std::string TagData = SupportData.getTagData();
        
        // check tag/type size
        if ( !IsValidSize(tag, type) ) {
            std::cerr << "Invalid size error" << std::endl;
            return false;
        }
        
        // check that storage type code is OK for T
        if ( !TagTypeHelper<T>::CanConvertTo(type.at(0)) ) {
            std::cerr << "Can convert to error - data provided for AddTag is probably wrong data type" << std::endl;
            assert(0);
            return false;
        }
        
        // localize the tag data
        char* pTagData = (char*)TagData.data();
        const unsigned int tagDataLength = TagData.size();
        unsigned int numBytesParsed = 0;
        
        // if tag already exists, return false
        // use EditTag explicitly instead
        if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
            std::cerr << "Find tag error" << std::endl;
            return false;
        }
        
        // otherwise, convert value to string
        union { T value; char valueBuffer[sizeof(T)]; } un;
        un.value = value;
        
        // copy original tag data to temp buffer
        const std::string newTag = tag + type;
        const size_t newTagDataLength = tagDataLength + newTag.size() + sizeof(T); // leave room for new T
        char * originalTagData = (char*) alloca(newTagDataLength);
        memcpy(originalTagData, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
        
        // append newTag
        strcat(originalTagData + tagDataLength, newTag.data());
        memcpy(originalTagData + tagDataLength + newTag.size(), un.valueBuffer, sizeof(T));
        
        // store temp buffer back in TagData
        const char* newTagData = (const char*)originalTagData;
        TagData.assign(newTagData, newTagDataLength);
        
        SupportData.setTagData(TagData);
        
        return true;
    }
    
    template<>
    inline bool BamAlignment::AddTag<std::string>(const std::string& tag,
                                                  const std::string& type,
                                                  const std::string& value)
    {
        std::string TagData = SupportData.getTagData();
        
        // check tag/type size
        if ( !IsValidSize(tag, type) ) {
            // TODO: set error string?
            return false;
        }
        
        // check that storage type code is OK for string
        if ( !TagTypeHelper<std::string>::CanConvertTo(type.at(0)) ) {
            // TODO: set error string?
            return false;
        }
        
        // localize the tag data
        char* pTagData = (char*)TagData.data();
        const unsigned int tagDataLength = TagData.size();
        unsigned int numBytesParsed = 0;
        
        // if tag already exists, return false
        // use EditTag explicitly instead
        if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
            // TODO: set error string?
            return false;
        }

        // otherwise, copy tag data to temp buffer
        const std::string newTag = tag + type + value;
        const size_t newTagDataLength = tagDataLength + newTag.size() + 1; // leave room for null-term
        char * originalTagData = (char *) alloca(newTagDataLength);
        memcpy(originalTagData, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
        
        // append newTag (removes original null-term, then appends newTag + null-term)
        strcat(originalTagData + tagDataLength, newTag.data());
        
        // store temp buffer back in TagData
        const char* newTagData = (const char*)originalTagData;
        TagData.assign(newTagData, newTagDataLength);

        SupportData.setTagData(TagData);
        return true;
    }
    
    /*! \fn template<typename T> bool AddTag(const std::string& tag, const std::vector<T>& values)
     \brief Adds a numeric array field to the BAM tags.
     
     Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.
     
     \param[in] tag    2-character tag name
     \param[in] values vector of data values to store
     \return \c true if the \b new tag was added successfully
     \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
     */
    template<typename T>
    inline bool BamAlignment::AddTag(const std::string& tag, const std::vector<T>& values) {
        
        std::string TagData = SupportData.getTagData();
        
        // check for valid tag name length
        if ( tag.size() != Constants::BAM_TAG_TAGSIZE )
            return false;
        
        // localize the tag data
        char* pTagData = (char*)TagData.data();
        const unsigned int tagDataLength = TagData.size();
        unsigned int numBytesParsed = 0;
        
        // if tag already exists, return false
        // use EditTag explicitly instead
        if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
            // TODO: set error string?
            return false;
        }
        
        // build new tag's base information
        char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
        memcpy( newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE );
        newTagBase[2] = Constants::BAM_TAG_TYPE_ARRAY;
        newTagBase[3] = TagTypeHelper<T>::TypeCode();
        
        // add number of array elements to newTagBase
        const int32_t numElements  = values.size();
        memcpy(newTagBase + 4, &numElements, sizeof(int32_t));
        
        // copy current TagData string to temp buffer, leaving room for new tag's contents
        const size_t newTagDataLength = tagDataLength +
        Constants::BAM_TAG_ARRAYBASE_SIZE +
        numElements*sizeof(T);
        char * originalTagData = (char *) alloca(newTagDataLength);
        memcpy(originalTagData, TagData.c_str(), tagDataLength+1); // '+1' for TagData's null-term
        
        // write newTagBase (removes old null term)
        strcat(originalTagData + tagDataLength, (const char*)newTagBase);
        
        // add vector elements to tag
        int elementsBeginOffset = tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE;
        for ( int i = 0 ; i < numElements; ++i ) {
            const T& value = values.at(i);
            memcpy(originalTagData + elementsBeginOffset + i*sizeof(T), &value, sizeof(T));
        }
        
        // store temp buffer back in TagData
        const char* newTagData = (const char*)originalTagData;
        TagData.assign(newTagData, newTagDataLength);
        
        SupportData.setTagData(TagData);

        return true;
    }
    
    /*! \fn template<typename T> bool EditTag(const std::string& tag, const std::string& type, const T& value)
     \brief Edits a BAM tag field.
     
     If \a tag does not exist, a new entry is created.
     
     \param tag[in]   2-character tag name
     \param type[in]  1-character tag type (must be "Z" or "H")
     \param value[in] new data value
     
     \return \c true if the tag was modified/created successfully
     
     \sa BamAlignment::RemoveTag()
     \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
     */
    template<typename T>
    inline bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const T& value) {
        
        // remove existing tag if present, then append tag with new value
        if ( HasTag(tag) )
            RemoveTag(tag);
        return AddTag(tag, type, value);
    }
    
    /*! \fn template<typename T> bool EditTag(const std::string& tag, const std::vector<T>& values)
     \brief Edits a BAM tag field containing a numeric array.
     
     If \a tag does not exist, a new entry is created.
     
     \param tag[in]   2-character tag name
     \param value[in] vector of data values
     
     \return \c true if the tag was modified/created successfully
     \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
     */
    template<typename T>
    inline bool BamAlignment::EditTag(const std::string& tag, const std::vector<T>& values) {
        
        // remove existing tag if present, then append tag with new values
        if ( HasTag(tag) )
            RemoveTag(tag);
        return AddTag(tag, values);
    }
    
    
    /*! \fn template<typename T> bool GetTag(const std::string& tag, T& destination) const
     \brief Retrieves the value associated with a BAM tag.
     
     \param tag[in]          2-character tag name
     \param destination[out] retrieved value
     \return \c true if found
     */
    template<typename T>
    inline bool BamAlignment::GetTag(const std::string& tag, T& destination) const {
        
        const TagDataView TagData = SupportData.getTagDataView();
        
        // skip if no tags present
        if ( TagData.empty() ) {
            // TODO: set error string?
            return false;
        }
        
        // localize the tag data
        char* pTagData = (char*)TagData.data();
        const unsigned int tagDataLength = TagData.size();
        unsigned int numBytesParsed = 0;
        
        // return failure if tag not found
        if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
            // TODO: set error string?
            return false;
        }
        
        // fetch data type
        const char type = *(pTagData - 1);
        if ( !TagTypeHelper<T>::CanConvertFrom(type) ) {
            // TODO: set error string ?
            return false;
        }
        
        // determine data length
        int destinationLength = 0;
        switch ( type ) {
                
                // 1 byte data
            case (Constants::BAM_TAG_TYPE_ASCII) :
            case (Constants::BAM_TAG_TYPE_INT8)  :
            case (Constants::BAM_TAG_TYPE_UINT8) :
                destinationLength = 1;
                break;
                
                // 2 byte data
            case (Constants::BAM_TAG_TYPE_INT16)  :
            case (Constants::BAM_TAG_TYPE_UINT16) :
                destinationLength = 2;
                break;
                
                // 4 byte data
            case (Constants::BAM_TAG_TYPE_INT32)  :
            case (Constants::BAM_TAG_TYPE_UINT32) :
            case (Constants::BAM_TAG_TYPE_FLOAT)  :
                destinationLength = 4;
                break;
                
                // var-length types not supported for numeric destination
            case (Constants::BAM_TAG_TYPE_STRING) :
            case (Constants::BAM_TAG_TYPE_HEX)    :
            case (Constants::BAM_TAG_TYPE_ARRAY)  :
                SetErrorString("BamAlignment::GetTag",
                               "cannot store variable length tag data into a numeric destination");
                return false;
                
                // unrecognized tag type
            default:
                const std::string message = std::string("invalid tag type: ") + type;
                SetErrorString("BamAlignment::GetTag", message);
                return false;
        }
        
        // store data in destination
        destination = 0;
        memcpy(&destination, pTagData, destinationLength);
        
        // return success
        return true;
    }
    
    template<>
    inline bool BamAlignment::GetTag<std::string>(const std::string& tag,
                                                  std::string& destination) const
    {
        const TagDataView TagData = SupportData.getTagDataView();
        
        // skip if no tags present
        if ( TagData.empty() ) {
            // TODO: set error string?
            return false;
        }
        
        // localize the tag data
        char* pTagData = (char*)TagData.data();
        const unsigned int tagDataLength = TagData.size();
        unsigned int numBytesParsed = 0;
        
        // return failure if tag not found
        if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
            // TODO: set error string?
            return false;
        }
        
        // otherwise copy data into destination
        const unsigned int dataLength = strlen(pTagData);
        destination.clear();
        destination.resize(dataLength);
        memcpy( (char*)destination.data(), pTagData, dataLength );
        
        // return success
        return true;
    }
    
    /*! \fn template<typename T> bool GetTag(const std::string& tag, std::vector<T>& destination) const
     \brief Retrieves the numeric array associated with a BAM tag.
     
     \param tag[in]          2-character tag name
     \param destination[out] retrieved values
     \return \c true if found
     */
    template<typename T>
    inline bool BamAlignment::GetTag(const std::string& tag, std::vector<T>& destination) const {
        
        const TagDataView TagData = SupportData.getTagDataView();
        
        // skip if no tags present
        if ( TagData.empty() ) {
            // TODO: set error string?
            return false;
        }
        
        // localize the tag data
        char* pTagData = (char*)TagData.data();
        const unsigned int tagDataLength = TagData.size();
        unsigned int numBytesParsed = 0;
        
        // return false if tag not found
        if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
            // TODO: set error string?
            return false;
        }
        
        // check that tag is array type
        const char tagType = *(pTagData - 1);
        if ( tagType != Constants::BAM_TAG_TYPE_ARRAY ) {
            SetErrorString("BamAlignment::GetTag", "cannot store a non-array tag in array destination");
            return false;
        }
        
        // fetch element type
        const char elementType = *pTagData;
        if ( !TagTypeHelper<T>::CanConvertFrom(elementType) ) {
            // TODO: set error string ?
            return false;
        }
        ++pTagData;
        
        // calculate length of each element in tag's array
        int elementLength = 0;
        switch ( elementType ) {
            case (Constants::BAM_TAG_TYPE_ASCII) :
            case (Constants::BAM_TAG_TYPE_INT8)  :
            case (Constants::BAM_TAG_TYPE_UINT8) :
                elementLength = sizeof(uint8_t);
                break;
                
            case (Constants::BAM_TAG_TYPE_INT16)  :
            case (Constants::BAM_TAG_TYPE_UINT16) :
                elementLength = sizeof(uint16_t);
                break;
                
            case (Constants::BAM_TAG_TYPE_INT32)  :
            case (Constants::BAM_TAG_TYPE_UINT32) :
            case (Constants::BAM_TAG_TYPE_FLOAT)  :
                elementLength = sizeof(uint32_t);
                break;
                
                // var-length types not supported for numeric destination
            case (Constants::BAM_TAG_TYPE_STRING) :
            case (Constants::BAM_TAG_TYPE_HEX)    :
            case (Constants::BAM_TAG_TYPE_ARRAY)  :
                SetErrorString("BamAlignment::GetTag",
                               "invalid array data, variable-length elements are not allowed");
                return false;
                
                // unknown tag type
            default:
                const std::string message = std::string("invalid array element type: ") + elementType;
                SetErrorString("BamAlignment::GetTag", message);
                return false;
        }
        
        // get number of elements
        int32_t numElements;
        memcpy(&numElements, pTagData, sizeof(int32_t));
        pTagData += 4;
        destination.clear();
        destination.reserve(numElements);
        
        // read in elements
        T value;
        for ( int i = 0 ; i < numElements; ++i ) {
            memcpy(&value, pTagData, sizeof(T));
            pTagData += sizeof(T);
            destination.push_back(value);
        }
        
        // return success
        return true;
    }
    
    typedef std::vector<BamAlignment> BamAlignmentVector;
    
} // namespace BamTools

#endif // BAMALIGNMENT_H
