// ***************************************************************************
// BamAux.h (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides data structures & utility methods that are used throughout the API.
// ***************************************************************************

#ifndef BAMAUX_H
#define BAMAUX_H

#include <cstring>
#include <fstream> 
#include <iostream>
#include <string>
#include <vector>
#include <stdint.h>

/*! \file BamAux.h

    Provides data structures & utility methods that are used throughout the API.
*/

/*! \namespace BamTools
    \brief Contains all BamTools classes & methods.

    The BamTools API contained in this namespace contains classes and methods
    for reading, writing, and manipulating BAM alignment files.
*/
namespace BamTools {

// ----------------------------------------------------------------
// BamRegion

/*! \struct BamTools::BamRegion
    \brief Represents a sequential genomic region

    Allowed to span multiple (sequential) references.

    \warning BamRegion now represents a zero-based, HALF-OPEN interval.
    In previous versions of BamTools (0.x & 1.x) all intervals were treated
    as zero-based, CLOSED.
*/
struct BamRegion {
  
    int LeftRefID;      //!< reference ID for region's left boundary
    int LeftPosition;   //!< position for region's left boundary
    int RightRefID;     //!< reference ID for region's right boundary
    int RightPosition;  //!< position for region's right boundary
    
    //! constructor
    BamRegion(const int& leftID   = -1, 
              const int& leftPos  = -1,
              const int& rightID  = -1,
              const int& rightPos = -1)
        : LeftRefID(leftID)
        , LeftPosition(leftPos)
        , RightRefID(rightID)
        , RightPosition(rightPos)
    { }
    
    //! copy constructor
    BamRegion(const BamRegion& other)
        : LeftRefID(other.LeftRefID)
        , LeftPosition(other.LeftPosition)
        , RightRefID(other.RightRefID)
        , RightPosition(other.RightPosition)
    { }
    
    //! Clears region boundaries
    void clear(void) {
        LeftRefID  = -1; LeftPosition  = -1;
        RightRefID = -1; RightPosition = -1;
    }

    //! Returns true if region has a left boundary
    bool isLeftBoundSpecified(void) const {
        return ( LeftRefID >= 0 && LeftPosition >= 0 );
    }

    //! Returns true if region boundaries are not defined
    bool isNull(void) const {
        return ( !isLeftBoundSpecified() && !isRightBoundSpecified() );
    }

    //! Returns true if region has a right boundary
    bool isRightBoundSpecified(void) const {
        return ( RightRefID >= 0 && RightPosition >= 1 );
    }
};

// ----------------------------------------------------------------
// General utility methods

/*! \fn bool SystemIsBigEndian(void)
    \brief checks host architecture's byte order
    \return \c true if system uses big-endian ordering
*/
inline bool SystemIsBigEndian(void) {
   const uint16_t one = 0x0001;
   return ((*(char*) &one) == 0 );
}

/*! \fn void PackUnsignedInt(char* buffer, unsigned int value)
    \brief stores unsigned integer value in a byte buffer

    \param[out] buffer destination buffer
    \param[in]  value  value to 'pack' in buffer
*/
inline void PackUnsignedInt(char* buffer, unsigned int value) {
    buffer[0] = (char)value;
    buffer[1] = (char)(value >> 8);
    buffer[2] = (char)(value >> 16);
    buffer[3] = (char)(value >> 24);
}

/*! \fn void PackUnsignedShort(char* buffer, unsigned short value)
    \brief stores unsigned short integer value in a byte buffer

    \param[out] buffer destination buffer
    \param[in]  value  value to 'pack' in buffer
*/
inline void PackUnsignedShort(char* buffer, unsigned short value) {
    buffer[0] = (char)value;
    buffer[1] = (char)(value >> 8);
}

/*! \fn double UnpackDouble(const char* buffer)
    \brief reads a double value from byte buffer

    \param[in] buffer source byte buffer
    \return the (double) value read from the buffer
*/
inline double UnpackDouble(const char* buffer) {
    union { double value; unsigned char valueBuffer[sizeof(double)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    un.valueBuffer[4] = buffer[4];
    un.valueBuffer[5] = buffer[5];
    un.valueBuffer[6] = buffer[6];
    un.valueBuffer[7] = buffer[7];
    return un.value;
}

/*! \fn double UnpackDouble(char* buffer)
    \brief reads a double value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (double) value read from the buffer
*/
inline double UnpackDouble(char* buffer) {
    return UnpackDouble( (const char*)buffer );
}

/*! \fn double UnpackFloat(const char* buffer)
    \brief reads a float value from byte buffer

    \param[in] buffer source byte buffer
    \return the (float) value read from the buffer
*/
inline float UnpackFloat(const char* buffer) {
    union { float value; unsigned char valueBuffer[sizeof(float)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

/*! \fn double UnpackFloat(char* buffer)
    \brief reads a float value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (float) value read from the buffer
*/
inline float UnpackFloat(char* buffer) {
    return UnpackFloat( (const char*)buffer );
}

/*! \fn signed int UnpackSignedInt(const char* buffer)
    \brief reads a signed integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (signed int) value read from the buffer
*/
inline signed int UnpackSignedInt(const char* buffer) {
    union { signed int value; unsigned char valueBuffer[sizeof(signed int)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

/*! \fn signed int UnpackSignedInt(char* buffer)
    \brief reads a signed integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (signed int) value read from the buffer
*/
inline signed int UnpackSignedInt(char* buffer) {
    return UnpackSignedInt( (const char*) buffer );
}

/*! \fn signed short UnpackSignedShort(const char* buffer)
    \brief reads a signed short integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (signed short) value read from the buffer
*/
inline signed short UnpackSignedShort(const char* buffer) {
    union { signed short value; unsigned char valueBuffer[sizeof(signed short)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    return un.value;
}

/*! \fn signed short UnpackSignedShort(char* buffer)
    \brief reads a signed short integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (signed short) value read from the buffer
*/
inline signed short UnpackSignedShort(char* buffer) {
    return UnpackSignedShort( (const char*)buffer );
}

/*! \fn unsigned int UnpackUnsignedInt(const char* buffer)
    \brief reads an unsigned integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (unsigned int) value read from the buffer
*/
inline unsigned int UnpackUnsignedInt(const char* buffer) {
    union { unsigned int value; unsigned char valueBuffer[sizeof(unsigned int)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

/*! \fn unsigned int UnpackUnsignedInt(char* buffer)
    \brief reads an unsigned integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (unsigned int) value read from the buffer
*/
inline unsigned int UnpackUnsignedInt(char* buffer) {
    return UnpackUnsignedInt( (const char*)buffer );
}

/*! \fn unsigned short UnpackUnsignedShort(const char* buffer)
    \brief reads an unsigned short integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (unsigned short) value read from the buffer
*/
inline unsigned short UnpackUnsignedShort(const char* buffer) {
    union { unsigned short value; unsigned char valueBuffer[sizeof(unsigned short)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    return un.value;
}

/*! \fn unsigned short UnpackUnsignedShort(char* buffer)
    \brief reads an unsigned short integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (unsigned short) value read from the buffer
*/
inline unsigned short UnpackUnsignedShort(char* buffer) {
    return UnpackUnsignedShort( (const char*)buffer );
}
} // namespace BamTools

#endif // BAMAUX_H
