// ***************************************************************************
// BamIndexFactory_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides interface for generating BamIndex implementations
// ***************************************************************************

#include "api/internal/index/BamIndexFactory_p.h"
#include "api/internal/index/BamStandardIndex_p.h"
using namespace BamTools;
using namespace BamTools::Internal;
using namespace std;

// generates index filename from BAM filename (depending on requested type)
// if type is unknown, returns empty string
const string BamIndexFactory::CreateIndexFilename(const string& bamFilename,
                                                  const BamIndex::IndexType& type)
{
    return ( bamFilename + BamStandardIndex::Extension() );
}

// creates a new BamIndex object, depending on extension of @indexFilename
BamIndex* BamIndexFactory::CreateIndexFromFilename(const string& indexFilename, BamReaderPrivate* reader) {

    // get file extension from index filename, including dot (".EXT")
    // if can't get file extension, return null index
    const string extension = FileExtension(indexFilename);
    if ( extension.empty() )
        return 0;

    // create index based on extension
    return new BamStandardIndex(reader);
}

// creates a new BamIndex, object of requested @type
BamIndex* BamIndexFactory::CreateIndexOfType(const BamIndex::IndexType& type,
                                             BamReaderPrivate* reader)
{
    return new BamStandardIndex(reader);
}

// retrieves file extension (including '.')
const string BamIndexFactory::FileExtension(const string& filename) {

    // if filename cannot contain valid path + extension, return empty string
    if ( filename.empty() || filename.length() <= 4 )
        return string();

    // look for last dot in filename
    const size_t lastDotPosition = filename.find_last_of('.');

    // if none found, return empty string
    if ( lastDotPosition == string::npos )
        return string();

    // return substring from last dot position
    return filename.substr(lastDotPosition);
}

// returns name of existing index file that corresponds to @bamFilename
// will defer to @preferredType if possible, if not will attempt to load any supported type
// returns empty string if not found
const string BamIndexFactory::FindIndexFilename(const string& bamFilename,
                                                const BamIndex::IndexType& preferredType)
{
    // skip if BAM filename provided is empty
    if ( bamFilename.empty() )
        return string();

    // try to find index of preferred type first
    // return index filename if found
    string indexFilename = CreateIndexFilename(bamFilename, preferredType);
    if ( !indexFilename.empty() )
        return indexFilename;

    // couldn't find preferred type, try the other supported types
    // return index filename if found
    if ( preferredType != BamIndex::STANDARD ) {
        indexFilename = CreateIndexFilename(bamFilename, BamIndex::STANDARD);
        if ( !indexFilename.empty() )
            return indexFilename;
    }

    // otherwise couldn't find any index matching this filename
    return string();
}
