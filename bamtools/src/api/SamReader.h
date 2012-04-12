//
//  SamReader.h
//  BamTools
//
//  Created by Lee Baker on 3/16/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#ifndef BamTools_SamReader_h
#define BamTools_SamReader_h
#include <string>
#include <vector>
#include <api/BamAlignment.h>
#include <api/SamHeader.h>
#include <iostream>
using namespace std;
using namespace BamTools;

// SamReader is capable of sequentially reading a SAM file. It doesn't support
// most of the features that BamReader does, only enough to support converting SAM
// files to the BAM format.
class SamReader
{
public:
  SamReader();
  bool Open(const string & filename);
  bool Close();
  bool IsLoaded();
  // retrieves next available alignment
  bool GetNextAlignment(BamAlignment& alignment);

  // returns the current file's header data
  SamHeader GetHeader(void) const;
  // get reference data
  const RefVector & GetRefData(void);
protected:
  // retrieves BAM alignment under file pointer
  // (does no overlap checking or character data parsing)
  bool LoadNextAlignment(BamAlignment& alignment);

  // retrieves header text from SAM file
  void LoadHeaderData(void);

  vector<BamAlignment> alignments;
  ifstream file;
  SamHeader header;
  RefVector m_refData;
  bool loaded;
};

#endif
