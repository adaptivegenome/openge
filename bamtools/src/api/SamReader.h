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
#include <queue>
#include "api/BamParallelismSettings.h"
using namespace std;
using namespace BamTools;

class BamThreadPool;

class SamReader;

struct SamLine_t {
    string line;
    BamAlignment * al;
    SamReader * reader;
    bool parsed;
    SamLine_t() : al(NULL), reader(NULL), parsed(false) {}
};

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
    BamAlignment * GetNextAlignment();

  // returns the current file's header data
  SamHeader GetHeader(void) const;
  // get reference data
    const RefVector & GetReferenceData(void);
    //read a single line of a SAM file
    BamAlignment * ParseAlignment(const string & line_s);
protected:
  // retrieves BAM alignment under file pointer
    // (does no overlap checking or character data parsing)
    BamAlignment * LoadNextAlignment();

  // retrieves header text from SAM file
  void LoadHeaderData(void);
    

  vector<BamAlignment> alignments;
  std::queue<SamLine_t *> jobs;
  BamThreadPool * pool;
  ifstream file;
  SamHeader header;
  RefVector m_refData;
  bool loaded;
  bool finished;
};

#endif
