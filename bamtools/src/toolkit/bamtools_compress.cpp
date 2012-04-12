// ***************************************************************************
// bamtools_Compress.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 October 2011
// ---------------------------------------------------------------------------
// Compresss between BAM and a number of other formats
// ***************************************************************************

#include "bamtools_compress.h"

#include <api/BamConstants.h>
#include <api/BamWriter.h>
#include <api/SamReader.h>

#include <utils/bamtools_options.h>
using namespace BamTools;

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

// ---------------------------------------------
// CompressSettings implementation

struct CompressTool::CompressSettings {
  
  // flag
  bool HasInput;
  bool HasOutput;
  bool HasCompressionLevel;
  
  unsigned int compression_level;
  string InputFilename;
  string OutputFilename;
  
  
  // constructor
  CompressSettings(void)
  : HasInput(false)
  , HasOutput(false)
  , HasCompressionLevel(false)
  { } 
};    

// ---------------------------------------------
// CompressToolPrivate implementation  

struct CompressTool::CompressToolPrivate {
  
  // ctor & dtor
public:
  CompressToolPrivate(CompressTool::CompressSettings* settings)
  : m_settings(settings)
  , m_out(cout.rdbuf())
  { }
  
  ~CompressToolPrivate(void) { }
  
  // interface
public:
  bool Run(void);
  
  // internal methods
private:
  
  // data members
private: 
  CompressTool::CompressSettings* m_settings;
  RefVector m_references;
  ostream m_out;
};

bool CompressTool::CompressToolPrivate::Run(void) {
  
  // ------------------------------------
  // initialize conversion input/output
  
  // set to default input if none provided
  if ( !m_settings->HasInput ) 
    m_settings->InputFilename = Options::StandardIn();
  if ( !m_settings->HasOutput ) 
    m_settings->InputFilename = Options::StandardOut();
  
  SamReader reader;
  reader.Open(m_settings->InputFilename);
  
  BamWriter writer;
  writer.Open( m_settings->OutputFilename, reader.GetHeader(), reader.GetRefData());
  
  int alignment_ct = 0;
  
  while(true) {
    BamAlignment alignment;
    
    if(!reader.GetNextAlignment(alignment))
      break;
    writer.SaveAlignment(alignment);
    alignment_ct++;
    
    //progress  indicator
    //if(alignment_ct % 500000 == 0)
    //  cerr << ".";
  }
  
  reader.Close();
  writer.Close();

  return true;
}


// ---------------------------------------------
// CompressTool implementation

CompressTool::CompressTool(void)
: AbstractTool()
, m_settings(new CompressSettings)
, m_impl(0)
{
  // set program details
  Options::SetProgramInfo("bamtools compress", "Compresses a SAM file to a BAM file", "[-in <filename>] [-out <filename>] ");
  
  // set up options 
  OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");
  Options::AddValueOption("-in",     "BAM filename", "the input BAM file", "", m_settings->HasInput,   m_settings->InputFilename,     IO_Opts, Options::StandardIn());
  Options::AddValueOption("-out",    "BAM filename", "the output BAM file",   "", m_settings->HasOutput,  m_settings->OutputFilename, IO_Opts, Options::StandardOut());
}

CompressTool::~CompressTool(void) {
  
  delete m_settings;
  m_settings = 0;
  
  delete m_impl;
  m_impl = 0;
}

int CompressTool::Help(void) {
  Options::DisplayHelp();
  return 0;
}

int CompressTool::Run(int argc, char* argv[]) {
  
  // parse command line arguments
  Options::Parse(argc, argv, 1);
  
  // initialize CompressTool with settings
  m_impl = new CompressToolPrivate(m_settings);
  
  // run CompressTool, return success/fail
  if ( m_impl->Run() ) 
    return 0;
  else 
    return 1;
}