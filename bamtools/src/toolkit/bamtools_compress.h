// ***************************************************************************
// bamtools_compress.h (c) 2012 Lee C. Baker
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 15 March 2012
// ---------------------------------------------------------------------------
// Compresss SAM files to BAM
// ***************************************************************************

#ifndef BAMTOOLS_COMPRESS_H
#define BAMTOOLS_COMPRESS_H

#include "bamtools_tool.h"

namespace BamTools { 
  
class CompressTool : public AbstractTool {
  
    public:
        CompressTool(void);
        ~CompressTool(void);

    public:
        int Help(void);
        int Run(int argc, char* argv[]); 
        
    private: 
        struct CompressSettings;
        CompressSettings* m_settings;
        
        struct CompressToolPrivate;
        CompressToolPrivate* m_impl;
};
  
} // namespace BamTools

#endif // BAMTOOLS_COMPRESS_H
