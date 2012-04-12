// ***************************************************************************
// bamtools_merge.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// Modified Lee C. Baker
// ---------------------------------------------------------------------------
// Last modified: 7 March 2012
// ---------------------------------------------------------------------------
// Merges multiple BAM files into one
// ***************************************************************************

#ifndef BAMTOOLS_MERGESORT_H
#define BAMTOOLS_MERGESORT_H

#include "bamtools_tool.h"

namespace BamTools {

class MergeSortTool : public AbstractTool {

    public:
        MergeSortTool(void);
        ~MergeSortTool(void);

    public:
        int Help(void);
        int Run(int argc, char* argv[]);

    private:
        struct MergeSortSettings;
        MergeSortSettings* m_settings;

        struct MergeSortToolPrivate;
        MergeSortToolPrivate* m_impl;
};

} // namespace BamTools

#endif // BAMTOOLS_MERGESORT_H
