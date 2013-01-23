#ifndef OGE_BAM_INDEX_H
#define OGE_BAM_INDEX_H
/*********************************************************************
 *
 * bam_index.h: Generate an index for a BAM file.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 28 Dec 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial
 * Purpose License. A copy of this license has been provided in
 * the openge/ directory.
 *
 *********************************************************************/

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <stdint.h>
#include "bam_header.h"
#include "oge_read.h"

#ifndef UINT64_MAX
#define UINT64_MAX        18446744073709551615ULL
#endif

class BgzfOutputStream;
class BgzfInputStream;

class BamIndex {
    typedef struct __metadata_t{
        uint64_t num_mapped_reads, num_unmapped_reads;
        uint64_t read_start_position, read_stop_position;
        __metadata_t()
        : num_mapped_reads(0)
        , num_unmapped_reads(0)
        , read_start_position(UINT64_MAX)
        , read_stop_position(0)
        {}
    } metadata_t;
    std::vector<metadata_t> metadata;
	class BamIndexSequence {
        
		class BamIndexBin {
			std::vector<std::pair<uint64_t, uint64_t> > chunks;
		public:
            BamIndexBin() {}
            BamIndexBin(uint64_t unmapped_reads, uint64_t mapped_reads, uint64_t data_start, uint64_t data_stop);  //special constructor for samtools' undocumented metadata bin :(
			void addRead(int start_pos, uint64_t file_start, uint64_t file_stop);
            void read(std::ifstream & stream);
			void write(std::ofstream & stream) const;
            void remap(BgzfOutputStream * remapper_stream);
		};
        std::vector<uint64_t> linear_index;
        std::map<uint32_t, BamIndexBin *> bins;
        void fillMissing();
	public:
		BamIndexSequence(const BamSequenceRecord & record);
        void setMetadataFrame(uint64_t unmapped_reads, uint64_t mapped_reads, uint64_t data_start, uint64_t data_stop);
		void addRead(const OGERead * read, int end_pos, int bin, uint64_t file_start, uint64_t file_stop);
        void read(std::ifstream & stream);
		void write(std::ofstream & stream) const;
        void remap(BgzfOutputStream * remapper_stream);
	};
    uint64_t num_coordless_reads;
	std::vector<BamIndexSequence *> sequences;
public:
	BamIndex(const BamHeader & h);
	~BamIndex();
	void addRead(const OGERead * read, int end_pos, int bin, uint64_t file_start, uint64_t file_stop);
    void readFile(const std::string & filename, BgzfInputStream * remapper_stream);
	void writeFile(const std::string & filename, BgzfOutputStream * remapper_stream) const;
};

#endif
