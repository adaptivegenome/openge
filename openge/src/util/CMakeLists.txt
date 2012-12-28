set(UTIL_DIR "util")

set(OPENGE_UTIL_SOURCES 
  ${UTIL_DIR}/bam_deserializer.h
  ${UTIL_DIR}/bam_header.h
  ${UTIL_DIR}/bam_header.cpp
  ${UTIL_DIR}/bam_serializer.h
  ${UTIL_DIR}/bgzf_input_stream.h
  ${UTIL_DIR}/bgzf_input_stream.cpp
  ${UTIL_DIR}/bgzf_output_stream.h
  ${UTIL_DIR}/bgzf_output_stream.cpp
  ${UTIL_DIR}/bpipe.h
  ${UTIL_DIR}/bpipe.cpp
  ${UTIL_DIR}/fasta_reader.h
  ${UTIL_DIR}/fasta_reader.cpp
  ${UTIL_DIR}/fastq_writer.h
  ${UTIL_DIR}/fastq_writer.cpp
  ${UTIL_DIR}/file_io.h
  ${UTIL_DIR}/oge_read.h
  ${UTIL_DIR}/oge_read.cpp
  ${UTIL_DIR}/picard_structures.h
  ${UTIL_DIR}/picard_structures.cpp
  ${UTIL_DIR}/read_stream_reader.h
  ${UTIL_DIR}/read_stream_reader.cpp
  ${UTIL_DIR}/sam_reader.h
  ${UTIL_DIR}/sam_reader.cpp
  ${UTIL_DIR}/sam_writer.h
  ${UTIL_DIR}/sam_writer.cpp
  ${UTIL_DIR}/sequential_reader_cache.h
  ${UTIL_DIR}/thread_pool.h
  ${UTIL_DIR}/thread_pool.cpp

  ${UTIL_DIR}/bamtools/BamAlignment.h
  ${UTIL_DIR}/bamtools/BamAlignment.cpp
  ${UTIL_DIR}/bamtools/BamAux.h
  ${UTIL_DIR}/bamtools/BamConstants.h
  PARENT_SCOPE )

set(OPENGE_UTIL_GATK_SOURCES 
  ${UTIL_DIR}/gatk/AlignmentBlock.h
  ${UTIL_DIR}/gatk/AlignmentUtils.h
  ${UTIL_DIR}/gatk/AlignmentUtils.cpp
  ${UTIL_DIR}/gatk/Allele.h
  ${UTIL_DIR}/gatk/Allele.cpp
  ${UTIL_DIR}/gatk/BaseUtils.h
  ${UTIL_DIR}/gatk/ConstrainedMateFixingManager.h
  ${UTIL_DIR}/gatk/ConstrainedMateFixingManager.cpp
  ${UTIL_DIR}/gatk/GATKFeature.h
  ${UTIL_DIR}/gatk/GenomeLoc.h
  ${UTIL_DIR}/gatk/GenomeLoc.cpp
  ${UTIL_DIR}/gatk/GenomeLocParser.h
  ${UTIL_DIR}/gatk/GenomeLocParser.cpp
  ${UTIL_DIR}/gatk/RODMetaDataContainer.h
  ${UTIL_DIR}/gatk/ReadMetaDataTracker.h
  ${UTIL_DIR}/gatk/ReadMetaDataTracker.cpp
  ${UTIL_DIR}/gatk/ReferenceContext.h
  ${UTIL_DIR}/gatk/ReferenceContext.cpp
  ${UTIL_DIR}/gatk/SequenceUtil.h
  ${UTIL_DIR}/gatk/SequenceUtil.cpp
  ${UTIL_DIR}/gatk/VariantContext.h
  ${UTIL_DIR}/gatk/VariantContext.cpp
  PARENT_SCOPE)

