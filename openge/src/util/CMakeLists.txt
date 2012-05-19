set(UTIL_DIR "util")

set(OPENGE_UTIL_SOURCES 
  ${UTIL_DIR}/fasta_reader.cpp
  ${UTIL_DIR}/picard_structures.cpp
  ${UTIL_DIR}/sam_reader.cpp
  ${UTIL_DIR}/thread_pool.cpp

  ${UTIL_DIR}/gatk/AlignmentUtils.cpp
  ${UTIL_DIR}/gatk/Allele.cpp
  ${UTIL_DIR}/gatk/BaseUtils.cpp
  ${UTIL_DIR}/gatk/ConstrainedMateFixingManager.cpp
  ${UTIL_DIR}/gatk/GenomeLoc.cpp
  ${UTIL_DIR}/gatk/GenomeLocParser.cpp
  ${UTIL_DIR}/gatk/ReadMetaDataTracker.cpp
  ${UTIL_DIR}/gatk/ReferenceContext.cpp
  ${UTIL_DIR}/gatk/SequenceUtil.cpp
  ${UTIL_DIR}/gatk/VariantContext.cpp
  PARENT_SCOPE)

