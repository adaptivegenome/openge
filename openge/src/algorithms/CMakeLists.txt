set(ALGORITHMS_DIR "algorithms")

set(OPENGE_ALGORITHMS_SOURCES
  ${ALGORITHMS_DIR}/algorithm_module.cpp
  ${ALGORITHMS_DIR}/filter.cpp
  ${ALGORITHMS_DIR}/file_reader.cpp
  ${ALGORITHMS_DIR}/file_writer.cpp
  ${ALGORITHMS_DIR}/local_realignment.cpp
  ${ALGORITHMS_DIR}/mark_duplicates.cpp
  ${ALGORITHMS_DIR}/measure_coverage.cpp
  ${ALGORITHMS_DIR}/print_history.cpp
  ${ALGORITHMS_DIR}/read_sorter.cpp
  ${ALGORITHMS_DIR}/repeatseq.cpp
  ${ALGORITHMS_DIR}/sorted_merge.cpp
  ${ALGORITHMS_DIR}/split_by_chromosome.cpp
  ${ALGORITHMS_DIR}/statistics.cpp
  PARENT_SCOPE)
