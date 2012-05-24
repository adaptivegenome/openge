set(ALGORITHMS_DIR "algorithms")

set(OPENGE_ALGORITHMS_SOURCES
  ${ALGORITHMS_DIR}/algorithm_module.cpp
  ${ALGORITHMS_DIR}/algorithm_module_adaptor.cpp
  ${ALGORITHMS_DIR}/filter.cpp
  ${ALGORITHMS_DIR}/file_reader.cpp
  ${ALGORITHMS_DIR}/file_writer.cpp
  ${ALGORITHMS_DIR}/local_realignment.cpp
  ${ALGORITHMS_DIR}/mark_duplicates.cpp
  ${ALGORITHMS_DIR}/read_sorter.cpp
  ${ALGORITHMS_DIR}/statistics.cpp
  PARENT_SCOPE)
