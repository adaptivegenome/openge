project(OpenGE)

cmake_minimum_required(VERSION 2.8)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

enable_testing()
add_subdirectory(bamtools)
add_subdirectory(openge)
add_dependencies(openge BamTools-static)

link_directories(${LIBRARY_OUTPUT_PATH})
