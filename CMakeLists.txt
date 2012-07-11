project(OpenGE)

cmake_minimum_required(VERSION 2.8)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

add_subdirectory(bamtools)
add_subdirectory( gtest-1.6.0 )
add_subdirectory(openge)
add_dependencies(openge BamTools-static)

link_directories(${LIBRARY_OUTPUT_PATH})
