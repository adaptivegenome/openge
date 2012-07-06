# some variables
project(openge_unit_tests)

#add_subdirectory( commands )
add_subdirectory( algorithms )
add_subdirectory( util )

set( OPENGE_TEST_SOURCES
     main.cpp
     ${OPENGE_ALGORITHMS_TEST_SOURCES}
     ${OPENGE_UTIL_TEST_SOURCES} )
set( TEST_EXE  openge_unit_tests )


# gtest and openge .h files location
include_directories(${gtest_SOURCE_DIR}/include
                    ${gtest_SOURCE_DIR}
                    ${OpenGE_SOURCE_DIR}/openge/src )


# final targets
add_executable(${TEST_EXE} ${OPENGE_TEST_SOURCES})
target_link_libraries(${TEST_EXE} gtest openge-lib)
