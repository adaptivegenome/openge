add_subdirectory( commands )
add_subdirectory( algorithms )
add_subdirectory( util )

set( OPENGE_SOURCES main.cpp ${OPENGE_ALGORITHMS_SOURCES} ${OPENGE_COMMANDS_SOURCES} ${OPENGE_UTIL_SOURCES})

set (OPENGE_VERSION 0.3)
set (OPENGE_BUILD_TYPE dev) # dev or final

configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/openge_constants.h.in" "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/openge_constants.h" )
include_directories("${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}")

add_executable(openge ${OPENGE_SOURCES})
target_link_libraries(openge ${Boost_LIBRARIES} ${BAMTOOLS_LIBRARIES} z pthread)
