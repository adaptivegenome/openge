add_subdirectory( commands )
add_subdirectory( algorithms )
add_subdirectory( util )

set( OPENGE_SOURCES main.cpp ${OPENGE_ALGORITHMS_SOURCES} ${OPENGE_COMMANDS_SOURCES} ${OPENGE_UTIL_SOURCES})

add_executable(openge ${OPENGE_SOURCES})
target_link_libraries(openge ${Boost_LIBRARIES} ${BAMTOOLS_LIBRARIES} z pthread)
