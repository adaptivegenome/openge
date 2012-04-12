add_subdirectory( commands)
add_subdirectory( util)

set( OPENGE_SOURCES main.cpp ${OPENGE_COMMANDS_SOURCES} ${OPENGE_UTIL_SOURCES})

message("BAMTOOLS_LIBRARIES " ${BAMTOOLS_LIBRARIES})
add_executable(openge ${OPENGE_SOURCES})
target_link_libraries(openge ${Boost_LIBRARIES} ${BAMTOOLS_LIBRARIES} z pthread)
