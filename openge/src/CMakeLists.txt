add_subdirectory( commands )
add_subdirectory( algorithms )
add_subdirectory( util )

set( OPENGE_LIB_SOURCES ${OPENGE_ALGORITHMS_SOURCES} ${OPENGE_COMMANDS_SOURCES} ${OPENGE_UTIL_SOURCES})
set( OPENGE_SOURCES main.cpp )

add_library(openge-lib ${OPENGE_LIB_SOURCES})
target_link_libraries(openge-lib ${Boost_LIBRARIES} ${BAMTOOLS_LIBRARIES} z pthread)

add_executable(openge ${OPENGE_SOURCES})
target_link_libraries(openge openge-lib)
