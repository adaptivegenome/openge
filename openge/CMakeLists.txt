
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost 1.45.0 REQUIRED COMPONENTS program_options)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )
set( LIBRARY_OUTPUT_PATH    "${CMAKE_SOURCE_DIR}/lib" )
set( DOCUMENTATION_OUTPUT_PATH    "${CMAKE_SOURCE_DIR}/doc" )

include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(src)
add_subdirectory(test)
#add_subdirectory(doc)
