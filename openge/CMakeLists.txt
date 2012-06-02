
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost 1.45.0 REQUIRED COMPONENTS program_options)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

if( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR )
  # stand-alone build
  find_package(BamTools REQUIRED)
else()
  cmake_policy(SET CMP0015 NEW)
  include_directories(../include)
  get_property(BAMTOOLS_LIBRARIES TARGET BamTools-static PROPERTY LOCATION)
endif()


set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )
set( LIBRARY_OUTPUT_PATH    "${CMAKE_SOURCE_DIR}/lib" )

include_directories(${Boost_INCLUDE_DIRS} ${BAMTOOLS_INCLUDE_DIR})

add_subdirectory(src)
