# - Try to find bamtools
# Once done this will define
#  BAMTOOLS_FOUND - System has BAMTOOLS
#  BAMTOOLS_INCLUDE_DIRS - The BAMTOOLS include directories
#  BAMTOOLS_LIBRARIES - The libraries needed to use BAMTOOLS
#  BAMTOOLS_DEFINITIONS - Compiler switches required for using BAMTOOLS

find_package(PkgConfig)
pkg_check_modules(BAMTOOLS QUIET BamTools)
set(BAMTOOLS_DEFINITIONS ${PC_BAMTOOLS_CFLAGS_OTHER})

find_path(BAMTOOLS_INCLUDE_DIR api/BamReader.h
          HINTS ../bamtools ${PC_BAMTOOLS_INCLUDEDIR} ${PC_BAMTOOLS_INCLUDE_DIRS}
          PATH_SUFFIXES include)

find_library(BAMTOOLS_LIBRARY NAMES bamtools
             HINTS ../bamtools ${PC_BAMTOOLS_LIBDIR} ${PC_BAMTOOLS_LIBRARY_DIRS}
             PATH_SUFFIXES lib)

set(BAMTOOLS_LIBRARIES ${BAMTOOLS_LIBRARY} )
set(BAMTOOLS_INCLUDE_DIRS ${BAMTOOLS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set BAMTOOLS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(bamtools  DEFAULT_MSG
                                  BAMTOOLS_LIBRARY BAMTOOLS_INCLUDE_DIR)

mark_as_advanced(BAMTOOLS_INCLUDE_DIR BAMTOOLS_LIBRARY )
