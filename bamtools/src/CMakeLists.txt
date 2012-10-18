# ==========================
# BamTools CMakeLists.txt
# (c) 2010 Derek Barnett
#
# src/
# ==========================

add_subdirectory( api )

# export shared headers
include( ExportHeader.cmake )
set( SharedIncludeDir "shared" )
ExportHeader( SharedHeaders shared/bamtools_global.h ${SharedIncludeDir} )
