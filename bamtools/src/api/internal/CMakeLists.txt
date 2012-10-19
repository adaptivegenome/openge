# ==========================
# BamTools CMakeLists.txt
# (c) 2011 Derek Barnett
#
# src/api/internal
# ==========================

set( InternalDir "internal" )

add_subdirectory( sam )
add_subdirectory( utils )

set( InternalSources
         ${InternalSamSources}
         ${InternalUtilsSources}

         PARENT_SCOPE # <-- leave this last
   )

