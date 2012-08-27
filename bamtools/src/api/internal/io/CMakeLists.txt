# ==========================
# BamTools CMakeLists.txt
# (c) 2011 Derek Barnett
#
# src/api/internal/io
# ==========================

set ( InternalIODir "${InternalDir}/io" )

#--------------------------
# platform-independent IO
#--------------------------
set ( CommonIOSources
        ${InternalIODir}/BamDeviceFactory_p.cpp
        ${InternalIODir}/BamFile_p.cpp
        ${InternalIODir}/BamPipe_p.cpp
        ${InternalIODir}/BgzfStream_p.cpp
        ${InternalIODir}/ByteArray_p.cpp
        ${InternalIODir}/ILocalIODevice_p.cpp
        ${InternalIODir}/RollingBuffer_p.cpp
)

#---------------------------
# make build-specific list
#---------------------------
set ( InternalIOSources 
        ${CommonIOSources} 

        PARENT_SCOPE # <-- leave this last
)

