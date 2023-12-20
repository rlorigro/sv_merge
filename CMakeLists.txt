######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.22)
project(sv_merge VERSION 0.0.0)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_SYSTEM: ${CMAKE_SYSTEM}")

# Compiler options.
set(OPENSSL_ROOT_DIR "/usr/local/opt/openssl/bin/openssl")
set(DOWNLOAD_EXTRACT_TIMESTAMP 1)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Commented until Ubuntu24/CMake3.24 becomes mainstream
#cmake_policy(SET CMP0135 NEW)
#set(CMAKE_POLICY_DEFAULT_CMP0135 NEW)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(LIB_SUFFIX "dylib")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LIB_SUFFIX "so")
endif()

# If using gperftools, run cmake with additional args:
# cmake -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg ..

if (debug)
    message(STATUS "---- Building with debug options ----")

    if (NOT no_asan OR NOT DEFINED no_asan)
        # Memory debug
        set(CMAKE_CXX_FLAGS "-fexceptions -fsanitize=address")
        set(ASAN_OPTIONS=check_initialization_order=1)
        set(ASAN_OPTIONS=detect_leaks=1)
    else()
        set(CMAKE_CXX_FLAGS "-fexceptions")
    endif()

    add_definitions(-ggdb3 -Og -Wall)       # Debugging + minimal optimization

else()
    message(STATUS "---- Building with optimization ----")

    # Standard compilation
    add_definitions(-O3 -Wall)              # Much optimization
endif()



#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test/executables.
set(SOURCES
        src/Alignment.cpp
        src/Authenticator.cpp
        src/bam.cpp
        src/bed.cpp
        src/fasta.cpp
        src/fetch.cpp
        src/gaf.cpp
        src/Graph.cpp
        src/HeteroGraph.cpp
        src/IntervalGraph.cpp
        src/misc.cpp
        src/pair_hash.cpp
        src/path_optimizer.cpp
        src/read_centroid_optimizer.cpp
        src/read_cluster_optimizer.cpp
        src/Region.cpp
        src/Sequence.cpp
        src/TransitiveMap.cpp
        src/Timer.cpp
        src/VcfReader.cpp
        src/windows.cpp
        )


project(sv_merge)
add_library(sv_merge STATIC ${SOURCES})

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(sv_merge PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(sv_merge PROPERTIES DEFINE_SYMBOL "")

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ------------------ WFA2 ------------------ #
# ------------------------------------------ #
##############################################


ExternalProject_Add(project_wfa
        GIT_REPOSITORY https://github.com/rlorigro/WFA2-lib.git
        GIT_TAG 4d0b9238857554b958ffa20048b155d29be92e04
#        URL https://github.com/smarco/WFA2-lib/archive/refs/tags/v2.3.4.tar.gz
        TIMEOUT 5
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/wfa2 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_RPATH=${CMAKE_SOURCE_DIR}/external/wfa2/lib/
        INSTALL_COMMAND make install
)

add_library(wfa2 SHARED IMPORTED)
set_property(TARGET wfa2
        PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/external/wfa2/lib/libwfa2cpp.${LIB_SUFFIX})
#                                            /sv_merge/external/wfa2/lib

add_dependencies(wfa2 project_wfa)
add_dependencies(sv_merge wfa2)

# Specify include dir
include_directories(
        ${CMAKE_SOURCE_DIR}/external/wfa2/include/wfa2lib
)

##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ----------------- curl ------------------- #
# ------------------------------------------ #
##############################################

find_package( ZLIB REQUIRED )


if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    message(STATUS "Fetching curl")

#    # This is not actually linked but it is required for Ubuntu Curl compilation to use http2
    find_library(NGHTTP2 REQUIRED NAMES nghttp2 libnghttp2 libnghttp2-dev)
#    find_package(NGHTTP2 REQUIRED)

#    message(STATUS "------------ found nghttp2 ${NGHTTP2_LIBRARIES}")


    include(FetchContent)

    FetchContent_Declare(
            project_curl
            URL      https://github.com/curl/curl/releases/download/curl-8_4_0/curl-8.4.0.tar.gz
    )

    FetchContent_MakeAvailable(project_curl)

    # Ensure that main library has access to primary dependencies' and secondary dependencies' headers
    include_directories(
            ${CMAKE_BINARY_DIR}/_deps/project_curl-src
            ${CMAKE_BINARY_DIR}/_deps/project_curl-src/curl/
    )

    message(STATUS "curl_SOURCE_DIR: ${project_curl_SOURCE_DIR}")
    message(STATUS "CMAKE_SOURCE_DIR/external/curl/: ${CMAKE_SOURCE_DIR}/external/curl/")


    #set(CURL_DIR _deps/project_curl-src)

    # Do not rebuild if -Ddev=1 is specified in cmake args
    if (NOT DEFINED dev OR NOT dev)
        # curl configure
        add_custom_target(
                BUILD_CURL
                ALL
                WORKING_DIRECTORY ${project_curl_SOURCE_DIR}
                COMMAND pwd
                COMMAND ./configure --with-openssl --with-nghttp2 --disable-ldap --disable-manual --disable-ldaps --disable-dict --disable-telnet --disable-tftp --disable-imap --disable-smb --disable-smtp --disable-gopher --disable-sspi --disable-ntlm-wb --disable-tls-srp --without-libidn2 --without-brotli --prefix=${CMAKE_SOURCE_DIR}/external/curl/
                COMMAND $(MAKE)
                COMMAND $(MAKE) install
        )

    endif()

    add_library(curl2 SHARED IMPORTED)
    set_property(TARGET curl2
            PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/external/curl/lib/libcurl.${LIB_SUFFIX})

    add_dependencies(curl2 BUILD_CURL zlib)
    add_dependencies(sv_merge curl2)
endif()


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- htslib ------------------ #
# ------------------------------------------ #
##############################################

message(STATUS "Fetching htslib")

find_library(BZ2LIB bz2)
if(${BZ2LIB} STREQUAL "BZ2LIB-NOTFOUND")
    message(WARNING "Couldn't find the 'bz2' library")
endif()

find_package( LibLZMA REQUIRED )

include(FetchContent)

FetchContent_Declare(
        project_htslib
        URL      https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2
)

FetchContent_MakeAvailable(project_htslib)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/htslib/
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/cram/
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/os/
        ${CMAKE_SOURCE_DIR}/external/
)

message(STATUS "htslib_SOURCE_DIR: ${project_htslib_SOURCE_DIR}")


set(HTS_DIR _deps/project_htslib-src)

# Do not rebuild if -Ddev=1 is specified in cmake args
if (NOT DEFINED dev OR NOT dev)
    # HTSLIB configure
    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set(CURL_FLAGS_1 "CPPFLAGS=-I${CMAKE_SOURCE_DIR}/external/curl/include/curl")
        set(CURL_FLAGS_2 "LDFLAGS='-L${CMAKE_SOURCE_DIR}/external/curl/lib'")
    else()
        set(CURL_FLAGS_1 "")
        set(CURL_FLAGS_2 "")
    endif()

    add_custom_target(
            BUILD_HTS
            ALL
            WORKING_DIRECTORY ${HTS_DIR}
            COMMAND pwd
            COMMAND autoreconf -i
            COMMAND ./configure --disable-lzma --without-libdeflate ${CURL_FLAGS_1} ${CURL_FLAGS_2}
            COMMAND $(MAKE) print-config
            COMMAND $(MAKE) prefix=${CMAKE_SOURCE_DIR}/external/htslib/ install
    )

endif()

add_library(htslib STATIC IMPORTED)
set_property(TARGET htslib
        PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/external/htslib/lib/libhts.a)


if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_dependencies(htslib BUILD_HTS curl2 zlib)
else()
    add_dependencies(htslib BUILD_HTS zlib)
endif()
add_dependencies(sv_merge htslib)


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# --------------- minimap2 ----------------- #
# ------------------------------------------ #
##############################################
#
#
#message(STATUS "Fetching minimap2")
#
#
## https://github.com/lh3/minimap2/archive/refs/tags/v2.24.tar.gz
#
#include(FetchContent)
#
#FetchContent_Declare(
#        project_minimap2
#        GIT_REPOSITORY https://github.com/lh3/minimap2.git
#        GIT_TAG v2.24
#)
#
#FetchContent_MakeAvailable(project_minimap2)
#
#set(MINIMAP2_SOURCES
#        ${project_minimap2_SOURCE_DIR}/ksw2_dispatch.c
#        ${project_minimap2_SOURCE_DIR}/lchain.c
#        ${project_minimap2_SOURCE_DIR}/format.c
##        ${project_minimap2_SOURCE_DIR}/example.c
#        ${project_minimap2_SOURCE_DIR}/ksw2_extz2_sse.c
#        ${project_minimap2_SOURCE_DIR}/ksw2_extd2_sse.c
#        ${project_minimap2_SOURCE_DIR}/kalloc.c
#        ${project_minimap2_SOURCE_DIR}/kthread.c
#        ${project_minimap2_SOURCE_DIR}/hit.c
#        ${project_minimap2_SOURCE_DIR}/ksw2_ll_sse.c
#        ${project_minimap2_SOURCE_DIR}/index.c
#        ${project_minimap2_SOURCE_DIR}/ksw2_exts2_sse.c
#        ${project_minimap2_SOURCE_DIR}/align.c
#        ${project_minimap2_SOURCE_DIR}/esterr.c
#        ${project_minimap2_SOURCE_DIR}/bseq.c
#        ${project_minimap2_SOURCE_DIR}/main.c
#        ${project_minimap2_SOURCE_DIR}/map.c
#        ${project_minimap2_SOURCE_DIR}/options.c
#        ${project_minimap2_SOURCE_DIR}/pe.c
#        ${project_minimap2_SOURCE_DIR}/sdust.c
#        ${project_minimap2_SOURCE_DIR}/sketch.c
#        ${project_minimap2_SOURCE_DIR}/splitidx.c
#        ${project_minimap2_SOURCE_DIR}/seed.c
#        ${project_minimap2_SOURCE_DIR}/misc.c
#        )
#
#message(STATUS "minimap2_SOURCE_DIR: ${project_minimap2_SOURCE_DIR}")
#
#
#set(MINIMAP2_DIR _deps/project_minimap2-src)
#
#add_custom_target(
#        BUILD_MINIMAP2
#        ALL
#        WORKING_DIRECTORY ${MINIMAP2_DIR}
#        COMMAND ""
#)
#
## Ensure that main library has access to primary dependencies' and secondary dependencies' headers
#include_directories(
#        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/
#        #        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/sse2neon/
#        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/
#        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/simde/
#        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/simde/simde/
#        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/simde/simde/arm/
#        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/simde/simde/arm/neon/
#        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/simde/simde/x86/
#)
#
##add_custom_target(HTSLIB_CONFIGURED DEPENDS "${CMAKE_BINARY_DIR}/_deps/project_htslib-src/config.h")
#
#add_library(minimap2 ${MINIMAP2_SOURCES})
#add_dependencies(minimap2 BUILD_MINIMAP2)
#add_dependencies(sv_merge minimap2)
#
#target_link_libraries(minimap2 pthread bz2 z)


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- ortools ----------------- #
# ------------------------------------------ #
##############################################

include(FetchContent)
FetchContent_Declare(
        ortools
        URL https://github.com/google/or-tools/archive/refs/tags/v9.7.tar.gz
)

set(BUILD_DEPS ON)
set(BUILD_DOC OFF)
set(BUILD_SAMPLES OFF)
set(BUILD_EXAMPLES OFF)
set(USE_COINOR OFF)
FetchContent_MakeAvailable(ortools)

include_directories(
        ${ortools_SOURCE_DIR}
        ${ortools_BINARY_DIR}
        ${absl_SOURCE_DIR}
        ${protobuf_SOURCE_DIR}/src
)

message(STATUS "ortools_SOURCE_DIR: " ${ortools_SOURCE_DIR})
message(STATUS "ortools_BINARY_DIR: " ${ortools_BINARY_DIR})
message(STATUS "absl_SOURCE_DIR: " ${absl_SOURCE_DIR})
message(STATUS "protobuf_SOURCE_DIR: /src: " ${protobuf_SOURCE_DIR}/src)

add_dependencies(sv_merge ortools)

##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ------------------ bdsg ------------------ #
# ------------------------------------------ #
##############################################

message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

if (dev)
    # Download or update library as an external project
    ExternalProject_Add(project_bdsg
            GIT_REPOSITORY https://github.com/vgteam/libbdsg.git
            DOWNLOAD_COMMAND ""
            UPDATE_COMMAND ""
            PREFIX ${CMAKE_SOURCE_DIR}/external/bdsg/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bdsg/ -DRUN_DOXYGEN=OFF -DBUILD_PYTHON_BINDINGS=OFF
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/bdsg/
            INSTALL_COMMAND make install
    )
else()
    # Download or update library as an external project
    ExternalProject_Add(project_bdsg
            GIT_REPOSITORY https://github.com/vgteam/libbdsg.git
            PREFIX ${CMAKE_SOURCE_DIR}/external/bdsg/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bdsg/ -DRUN_DOXYGEN=OFF -DBUILD_PYTHON_BINDINGS=OFF
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/bdsg/
            INSTALL_COMMAND make install
    )
endif()

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_bdsg INSTALL_DIR)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(LIB_SUFFIX "dylib")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LIB_SUFFIX "so")
endif()


# Create new library for external project (so it can be linked with main library)
add_library(bdsg SHARED IMPORTED)
set_property(TARGET bdsg
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libbdsg.${LIB_SUFFIX})

add_library(divsufsort SHARED IMPORTED)
set_property(TARGET divsufsort
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libdivsufsort64.${LIB_SUFFIX})

add_library(libhandlegraph SHARED IMPORTED)
set_property(TARGET libhandlegraph
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libhandlegraph.${LIB_SUFFIX})

add_library(libsdsl SHARED IMPORTED)
set_property(TARGET libsdsl
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libsdsl.${LIB_SUFFIX})

# Define library as dependent on the downloaded project
add_dependencies(bdsg
        project_bdsg
        libsdsl
        libhandlegraph
        divsufsort)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(sv_merge bdsg)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/bdsg/include/
        external/bdsg/include/handlegraph/
        external/bdsg/include/bdsg/)

# Add runtime path for main library so it can access external library
#set_property(TARGET GFAse PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/bdsg")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")

# ------- Misc -------

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_package(OpenSSL REQUIRED)

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    message(STATUS "using compiled curl")
else()
    find_package(CURL REQUIRED)
endif()


#message(STATUS "CURL_LIBRARIES: ${CURL_LIBRARIES}")


# -------- TESTS --------

set(TESTS
        test_bam
        test_bed
        test_fasta
        test_gaf
        test_graph
        test_graphaligner
        test_heterograph
        test_htslib
        test_hts_gcp
        test_hts_gcp_speed
        test_interval_graph
#        test_minimap2
#        test_minimap2_no_io
        test_ortools
        test_read_clustering
        test_transmap
        test_triangle_inequality
        test_wfa2
        test_vcfreader
        )

foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)

    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        target_link_libraries(${FILENAME_PREFIX}
                sv_merge
                Threads::Threads
                htslib
                bdsg
                ${ZLIB_LIBRARIES}
                curl2
                ${OPENSSL_LIBRARIES}
                bz2
                lzma
#                minimap2
                ortools::ortools
                wfa2
        )
    else()
        target_link_libraries(${FILENAME_PREFIX}
                sv_merge
                Threads::Threads
                htslib
                bdsg
                ${ZLIB_LIBRARIES}
                curl
                ${OPENSSL_LIBRARIES}
                bz2
                lzma
#                minimap2
                ortools::ortools
                wfa2
        )
    endif()

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        a
        hapestry
        find_windows
)

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        target_link_libraries(${FILENAME_PREFIX}
                sv_merge
                Threads::Threads
                htslib
                bdsg
                ${ZLIB_LIBRARIES}
                curl2
                ${OPENSSL_LIBRARIES}
                bz2
                lzma
#                minimap2
                ortools::ortools
                wfa2
                )
    else()
        target_link_libraries(${FILENAME_PREFIX}
                sv_merge
                Threads::Threads
                htslib
                bdsg
                ${ZLIB_LIBRARIES}
                curl
                ${OPENSSL_LIBRARIES}
                bz2
                lzma
#                minimap2
                ortools::ortools
                wfa2
        )
    endif()

endforeach()

#set_target_properties(sv_merge PROPERTIES LINK_FLAGS "-static" )
#SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")


# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

## The shared library goes to the bin directory.
#install(TARGETS
#        sv_merge
#        sv_merge
#        DESTINATION ${BINARY_INSTALL_DIR})

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

