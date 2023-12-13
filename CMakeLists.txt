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
        src/gaf.cpp
        src/Graph.cpp
        src/HeteroGraph.cpp
        src/IntervalGraph.cpp
        src/misc.cpp
        src/pair_hash.cpp
        src/Region.cpp
        src/Sequence.cpp
        src/TransitiveMap.cpp
        src/Timer.cpp
        inc/Sequence.hpp
        src/VcfReader.cpp
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
# ----------------- curl ------------------- #
# ------------------------------------------ #
##############################################

find_package( ZLIB REQUIRED )
find_package( LibLZMA REQUIRED )


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
        test_transmap
        test_vcfreader
        )

foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)

    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        target_link_libraries(${FILENAME_PREFIX}
                sv_merge
                Threads::Threads
                htslib
                ${ZLIB_LIBRARIES}
                curl2
                lzma
                bz2
                ${OPENSSL_LIBRARIES}
        )
    else()
        target_link_libraries(${FILENAME_PREFIX}
                sv_merge
                Threads::Threads
                htslib
                ${ZLIB_LIBRARIES}
                curl
                lzma
                bz2
                ${OPENSSL_LIBRARIES}
        )
    endif()

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        a
)

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        target_link_libraries(${FILENAME_PREFIX}
                sv_merge
                Threads::Threads
                htslib
                ${ZLIB_LIBRARIES}
                curl2
                ${OPENSSL_LIBRARIES}
                )
    else()
        target_link_libraries(${FILENAME_PREFIX}
                sv_merge
                Threads::Threads
                htslib
                ${ZLIB_LIBRARIES}
                curl
                ${OPENSSL_LIBRARIES}
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

