cmake_minimum_required(VERSION 3.27)

project(Ringli)

find_program(FFMPEG ffmpeg)
if ("${FFMPEG}" STREQUAL "FFMPEG-NOTFOUND")
    message(FATAL_ERROR "ffmpeg not found!")
endif()

find_program(XXD xxd)
if ("${XXD}" STREQUAL "XXD-NOTFOUND")
    message(FATAL_ERROR "xxd not found!")
endif()

include(FetchContent)

include(cmake/googletest.cmake)
include(cmake/highway.cmake)
include(cmake/eigen.cmake)
include(cmake/protobuf.cmake)
include(cmake/libsvm.cmake)
include(cmake/pffft.cmake)
include(cmake/visqol.cmake)
include(cmake/armadillo.cmake)
include(cmake/opus.cmake)
include(cmake/flac.cmake)
include(cmake/samplerate.cmake)
include(cmake/zimtohrli.cmake)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

enable_testing()
include(GoogleTest)

set(CMAKE_CXX_CLANG_TIDY clang-tidy)

add_subdirectory(common)
add_subdirectory(encode)
add_subdirectory(decode)
add_subdirectory(analysis)
