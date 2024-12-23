# setting up cmake
cmake_minimum_required(VERSION 2.6)

# enabling compiler suite on osx
#if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
#    set(CMAKE_C_COMPILER "/opt/local/bin/gcc-mp-4.5")
#    set(CMAKE_CXX_COMPILER "/opt/local/bin/gcc-mp-4.5")
#endif()

# main project name defined
project(safe_pbrt)

# setting up macros
include(macros.cmake)

# set default build to release
default_release()

# setting up project
output_directories_default()
compiler_flags_default()
proj_directories_default()

# include external libraries
#include(libs.cmake)

# libs subdir list
add_subdirectory(pbrt)
