# setting up cmake
cmake_minimum_required(VERSION 2.6)

# enabling compiler suite on osx
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    #set(CMAKE_C_COMPILER "gcc-mp-4.5")
    #set(CMAKE_CXX_COMPILER "gcc-mp-4.5")
endif()

if(NOT WIN32)
add_definitions(-std=c++0x)
endif()

# main project name defined
project(safegi_untyped)


# setting up macros and library paths
include(macros.cmake)
include(libs.cmake)

# set default build to release
default_release()

# setting up project
output_directories_default()
compiler_flags_default()
proj_directories_default()

# set renderer sources


add_subdirectory(libs/usafegi/)
add_subdirectory(libs/gi_aux)

add_subdirectory(apps/usafegi_gl)
add_subdirectory(apps/usafegi_ray)
