# setting cmake configuration
macro(verbose_on)
    set(CMAKE_VERBOSE_MAKEFILE on)
endmacro()

# setting up output paths
macro(output_directories_default)
    if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
        set(BIN_DIR_SUFFIX mk)
    endif()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/../bin/${BIN_DIR_SUFFIX})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/../bin/${BIN_DIR_SUFFIX})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/../bin/${BIN_DIR_SUFFIX})
endmacro()

macro(compiler_flags_default)
    # setting up os specific paths
    if(APPLE)
        add_definitions(-Wall)
        add_definitions(-Wno-non-virtual-dtor)
        add_definitions(-Wno-unknown-pragmas)
        add_definitions(-std=c++0x)
    endif()
    if(WIN32)
        add_definitions(-D_SCL_SECURE_NO_WARNINGS)
        add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    endif()
endmacro()

# macro to add other libs directory
macro(proj_directories dir)
	#set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${DEP_DIR}/lib)
	#set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${DEP_DIR}/${CMAKE_CFG_INTDIR}/lib)
endmacro()

# macro to add default libs directory
macro(proj_directories_default)
    proj_directories(${PROJECT_SOURCE_DIR})
    foreach(dir ${ARGV})
        proj_directories(${dir})
    endforeach()
endmacro()

# macro that sets the default build type to release
macro(default_release)
    set(CMAKE_BUILD_TYPE Release)
endmacro()

macro(default_debug)
    set(CMAKE_BUILD_TYPE Debug)
endmacro()
