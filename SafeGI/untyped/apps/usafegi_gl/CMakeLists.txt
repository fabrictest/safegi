# AUX_SOURCE_DIRECTORY(. SOURCES)

SET(SOURCES
main.cpp
shader_src.h
shader_src.cpp
gl_renderer.h
gl_renderer.cpp
)

if(NOT WIN32)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lGL -lGLU -lXmu -lX11 -lm")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lGL -lGLU -lXmu -lX11 -lm")
endif()

ADD_EXECUTABLE(usafegi_gl ${SOURCES})
include_directories(../../libs/)


TARGET_LINK_LIBRARIES(usafegi_gl gi_aux ${LIBS_glew} ${LIBS_glut})
