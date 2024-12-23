# setting up os specific paths
if(UNIX)
    find_library(LIBS_glew libGLEW.a ../dep/linux/lib)
    find_library(LIBS_glut libglut.a ../dep/linux/lib)
endif()
if(WIN32)
    find_library(LIBS_glew glew32 ../dep/w32/lib)
    find_library(LIBS_glut freeglut ../dep/w32/lib)
endif()

