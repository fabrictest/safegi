In this folder is the source code of untyped SafeGI CPU renderer and GPU renderer prototype.

Compilation Instruction
-----------------------------------------------------------
For compilation instruction, please refer to "../ReadMe.txt"

Folder Content
-----------------------------------------------------------
ReadMe.txt          --This file
[ libs ]            **  Library source code folder.
[ libs\usafegi ]    **  untyped SafeGI library
    real.h          --  basic real number and tuple type.
    u_units.h       --  units definition, such as meters, seconds, degrees..
    u_spectrum.h    --  spectrum type
    u_geom.h        --  geometry type
                        + point     contains x,y,z components stores Cartesian 
                                    coordinates.
                        + vector    contains x,y,z components.
                        + direction directional vector, which is normalized to 
                                    length 1.
                        + normal    specific direction type represent surface normal
                        + rigidmap  rigid mapping type used to transform point, vector,
                                    direction and normal.
    xform.h         --  transform class definition
    linalg.h        --  basic linear algebra operation definitions.
    lens.h          --  lens class used to generate ray samples.
    u_sampling.h    --  different sample types, such as shadow sample, lens sample
                        direction sample and brdf sample. 
    direct.h        --  direct lighting integrator. 
                        + spectrum l(..)
                            compute radiance spectrum for lens sample point p 
                            and direction w
                        + spectrum le(..)
                            compute emission radiance spectrum of the surface.
                        + spectrum ld(..)
                            compute direct illumination radiance spectrum of the 
                            surface.
    path.h          --  path tracer integrator.
                        +spectrum li(..)
                            compute indirect illumination radiance spectrum of the 
                            surface
    brdf.h          --  brdf for materials, including Lambert, Phong and mirror. 
                        The brdf class is defined similar to BxDF class in PBRT. 
                        +spectrum f(..)
                            return the brdf spectrum according of the brdf.
                        +spectrum rho(..)
                            return the rho spectrum of the brdf
                        +brdfSample sample(..)
                            return the brdf sample of the brdf. brdfSample<object> 
                            contains the value of incoming direction sample, pdf 
                            and sampled brdf spectrum.
                        +real pdf(..)
                            return the pdf.
    scene.h         --  definitions for scene objects, such camera, light and 
                        materials.
    source.h        --  point light source and area light source.
    shape.h         --  base shapes, including sphere, quad and triangle.
    mesh.h          --  triangle mesh shape.
    bbox.h          --  bounding box class.
    tracer.h        --  ray tracer class
    tracerutils.h   --  utility class for ray tracing.
    safe_gl.h       --  untyped OpenGL wrapper functions used by SafeGI GPU renderer.
    common.h        --  common header
    image.h         --  image class for rendering outputs.    
    std.h           --  standard C++ STL header
    stub.cpp        --  an empty cpp file
    

[ libs/gi_aux ]     --  SafeGI auxiliary classes.
    bunny.cpp       --  bunny data
    fileio_u.h      --  image I/O functions
    test_scene_u.h  --  test scene generation utilities
    timer.h         --  timer class.
    
[ apps ]                **  Application source code folder.
[ apps/usafegi_gl ]     **  Source code folder for untyped SafeGI GPU renderer.
    gl_renderer.cpp     --  Source file for SafeGI GPU renderer class used as a reference 
                            implementation of untyped GPU renderer.
    gl_renderer.h       --  untyped GPU renderer class header
    shader_src.cpp      --  untyped GLSL source code for point light, area light and shadow mapping.
    shader_src.h        --  untyped GLSL source code header
    main.cpp            --  Main entry point of untyped GPU renderer program.
        
[ apps/tsafegi_ray ]    **  Source code folder for SafeGI CPU renderer.
    main.cpp            --  Main entry point of untyped ray tracer program.

