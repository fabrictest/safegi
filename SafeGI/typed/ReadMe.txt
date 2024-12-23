In this folder is the source code of typed SafeGI CPU renderer and GPU renderer prototype.

Compilation Instruction
-----------------------------------------------------------
For compilation instruction, please refer to "../ReadMe.txt"

Folder Content
-----------------------------------------------------------
ReadMe.txt          --This file
[ libs ]            **  Library source code folder.
[ libs\tsafegi ]    **  typed SafeGI library
    dimensions.h    --  contains the dimension definition and dimension analysis 
                        implementation mentioned in section 3.1.1 Dimensional 
                        Analysis.
    real.h          --  basic real number and tuple type without dimension.
    m_real.h        --  real number type and tuple types with dimensional 
                        constraint.
    m_units.h       --  units definition, such as meters, seconds, degrees..
    m_spectrum.h    --  spectrum type with dimensional analysis
    spaces.h        --  basic space definition, such as shape_s, world_s, local_s, 
                        mentioned in section 3.1.1 Space Analysis
    m_geom.h        --  geometry type with space analysis implementation mentioned
                        in section 3.2 Geometric Space Analysis.
                        + point     contains x,y,z components stores Cartesian 
                                    coordinates in length dimension.
                        + vector    contains x,y,z components in length dimension.
                        + direction directional vector, which is normalized to 
                                    length 1. Contains x,y,z components in unit 
                                    dimension.
                        + normal    specific direction type represent surface normal
                        + rigidmap  rigid mapping type used to transform point, vector,
                                    direction and normal from one space to another.
    xform.h         --  transform class definition
    linalg.h        --  basic linear algebra operation definitions.
    lens.h          --  lens class used to generate ray samples.
    m_sampling.h    --  different sample types, such as shadow sample, lens sample
                        direction sample and brdf sample. 
    direct.h        --  direct lighting integrator. 
                        + spectrum<radiance_d> l(..)
                            compute radiance spectrum for lens sample point p 
                            and direction w
                        + spectrum<radiance_d> le(..)
                            compute emission radiance spectrum of the surface.
                        + spectrum<radiance_d> ld(..)
                            compute direct illumination radiance spectrum of the 
                            surface.
    path.h          --  path tracer integrator.
                        +spectrum<radiance_d> li(..)
                            compute indirect illumination radiance spectrum of the 
                            surface
    brdf.h          --  brdf for materials, including Lambert, Phong and mirror. 
                        The brdf class is defined similar to BxDF class in PBRT. 
                        +spectrum<brdf_d> f(..)
                            return the brdf spectrum according of the brdf.
                        +spectrum<rho_d> rho(..)
                            return the rho spectrum of the brdf
                        +brdfSample<local_s> sample(..)
                            return the brdf sample of the brdf. brdfSample<object> 
                            contains the value of incoming direction sample, pdf 
                            and sampled brdf spectrum.
                        +mreal<invsolidangle_d> pdf(..)
                            return the pdf.
    scene.h         --  definitions for scene objects, such camera, light and 
                        materials.
    source.h        --  point light source and area light source.
    shape.h         --  base shapes, including sphere, quad and triangle.
    mesh.h          --  triangle mesh shape.
    bbox.h          --  bounding box class.
    tracer.h        --  ray tracer class
    tracerutils.h   --  utility class for ray tracing.
    type_trait.h    --  type trait deduction facilities used to impose type 
                        constraint between CPU-end and GPU-end of the SafeGI GPU 
                        renderer.
    safe_gl.h       --  typed OpenGL wrapper functions used by SafeGI GPU renderer.
    common.h        --  common header
    image.h         --  image class for rendering outputs.    
    std.h           --  standard C++ STL header
    stub.cpp        --  an empty cpp file
    

[ libs/gi_aux ]     --  SafeGI auxiliary classes.
    bunny.cpp       --  bunny data
    fileio_t.h      --  image I/O functions
    test_scene_t.h  --  test scene generation utilities
    timer.h         --  timer class.
    
[ libs/sparser ]        --  SafeGI GLSL parser for type checking.
    aux_data.cpp        --  parser auxiliary data structures.
    aux_data.h          --  parser auxiliary data structures.
    context.cpp         --  parsing context
    context.h           --  parsing context
    parser.y            --  bison file for parser.
    predef_symbols.h    --  pre-define file for dimensions, spaces and built-in 
                            functions mentioned in section 3.1.2 in the paper.
    scanner.l           --  flex file for scanner.
    type_checker.cpp    --  SafeGI GLSL type check interface
    type_checker.h      --  type check interface
    unistd.h

[ apps ]                **  Application source code folder.
[ apps/tsafegi_gl ]     **  Source code folder for typed SafeGI GPU renderer.
    gl_renderer.cpp     --  Source file for GPU renderer class. This class demonstrate how 
                            our SafeGI interface is used to create type safe GPU rendering 
                            program. 
    gl_renderer.h       --  GPU renderer class header
    shader_src.cpp      --  Typed GLSL source code for point light, area light and shadow mapping.
    shader_src.h        --  Typed GLSL source code header
    main.cpp            --  Main entry point of the typed SafeGI GPU renderer the program.
        
[ apps/tsafegi_ray ]    **  Source code folder for SafeGI CPU renderer.
    main.cpp            --  Main entry point of the typed SafeGI CPU ray tracer program.
