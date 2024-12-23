
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


// core/api.cpp*
#include "api.h"
#include "parallel.h"
#include "paramset.h"
#include "spectrum.h"
#include "scene.h"
#include "renderer.h"
#include "film.h"
#include "volume.h"

// API Additional Headers
#include "accelerators/bvh.h"
#include "cameras/orthographic.h"
#include "cameras/perspective.h"
#include "film/image.h"
#include "filters/box.h"
#include "integrators/directlighting.h"
#include "integrators/emission.h"
#include "integrators/path.h"
#include "integrators/photonmap.h"
#include "lights/diffuse.h"
#include "lights/infinite.h"
#include "lights/point.h"
#include "materials/glass.h"
#include "materials/matte.h"
#include "materials/metal.h"
#include "materials/mirror.h"
#include "materials/plastic.h"
#include "materials/uber.h"
#include "renderers/samplerrenderer.h"
#include "samplers/lowdiscrepancy.h"
#include "samplers/random.h"
#include "shapes/disk.h"
#include "shapes/sphere.h"
#include "shapes/trianglemesh.h"
#include "textures/constant.h"
#include "textures/imagemap.h"
#include "textures/scale.h"
#include "volumes/homogeneous.h"

#include <map>
 #if (_MSC_VER >= 1400)
 #include <stdio.h>
 #define snprintf _snprintf
 #endif
using std::map;

// API Global Variables
Options PbrtOptions;

// API Local Classes
#define MAX_TRANSFORMS 2
#define START_TRANSFORM_BITS (1 << 0)
#define END_TRANSFORM_BITS   (1 << 1)
#define ALL_TRANSFORMS_BITS ((1 << MAX_TRANSFORMS) - 1)

struct TransformSet {
   // TransformSet Public Methods
   template<typename S1, typename S2>
   Transform<S1, S2> &get(int i)
   {
       Assert(i >= 0 && i < MAX_TRANSFORMS);
       return *reinterpret_cast<Transform<S1, S2>*>(&t[i]);
   }

   template<typename S1, typename S2>
   const Transform<S1, S2> &get(int i) const
   {
       Assert(i >= 0 && i < MAX_TRANSFORMS);
       return *reinterpret_cast<const Transform<S1, S2>*>(&t[i]);
   }

   friend TransformSet Inverse(const TransformSet &ts) {
       TransformSet t2;
       for (int i = 0; i < MAX_TRANSFORMS; ++i)
           t2.t[i] = Inverse(ts.t[i]);
       return t2;
   }
   bool IsAnimated() const {
       for (int i = 0; i < MAX_TRANSFORMS-1; ++i)
           if (t[i] != t[i+1]) return true;
       return false;
   }
private:
    Transform<unknown_s, unknown_s> t[MAX_TRANSFORMS];
};


struct RenderOptions {
    // RenderOptions Public Methods
    RenderOptions();
    Scene *MakeScene();
    Camera *MakeCamera() const;
    Renderer *MakeRenderer() const;

    // RenderOptions Public Data
    mfloat<time_d> transformStartTime, transformEndTime;
    string FilterName;
    ParamSet FilterParams;
    string FilmName;
    ParamSet FilmParams;
    string SamplerName;
    ParamSet SamplerParams;
    string AcceleratorName;
    ParamSet AcceleratorParams;
    string RendererName;
    string SurfIntegratorName, VolIntegratorName;
    ParamSet RendererParams;
    ParamSet SurfIntegratorParams, VolIntegratorParams;
    string CameraName;
    ParamSet CameraParams;
    TransformSet CameraToWorld;
    std::vector<Light *> lights;
    std::vector<Reference<Primitive> > primitives;
    mutable std::vector<VolumeRegion *> volumeRegions;
    map<string, std::vector<Reference<Primitive> > > instances;
    std::vector<Reference<Primitive> > *currentInstance;
};


RenderOptions::RenderOptions() {
    // RenderOptions Constructor Implementation
    transformStartTime = mfloat<time_d>(0.f);
    transformEndTime = mfloat<time_d>(1.f);
    FilterName = "box";
    FilmName = "image";
    SamplerName = "lowdiscrepancy";
    AcceleratorName = "bvh";
    RendererName = "sample";
    SurfIntegratorName = "directlighting";
    VolIntegratorName = "emission";
    CameraName = "perspective";
    currentInstance = NULL;
}


struct GraphicsState {
    // Graphics State Methods
    GraphicsState();
    Reference<Material> CreateMaterial(const ParamSet &params);

    // Graphics State
    map<string, Reference<Texture<mfloat<unit_d> > > > floatTextures;
    map<string, Reference<Texture<Spectrum<rho_d> > > > spectrumTextures;
    ParamSet materialParams;
    string material;
    map<string, Reference<Material> > namedMaterials;
    string currentNamedMaterial;
    ParamSet afloatightParams;
    string areaLight;
    bool reverseOrientation;
};


GraphicsState::GraphicsState() {
    // GraphicsState Constructor Implementation
    material = "matte";
    reverseOrientation = false;
}

typedef Transform<unknown_s, unknown_s> UnknownTransform;
class TransformCache {
public:
    // TransformCache Public Methods
    template<typename S1, typename S2>
    void Lookup(const Transform<S1, S2> &t, Transform<S1, S2> **tCached = NULL,
            Transform<S2, S1> **tCachedInverse = NULL) 
    {
        const Transform<unknown_s, unknown_s> &tt = 
            *reinterpret_cast<const UnknownTransform*>(&t);
        map<UnknownTransform, 
            std::pair<UnknownTransform *, UnknownTransform *> >::iterator iter;
        iter = cache.find(tt);
        if (iter == cache.end()) {
            Transform<S1, S2> *tr = arena.Alloc<Transform<S1, S2> >();
            *tr = t;
            Transform<S2, S1> *tinv = arena.Alloc<Transform<S2, S1>>();
            *tinv = Inverse(t);
            cache[tt] = 
                std::make_pair(reinterpret_cast<UnknownTransform*>(tr), reinterpret_cast<UnknownTransform*>(tinv));
            iter = cache.find(tt);
        }
        if (tCached) *tCached = 
            reinterpret_cast<Transform<S1, S2> *>(iter->second.first);
        if (tCachedInverse) *tCachedInverse = 
            reinterpret_cast<Transform<S2, S1> *>(iter->second.second);
        PBRT_ALLOCATED_CACHED_TRANSFORM();
    }
    void Clear() {
        cache.erase(cache.begin(), cache.end());
    }
private:
    // TransformCache Private Data
    map<UnknownTransform, std::pair<UnknownTransform *, UnknownTransform *> > cache;
    MemoryArena arena;
};



// API Static Data
#define STATE_UNINITIALIZED  0
#define STATE_OPTIONS_BLOCK  1
#define STATE_WORLD_BLOCK    2
static int currentApiState = STATE_UNINITIALIZED;
static TransformSet curTransform;
static int activeTransformBits = ALL_TRANSFORMS_BITS;
static map<string, TransformSet> namedCoordinateSystems;
static RenderOptions *renderOptions = NULL;
static GraphicsState graphicsState;
static std::vector<GraphicsState> pushedGraphicsStates;
static std::vector<TransformSet> pushedTransforms;
static std::vector<uint32_t> pushedActiveTransformBits;
static TransformCache transformCache;

// API Macros
#define VERIFY_INITIALIZED(func) \
if (currentApiState == STATE_UNINITIALIZED) { \
    Error("pbrtInit() must be before calling \"%s()\". " \
          "Ignoring.", func); \
    return; \
} else /* swallow trailing semicolon */
#define VERIFY_OPTIONS(func) \
VERIFY_INITIALIZED(func); \
if (currentApiState == STATE_WORLD_BLOCK) { \
    Error("Options cannot be set inside world block; " \
          "\"%s\" not allowed.  Ignoring.", func); \
    return; \
} else /* swallow trailing semicolon */
#define VERIFY_WORLD(func) \
VERIFY_INITIALIZED(func); \
if (currentApiState == STATE_OPTIONS_BLOCK) { \
    Error("Scene description must be inside world block; " \
          "\"%s\" not allowed. Ignoring.", func); \
    return; \
} else /* swallow trailing semicolon */
#define FOR_ACTIVE_TRANSFORMS \
    for (int i = 0; i < MAX_TRANSFORMS; ++i) \
        if (activeTransformBits & (1 << i))
#define WARN_IF_ANIMATED_TRANSFORM(func) \
do { if (curTransform.IsAnimated()) \
         Warning("Animated transformations set; ignoring for \"%s\"" \
                 "and using the start transform only", func); \
} while (false)

// Object Creation Function Definitions
Reference<Shape> MakeShape(const string &name,
        const Transform<object_s, world_s> *object2world, const Transform<world_s, object_s> *world2object,
        bool reverseOrientation, const ParamSet &paramSet) {
    Shape *s = NULL;
    if (name == "sphere")
        s = CreateSphereShape(object2world, world2object, reverseOrientation,
                              paramSet);
    // Create remaining \use{Shape} types
    else if (name == "disk")
        s = CreateDiskShape(object2world, world2object, reverseOrientation,
                            paramSet);
    else if (name == "trianglemesh")
        s = CreateTriangleMeshShape(object2world, world2object, reverseOrientation,
                                    paramSet);
    else
        Warning("Shape \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return s;
}


Reference<Material> MakeMaterial(const string &name,
        const Transform<material_s, world_s> &mtl2world,
        const TextureParams &mp) {
    Material *material = NULL;
    if (name == "matte")
        material = CreateMatteMaterial(mtl2world, mp);
    else if (name == "plastic")
        material = CreatePlasticMaterial(mtl2world, mp);
    else if (name == "glass")
        material = CreateGlassMaterial(mtl2world, mp);
    else if (name == "mirror")
        material = CreateMirrorMaterial(mtl2world, mp);
    else if (name == "metal")
        material = CreateMetalMaterial(mtl2world, mp);
    else if (name == "uber")
        material = CreateUberMaterial(mtl2world, mp);
    else
        Warning("Material \"%s\" unknown.", name.c_str());
    mp.ReportUnused();
    if (!material) Error("Unable to create material \"%s\"\n", name.c_str());
    return material;
}


Reference<Texture<mfloat<unit_d> > > MakeFloatTexture(const string &name,
        const Transform<texture_s, world_s> &tex2world, const TextureParams &tp) {
    Texture<mfloat<unit_d> > *tex = NULL;
    if (name == "constant")
        tex = CreateConstantFloatTexture(tex2world, tp);
    else if (name == "scale")
        tex = CreateScaleFloatTexture(tex2world, tp);
    else if (name == "imagemap")
        tex = CreateImageFloatTexture(tex2world, tp);
    else
        Warning("Float texture \"%s\" unknown.", name.c_str());
    tp.ReportUnused();
    return tex;
}


Reference<Texture<Spectrum<rho_d> > > MakeSpectrumTexture(const string &name,
        const Transform<texture_s, world_s> &tex2world, const TextureParams &tp) {
    Texture<Spectrum<rho_d> > *tex = NULL;
    if (name == "constant")
        tex = CreateConstantSpectrumTexture(tex2world, tp);
    else if (name == "scale")
        tex = CreateScaleSpectrumTexture(tex2world, tp);
    else if (name == "imagemap")
        tex = CreateImageSpectrumTexture(tex2world, tp);
    else
        Warning("Spectrum texture \"%s\" unknown.", name.c_str());
    tp.ReportUnused();
    return tex;
}


Light *MakeLight(const string &name,
        const Transform<light_s, world_s> &light2world, const ParamSet &paramSet) {
    Light *light = NULL;
    if (name == "point")
        light = CreatePointLight(light2world, paramSet);
    else if (name == "infinite" || name == "exinfinite")
        light = CreateInfiniteLight(light2world, paramSet);
    else
        Warning("Light \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return light;
}


AreaLight *MakeAreaLight(const string &name,
        const Transform<light_s, world_s> &light2world, const ParamSet &paramSet,
        const Reference<Shape> &shape) {
    AreaLight *area = NULL;
    if (name == "area" || name == "diffuse")
        area = CreateDiffuseAreaLight(light2world, paramSet, shape);
    else
        Warning("Area light \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return area;
}


VolumeRegion *MakeVolumeRegion(const string &name,
        const Transform<volume_s, world_s> &volume2world, const ParamSet &paramSet) {
    VolumeRegion *vr = NULL;
    if (name == "homogeneous")
        vr = CreateHomogeneousVolumeDensityRegion(volume2world, paramSet);
    else
        Warning("Volume region \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return vr;
}


SurfaceIntegrator *MakeSurfaceIntegrator(const string &name,
        const ParamSet &paramSet) {
    SurfaceIntegrator *si = NULL;
    if (name == "directlighting")
        si = CreateDirectLightingIntegrator(paramSet);
    else if (name == "path")
        si = CreatePathSurfaceIntegrator(paramSet);
    else if (name == "photonmap" || name == "exphotonmap")
        si = CreatePhotonMapSurfaceIntegrator(paramSet);
    else
        Warning("Surface integrator \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return si;
}


VolumeIntegrator *MakeVolumeIntegrator(const string &name,
        const ParamSet &paramSet) {
    VolumeIntegrator *vi = NULL;
    if (name == "null")
        vi = NULL;
    else if (name == "emission")
        vi = CreateEmissionVolumeIntegrator(paramSet);
    else
        Warning("Volume integrator \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return vi;
}


Primitive *MakeAccelerator(const string &name,
        const std::vector<Reference<Primitive> > &prims,
        const ParamSet &paramSet) {
    Primitive *accel = NULL;
    if (name == "bvh")
        accel = CreateBVHAccelerator(prims, paramSet);
    else
        Warning("Accelerator \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return accel;
}


Camera *MakeCamera(const string &name,
        const ParamSet &paramSet,
        const TransformSet &cam2worldSet, const mfloat<time_d> transformStart,
        const mfloat<time_d> transformEnd, Film *film) {
    Camera *camera = NULL;
    Assert(MAX_TRANSFORMS == 2);
    Transform<camera_s, world_s> *cam2world[2];
    transformCache.Lookup(cam2worldSet.get<camera_s, world_s>(0), &cam2world[0]);
    transformCache.Lookup(cam2worldSet.get<camera_s, world_s>(1), &cam2world[1]);
    AnimatedTransform<camera_s, world_s> animatedCam2World(cam2world[0], transformStart,
        cam2world[1], transformEnd);
    if (name == "perspective")
        camera = CreatePerspectiveCamera(paramSet, animatedCam2World, film);
    else if (name == "orthographic")
        camera = CreateOrthographicCamera(paramSet, animatedCam2World, film);
    else
        Warning("Camera \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return camera;
}


Sampler *MakeSampler(const string &name,
        const ParamSet &paramSet, const Film *film, const Camera *camera) {
    Sampler *sampler = NULL;
    if (name == "lowdiscrepancy")
        sampler = CreateLowDiscrepancySampler(paramSet, film, camera);
    else if (name == "random")
        sampler = CreateRandomSampler(paramSet, film, camera);
    else
        Warning("Sampler \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return sampler;
}


Filter *MakeFilter(const string &name,
    const ParamSet &paramSet) {
    Filter *filter = NULL;
    if (name == "box")
        filter = CreateBoxFilter(paramSet);
    else
        Warning("Filter \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return filter;
}


Film *MakeFilm(const string &name,
    const ParamSet &paramSet, Filter *filter) {
    Film *film = NULL;
    if (name == "image")
        film = CreateImageFilm(paramSet, filter);
    else
        Warning("Film \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return film;
}



// API Function Definitions
void pbrtInit(const Options &opt) {
    PbrtOptions = opt;
    // API Initialization
    if (currentApiState != STATE_UNINITIALIZED)
        Error("pbrtInit() has already been called.");
    currentApiState = STATE_OPTIONS_BLOCK;
    renderOptions = new RenderOptions;
    graphicsState = GraphicsState();
}


void pbrtCleanup() {
    ProbesCleanup();
    // API Cleanup
    if (currentApiState == STATE_UNINITIALIZED)
        Error("pbrtCleanup() called without pbrtInit().");
    else if (currentApiState == STATE_WORLD_BLOCK)
        Error("pbrtCleanup() called while inside world block.");
    currentApiState = STATE_UNINITIALIZED;
    delete renderOptions;
    renderOptions = NULL;
}


void pbrtIdentity() {
    VERIFY_INITIALIZED("Identity");
    FOR_ACTIVE_TRANSFORMS
    {
        curTransform.get<unknown_s, unknown_s>(i) = Transform<unknown_s, unknown_s>();
    }
}


void pbrtTranslate(float dx, float dy, float dz) {
    VERIFY_INITIALIZED("Translate");
    FOR_ACTIVE_TRANSFORMS
    {
        curTransform.get<unknown_s, unknown_s>(i) =
        curTransform.get<unknown_s, unknown_s>(i) * 
        Translate<unknown_s, unknown_s>(Vector<unknown_s>(dx*meters, dy*meters, dz*meters));
    }
}


void pbrtTransform(float tr[16]) {
    VERIFY_INITIALIZED("Transform");
    FOR_ACTIVE_TRANSFORMS
    {
        curTransform.get<unknown_s, unknown_s>(i) = Transform<unknown_s, unknown_s>(Matrix4x4(
        tr[0], tr[4], tr[8], tr[12],
        tr[1], tr[5], tr[9], tr[13],
        tr[2], tr[6], tr[10], tr[14],
        tr[3], tr[7], tr[11], tr[15]));
    }
}


void pbrtConcatTransform(float tr[16]) {
    VERIFY_INITIALIZED("ConcatTransform");
    FOR_ACTIVE_TRANSFORMS
    {
        curTransform.get<unknown_s, unknown_s>(i) = curTransform.get<unknown_s, unknown_s>(i) * Transform<unknown_s, unknown_s>(
                Matrix4x4(tr[0], tr[4], tr[8], tr[12],
                          tr[1], tr[5], tr[9], tr[13],
                          tr[2], tr[6], tr[10], tr[14],
                          tr[3], tr[7], tr[11], tr[15]));
    }
}


void pbrtRotate(float angle, float dx, float dy, float dz) {
    VERIFY_INITIALIZED("Rotate");
    FOR_ACTIVE_TRANSFORMS
    {
        curTransform.get<unknown_s, unknown_s>(i) = 
            curTransform.get<unknown_s, unknown_s>(i) * 
            Rotate<unknown_s, unknown_s>(angle, Direction<unknown_s>(dx, dy, dz));
    }
}


void pbrtScale(float sx, float sy, float sz) {
    VERIFY_INITIALIZED("Scale");
    FOR_ACTIVE_TRANSFORMS
    {
        curTransform.get<unknown_s, unknown_s>(i) = 
        curTransform.get<unknown_s, unknown_s>(i) * Scale<unknown_s, unknown_s>(sx, sy, sz);
    }
}


void pbrtLookAt(float ex, float ey, float ez, float lx, float ly,
        float lz, float ux, float uy, float uz) {
    VERIFY_INITIALIZED("LookAt");
    FOR_ACTIVE_TRANSFORMS{
        { Warning("This version of pbrt fixes a bug in the LookAt transformation.\n"
                                    "If your rendered images unexpectedly change, add a \"Scale -1 1 1\"\n"
                                    "to the start of your scene file."); break; }
    }
    FOR_ACTIVE_TRANSFORMS
    {
        curTransform.get<unknown_s, unknown_s>(i) =
        curTransform.get<unknown_s, unknown_s>(i) * LookAt<unknown_s, unknown_s>(Point<unknown_s>(ex*meters, ey*meters, ez*meters), 
        Point<unknown_s>(lx*meters, ly*meters, lz*meters), Direction<unknown_s>(ux, uy, uz));
    }
}


void pbrtCoordinateSystem(const string &name) {
    VERIFY_INITIALIZED("CoordinateSystem");
    namedCoordinateSystems[name] = curTransform;
}


void pbrtCoordSysTransform(const string &name) {
    VERIFY_INITIALIZED("CoordSysTransform");
    if (namedCoordinateSystems.find(name) !=
        namedCoordinateSystems.end())
        curTransform = namedCoordinateSystems[name];
    else
        Warning("Couldn't find named coordinate system \"%s\"",
                name.c_str());
}


void pbrtActiveTransformAll() {
    activeTransformBits = ALL_TRANSFORMS_BITS;
}


void pbrtActiveTransformEndTime() {
    activeTransformBits = END_TRANSFORM_BITS;
}


void pbrtActiveTransformStartTime() {
    activeTransformBits = START_TRANSFORM_BITS;
}


void pbrtTransformTimes(float start, float end) {
    VERIFY_OPTIONS("TransformTimes");
    renderOptions->transformStartTime = mfloat<time_d>(start);
    renderOptions->transformEndTime = mfloat<time_d>(end);
}


void pbrtPixelFilter(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("PixelFilter");
    renderOptions->FilterName = name;
    renderOptions->FilterParams = params;
}


void pbrtFilm(const string &type, const ParamSet &params) {
    VERIFY_OPTIONS("Film");
    renderOptions->FilmParams = params;
    renderOptions->FilmName = type;
}


void pbrtSampler(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Sampler");
    renderOptions->SamplerName = name;
    renderOptions->SamplerParams = params;
}


void pbrtAccelerator(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Accelerator");
    renderOptions->AcceleratorName = name;
    renderOptions->AcceleratorParams = params;
}


void pbrtSurfaceIntegrator(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("SurfaceIntegrator");
    renderOptions->SurfIntegratorName = name;
    renderOptions->SurfIntegratorParams = params;
}


void pbrtVolumeIntegrator(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("VolumeIntegrator");
    renderOptions->VolIntegratorName = name;
    renderOptions->VolIntegratorParams = params;
}


void pbrtRenderer(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Renderer");
    renderOptions->RendererName = name;
    renderOptions->RendererParams = params;
}


void pbrtCamera(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Camera");
    renderOptions->CameraName = name;
    renderOptions->CameraParams = params;
    renderOptions->CameraToWorld = Inverse(curTransform);
    namedCoordinateSystems["camera"] = renderOptions->CameraToWorld;
}


void pbrtWorldBegin() {
    VERIFY_OPTIONS("WorldBegin");
    currentApiState = STATE_WORLD_BLOCK;
    for (int i = 0; i < MAX_TRANSFORMS; ++i)
        curTransform.get<unknown_s, unknown_s>(i) = Transform<unknown_s, unknown_s>();
    activeTransformBits = ALL_TRANSFORMS_BITS;
    namedCoordinateSystems["world"] = curTransform;
}


void pbrtAttributeBegin() {
    VERIFY_WORLD("AttributeBegin");
    pushedGraphicsStates.push_back(graphicsState);
    pushedTransforms.push_back(curTransform);
    pushedActiveTransformBits.push_back(activeTransformBits);
}


void pbrtAttributeEnd() {
    VERIFY_WORLD("AttributeEnd");
    if (!pushedGraphicsStates.size()) {
        Error("Unmatched pbrtAttributeEnd() encountered. "
              "Ignoring it.");
        return;
    }
    graphicsState = pushedGraphicsStates.back();
    pushedGraphicsStates.pop_back();
    curTransform = pushedTransforms.back();
    pushedTransforms.pop_back();
    activeTransformBits = pushedActiveTransformBits.back();
    pushedActiveTransformBits.pop_back();
}


void pbrtTransformBegin() {
    VERIFY_WORLD("TransformBegin");
    pushedTransforms.push_back(curTransform);
    pushedActiveTransformBits.push_back(activeTransformBits);
}


void pbrtTransformEnd() {
    VERIFY_WORLD("TransformEnd");
    if (!pushedTransforms.size()) {
        Error("Unmatched pbrtTransformEnd() encountered. "
            "Ignoring it.");
        return;
    }
    curTransform = pushedTransforms.back();
    pushedTransforms.pop_back();
    activeTransformBits = pushedActiveTransformBits.back();
    pushedActiveTransformBits.pop_back();
}


void pbrtTexture(const string &name, const string &type,
                 const string &texname, const ParamSet &params) {
    VERIFY_WORLD("Texture");
    TextureParams tp(params, params, graphicsState.floatTextures,
                     graphicsState.spectrumTextures);
    if (type == "float")  {
        // Create _float_ texture and store in _floatTextures_
        if (graphicsState.floatTextures.find(name) !=
            graphicsState.floatTextures.end())
            Warning("Texture \"%s\" being redefined", name.c_str());
        WARN_IF_ANIMATED_TRANSFORM("Texture");
        Reference<Texture<mfloat<unit_d> > > ft = MakeFloatTexture(texname,
                                                         curTransform.get<texture_s, world_s>(0), tp);
        if (ft) graphicsState.floatTextures[name] = ft;
    }
    else if (type == "color")  {
        // Create _color_ texture and store in _spectrumTextures_
        if (graphicsState.spectrumTextures.find(name) != graphicsState.spectrumTextures.end())
            Warning("Texture \"%s\" being redefined", name.c_str());
        WARN_IF_ANIMATED_TRANSFORM("Texture");
        Reference<Texture<Spectrum<rho_d> > > st = MakeSpectrumTexture(texname,
            curTransform.get<texture_s, world_s>(0), tp);
        if (st) graphicsState.spectrumTextures[name] = st;
    }
    else
        Error("Texture type \"%s\" unknown.", type.c_str());
}


void pbrtMaterial(const string &name, const ParamSet &params) {
    VERIFY_WORLD("Material");
    graphicsState.material = name;
    graphicsState.materialParams = params;
    graphicsState.currentNamedMaterial = "";
}


void pbrtMakeNamedMaterial(const string &name,
        const ParamSet &params) {
    VERIFY_WORLD("MakeNamedMaterial");
    // error checking, warning if replace, what to use for transform?
    TextureParams mp(params, graphicsState.materialParams,
                     graphicsState.floatTextures,
                     graphicsState.spectrumTextures);
    string matName = mp.FindString("type");
    WARN_IF_ANIMATED_TRANSFORM("MakeNamedMaterial");
    if (matName == "") Error("No parameter string \"type\" found in MakeNamedMaterial");
    else {
        Reference<Material> mtl = MakeMaterial(matName, curTransform.get<material_s, world_s>(0), mp);
        if (mtl) graphicsState.namedMaterials[name] = mtl;
    }
}



void pbrtNamedMaterial(const string &name) {
    VERIFY_WORLD("NamedMaterial");
    graphicsState.currentNamedMaterial = name;
}


void pbrtLightSource(const string &name, const ParamSet &params) {
    VERIFY_WORLD("LightSource");
    WARN_IF_ANIMATED_TRANSFORM("LightSource");
    Light *lt = MakeLight(name, curTransform.get<light_s, world_s>(0), params);
    if (lt == NULL)
        Error("pbrtLightSource: light type \"%s\" unknown.", name.c_str());
    else
        renderOptions->lights.push_back(lt);
}


void pbrtAreaLightSource(const string &name,
                         const ParamSet &params) {
    VERIFY_WORLD("AfloatightSource");
    graphicsState.areaLight = name;
    graphicsState.afloatightParams = params;
}


void pbrtShape(const string &name, const ParamSet &params) {
    VERIFY_WORLD("Shape");
    Reference<Primitive> prim;
    AreaLight *area = NULL;
    if (!curTransform.IsAnimated()) {
        // Create primitive for static shape
        Transform<object_s, world_s> *obj2world;
        Transform<world_s, object_s> *world2obj;
        transformCache.Lookup(curTransform.get<object_s, world_s>(0), &obj2world, &world2obj);
        Reference<Shape> shape = MakeShape(name, obj2world, world2obj,
            graphicsState.reverseOrientation, params);
        if (!shape) return;
        Reference<Material> mtl = graphicsState.CreateMaterial(params);
        params.ReportUnused();

        // Possibly create area light for shape
        if (graphicsState.areaLight != "") {
            area = MakeAreaLight(graphicsState.areaLight, curTransform.get<light_s, world_s>(0),
                                 graphicsState.afloatightParams, shape);
        }
        prim = new GeometricPrimitive(shape, mtl, area);
    }
    else {
        // Create primitive for animated shape
        // Create initial _Shape_ for animated shape
        if (graphicsState.areaLight != "")
            Warning("Ignoring currently set area light when creating "
                    "animated shape");
        Transform<object_s, world_s> *o2w_identity;
        Transform<world_s, object_s> *w2o_identity;
        transformCache.Lookup(Transform<object_s, world_s>(), &o2w_identity, &w2o_identity);
        Reference<Shape> shape = MakeShape(name, o2w_identity, w2o_identity,
            graphicsState.reverseOrientation, params);
        if (!shape) return;
        Reference<Material> mtl = graphicsState.CreateMaterial(params);
        params.ReportUnused();

        // Get _animatedWorldToObject_ transform for shape
        Assert(MAX_TRANSFORMS == 2);
        Transform<world_s, world_s> *world2obj[2];
        transformCache.Lookup(curTransform.get<world_s, world_s>(0), reinterpret_cast<Transform<world_s, world_s>**>(NULL), &world2obj[0]);
        transformCache.Lookup(curTransform.get<world_s, world_s>(1), reinterpret_cast<Transform<world_s, world_s>**>(NULL), &world2obj[1]);
        AnimatedTransform<world_s, world_s>
             animatedWorldToObject(world2obj[0], renderOptions->transformStartTime,
                                   world2obj[1], renderOptions->transformEndTime);
        Reference<Primitive> baseprim = new GeometricPrimitive(shape, mtl, NULL);
        if (!baseprim->CanIntersect()) {
            // Refine animated shape and create BVH if more than one shape created
            std::vector<Reference<Primitive> > refinedPrimitives;
            baseprim->FullyRefine(refinedPrimitives);
            if (refinedPrimitives.size() == 0) return;
            if (refinedPrimitives.size() > 1)
                baseprim = new BVHAccel(refinedPrimitives);
            else
                baseprim = refinedPrimitives[0];
        }
        prim = new TransformedPrimitive(baseprim, animatedWorldToObject);
    }
    // Add primitive to scene or current instance
    if (renderOptions->currentInstance) {
        if (area)
            Warning("Area lights not supported with object instancing");
        renderOptions->currentInstance->push_back(prim);
    }
    else {
        renderOptions->primitives.push_back(prim);
        if (area != NULL) {
            renderOptions->lights.push_back(area);
        }
    }
}


Reference<Material> GraphicsState::CreateMaterial(const ParamSet &params) {
    TextureParams mp(params, materialParams,
                     floatTextures,
                     spectrumTextures);
    Reference<Material> mtl;
    if (currentNamedMaterial != "" &&
        namedMaterials.find(currentNamedMaterial) != namedMaterials.end())
        mtl = namedMaterials[graphicsState.currentNamedMaterial];
    if (!mtl)
        mtl = MakeMaterial(material, curTransform.get<material_s, world_s>(0), mp);
    if (!mtl)
        mtl = MakeMaterial("matte", curTransform.get<material_s, world_s>(0), mp);
    if (!mtl)
        Severe("Unable to create \"matte\" material?!");
    return mtl;
}


void pbrtReverseOrientation() {
    VERIFY_WORLD("ReverseOrientation");
    graphicsState.reverseOrientation =
        !graphicsState.reverseOrientation;
}


void pbrtVolume(const string &name, const ParamSet &params) {
    VERIFY_WORLD("Volume");
    WARN_IF_ANIMATED_TRANSFORM("Volume");
    VolumeRegion *vr = MakeVolumeRegion(name, curTransform.get<volume_s, world_s>(0), params);
    if (vr) renderOptions->volumeRegions.push_back(vr);
}


void pbrtObjectBegin(const string &name) {
    VERIFY_WORLD("ObjectBegin");
    pbrtAttributeBegin();
    if (renderOptions->currentInstance)
        Error("ObjectBegin called inside of instance definition");
    renderOptions->instances[name] = std::vector<Reference<Primitive> >();
    renderOptions->currentInstance = &renderOptions->instances[name];
}


void pbrtObjectEnd() {
    VERIFY_WORLD("ObjectEnd");
    if (!renderOptions->currentInstance)
        Error("ObjectEnd called outside of instance definition");
    renderOptions->currentInstance = NULL;
    pbrtAttributeEnd();
}


void pbrtObjectInstance(const string &name) {
    VERIFY_WORLD("ObjectInstance");
    // Object instance error checking
    if (renderOptions->currentInstance) {
        Error("ObjectInstance can't be called inside instance definition");
        return;
    }
    if (renderOptions->instances.find(name) == renderOptions->instances.end()) {
        Error("Unable to find instance named \"%s\"", name.c_str());
        return;
    }
    std::vector<Reference<Primitive> > &in = renderOptions->instances[name];
    if (in.size() == 0) return;
    if (in.size() > 1 || !in[0]->CanIntersect()) {
        // Refine instance _Primitive_s and create aggregate
        Reference<Primitive> accel =
             MakeAccelerator(renderOptions->AcceleratorName,
                             in, renderOptions->AcceleratorParams);
        if (!accel) accel = MakeAccelerator("bvh", in, ParamSet());
        if (!accel) Severe("Unable to create \"bvh\" accelerator");
        in.erase(in.begin(), in.end());
        in.push_back(accel);
    }
    Assert(MAX_TRANSFORMS == 2);
    Transform<world_s, world_s> *world2instance[2];
    transformCache.Lookup(curTransform.get<world_s, world_s>(0), reinterpret_cast<Transform<world_s, world_s>**>(NULL), &world2instance[0]);
    transformCache.Lookup(curTransform.get<world_s, world_s>(0), reinterpret_cast<Transform<world_s, world_s>**>(NULL), &world2instance[1]);
    AnimatedTransform<world_s, world_s> animatedWorldToInstance(world2instance[0],
        renderOptions->transformStartTime,
        world2instance[1], renderOptions->transformEndTime);
    Reference<Primitive> prim =
        new TransformedPrimitive(in[0], animatedWorldToInstance);
    renderOptions->primitives.push_back(prim);
}


void pbrtWorldEnd() {
    VERIFY_WORLD("WorldEnd");
    // Ensure there are no pushed graphics states
    while (pushedGraphicsStates.size()) {
        Warning("Missing end to pbrtAttributeBegin()");
        pushedGraphicsStates.pop_back();
        pushedTransforms.pop_back();
    }
    while (pushedTransforms.size()) {
        Warning("Missing end to pbrtTransformBegin()");
        pushedTransforms.pop_back();
    }

    // Create scene and render
    Renderer *renderer = renderOptions->MakeRenderer();
    Scene *scene = renderOptions->MakeScene();
    if (scene && renderer) renderer->Render(scene);
    TasksCleanup();
    delete renderer;
    delete scene;

    // Clean up after rendering
    graphicsState = GraphicsState();
    transformCache.Clear();
    currentApiState = STATE_OPTIONS_BLOCK;
    ProbesPrint(stdout);
    for (int i = 0; i < MAX_TRANSFORMS; ++i)
        curTransform.get<unknown_s, unknown_s>(i) = Transform<unknown_s, unknown_s>();
    activeTransformBits = ALL_TRANSFORMS_BITS;
    namedCoordinateSystems.erase(namedCoordinateSystems.begin(),
                                 namedCoordinateSystems.end());
}


Scene *RenderOptions::MakeScene() {
    // Initialize _volumeRegion_ from volume region(s)
    VolumeRegion *volumeRegion = NULL;
    if (volumeRegions.size() == 0)
        volumeRegion = NULL;
    else if (volumeRegions.size() == 1)
        volumeRegion = volumeRegions[0];
    else
        volumeRegion = new AggregateVolume(volumeRegions);
    Primitive *accelerator = MakeAccelerator(AcceleratorName,
        primitives, AcceleratorParams);
    if (!accelerator)
        accelerator = MakeAccelerator("bvh", primitives, ParamSet());
    if (!accelerator)
        Severe("Unable to create \"bvh\" accelerator.");
    Scene *scene = new Scene(accelerator, lights, volumeRegion);
    // Erase primitives, lights, and volume regions from _RenderOptions_
    primitives.erase(primitives.begin(), primitives.end());
    lights.erase(lights.begin(), lights.end());
    volumeRegions.erase(volumeRegions.begin(), volumeRegions.end());
    return scene;
}


Renderer *RenderOptions::MakeRenderer() const {
    Renderer *renderer = NULL;
    Camera *camera = MakeCamera();
    if (RendererName != "sample")
        Warning("Renderer type \"%s\" unknown.  Using standard.",
        RendererName.c_str());
    RendererParams.ReportUnused();
    Sampler *sampler = MakeSampler(SamplerName, SamplerParams, camera->film, camera);
    if (!sampler) Severe("Unable to create sampler.");
    // Create surface and volume integrators
    SurfaceIntegrator *surfaceIntegrator = MakeSurfaceIntegrator(SurfIntegratorName,
        SurfIntegratorParams);
    if (!surfaceIntegrator) Severe("Unable to create surface integrator.");
    VolumeIntegrator *volumeIntegrator = MakeVolumeIntegrator(VolIntegratorName,
        VolIntegratorParams);
    if (!volumeIntegrator) Severe("Unable to create volume integrator.");
    renderer = new SamplerRenderer(sampler, camera, surfaceIntegrator,
        volumeIntegrator);
    return renderer;
}


Camera *RenderOptions::MakeCamera() const {
    Filter *filter = MakeFilter(FilterName, FilterParams);
    Film *film = MakeFilm(FilmName, FilmParams, filter);
    if (!film) Severe("Unable to create film.");
    Camera *camera = ::MakeCamera(CameraName, CameraParams,
        CameraToWorld, renderOptions->transformStartTime,
        renderOptions->transformEndTime, film);
    if (!camera) Severe("Unable to create camera.");
    return camera;
}


