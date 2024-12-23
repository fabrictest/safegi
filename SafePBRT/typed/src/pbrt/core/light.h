
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

#ifndef PBRT_CORE_LIGHT_H
#define PBRT_CORE_LIGHT_H

// core/light.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
#include "spectrum.h"
#include "rng.h"
#include "memory.h"
#include "scene.h"
#include "renderer.h"

// Light Declarations
class Light {
public:
    // Light Interface
    virtual ~Light();
    Light(const Transform<light_s, world_s> &l2w, int ns = 1)
        : nSamples(max(1, ns)), LightToWorld(l2w),
          WorldToLight(Inverse(l2w)) {
        // Warn if light has transformation with scale
        if (WorldToLight.HasScale())
            Warning("Scaling detected in world to light transformation!\n"
                    "The system has numerous assumptions, implicit and explicit,\n"
                    "that this transform will have no scale factors in it.\n"
                    "Proceed at your own risk; your image may have errors or\n"
                    "the system may crash as a result of this.");
    }
    virtual Spectrum<radiance_d> Sample_L(const Point<world_s> &p, const mfloat<length_d> &pEpsilon,
        const LightSample &ls, const mfloat<time_d> &time, Direction<world_s> *wi, mfloat<invsolidangle_d> *pdf,
        VisibilityTester<world_s> *vis) const = 0;
    virtual Spectrum<power_d> Power(const Scene *) const = 0;
    virtual bool IsDeltaLight() const = 0;
    virtual Spectrum<radiance_d> Le(const RayDifferential<world_s> &r) const;
    virtual mfloat<invsolidangle_d> Pdf(const Point<world_s> &p, const Direction<world_s> &wi) const = 0;
    virtual Spectrum<radiance_d> Sample_L(const Scene *scene, const LightSample &ls,
                              float u1, float u2, const mfloat<time_d> &time, Ray<world_s> *ray,
                              Normal<world_s> *Ns, mfloat<invsolidanglearea_d> *pdf) const = 0;

    // Light Public Data
    const int nSamples;
protected:
    // Light Protected Data
    const Transform<light_s, world_s> LightToWorld;
    const Transform<world_s, light_s> WorldToLight;
};

template<typename S>
struct VisibilityTester {
    // VisibilityTester Public Methods
    void SetSegment(const Point<S> &p1, const mfloat<length_d> &eps1,
                    const Point<S> &p2, const float &eps2, const mfloat<time_d> &time) {
        mfloat<length_d> dist = Distance(p1, p2);
        Vector<S> rv = p2 - p1;
        r = Ray<S>(p1, Normalize(rv), eps1, rv.Length() * (1.f - eps2), time);
        Assert(!r.HasNaNs());
    }
    void SetRay(const Point<S> &p, const mfloat<length_d> & eps, const Direction<S> &w, const mfloat<time_d> &time) {
        r = Ray<S>(p, w, eps, mfloat<length_d>(INFINITY), time);
        Assert(!r.HasNaNs());
    }

    bool Unoccluded(const Scene *scene) const {
        return !scene->IntersectP(r);
    }

    Spectrum<rho_d> Transmittance(const Scene *scene,
        const Renderer *renderer, const Sample *sample,
        RNG &rng, MemoryArena &arena) const {
            return renderer->Transmittance(scene, RayDifferential<S>(r), sample,
                rng, arena);
    }
    Ray<S> r;
};


class AreaLight : public Light {
public:
    // AreaLight Interface
    AreaLight(const Transform<light_s, world_s> &l2w, int ns) : Light(l2w, ns) { }
    virtual Spectrum<radiance_d> L(const Point<world_s> &p, const Normal<world_s> &n,
                       const Direction<world_s> &w) const = 0;
};


struct LightSample {
   // LightSample Public Methods
   LightSample() { }
   LightSample(const Sample *sample, const LightSampleOffsets &offsets, uint32_t num);
   LightSample(RNG &rng) {
       uPos[0] = rng.RandomFloat();
       uPos[1] = rng.RandomFloat();
       uComponent = rng.RandomFloat();
   }
   LightSample(float up0, float up1, float ucomp) {
       uPos[0] = up0; uPos[1] = up1;
       uComponent = ucomp;
   }
   float uPos[2], uComponent;
};


struct LightSampleOffsets {
    LightSampleOffsets() { }
    LightSampleOffsets(int count, Sample *sample);
    int nSamples, componentOffset, posOffset;
};



// ShapeSet Declarations
class ShapeSet {
public:
    // ShapeSet Public Methods
    ShapeSet(const Reference<Shape> &s);
    mfloat<area_d> Area() const { return sumArea; }
    ~ShapeSet();
    Point<world_s> Sample(const Point<world_s> &p, const LightSample &ls, Normal<world_s> *Ns) const;
    Point<world_s> Sample(const LightSample &ls, Normal<world_s> *Ns) const;
    mfloat<invsolidangle_d> Pdf(const Point<world_s> &p, const Direction<world_s> &wi) const;
    mfloat<invarea_d> Pdf(const Point<world_s> &p) const;
private:
    // ShapeSet Private Data
    std::vector<Reference<Shape> > shapes;
    mfloat<area_d> sumArea;
    std::vector<mfloat<area_d> > areas;
    Distribution1D<area_d> *areaDistribution;
};



#endif // PBRT_CORE_LIGHT_H
