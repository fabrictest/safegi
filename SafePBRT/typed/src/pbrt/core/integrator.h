
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

#ifndef PBRT_CORE_INTEGRATOR_H
#define PBRT_CORE_INTEGRATOR_H

// core/integrator.h*
#include "pbrt.h"
#include "primitive.h"
#include "spectrum.h"
#include "reflection.h"
#include "sampler.h"
#include "material.h"
#include "renderer.h"

// Integrator Declarations
class Integrator {
public:
    // Integrator Interface
    virtual ~Integrator();
    virtual void Preprocess(const Scene *scene, const Camera *camera,
                            const Renderer *renderer) {
    }
    virtual void RequestSamples(Sampler *sampler, Sample *sample,
                                const Scene *scene) {
    }
};


class SurfaceIntegrator : public Integrator {
public:
    // SurfaceIntegrator Interface
    virtual Spectrum<radiance_d> Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential<world_s> &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const = 0;
};


Spectrum<radiance_d> UniformSampleAllLights(const Scene *scene, const Renderer *renderer,
    MemoryArena &arena, const Point<world_s> &p, const Normal<world_s> &n, const Direction<world_s> &wo,
    const mfloat<length_d> &rayEpsilon, const mfloat<time_d> &time, BSDF *bsdf, const Sample *sample, RNG &rng,
    const LightSampleOffsets *lightSampleOffsets,
    const BSDFSampleOffsets *bsdfSampleOffsets);
Spectrum<radiance_d> UniformSampleOneLight(const Scene *scene, const Renderer *renderer,
    MemoryArena &arena, const Point<world_s> &p, const Normal<world_s> &n, const Direction<world_s> &wo,
    const mfloat<length_d> &rayEpsilon, const mfloat<time_d> &time, BSDF *bsdf,
    const Sample *sample, RNG &rng, int lightNumOffset = -1,
    const LightSampleOffsets *lightSampleOffset = NULL,
    const BSDFSampleOffsets *bsdfSampleOffset = NULL);
Spectrum<radiance_d> EstimateDirect(const Scene *scene, const Renderer *renderer,
    MemoryArena &arena, const Light *light, const Point<world_s> &p,
    const Normal<world_s> &n, const Direction<world_s> &wo, const mfloat<length_d> &rayEpsilon, const mfloat<time_d> &time, BSDF *bsdf,
    RNG &rng, const LightSample &lightSample, const BSDFSample &bsdfSample);
Spectrum<radiance_d> SpecularReflect(const RayDifferential<world_s> &ray, BSDF *bsdf, RNG &rng,
    const Intersection &isect, const Renderer *renderer, const Scene *scene,
    const Sample *sample, MemoryArena &arena);
Spectrum<radiance_d> SpecularTransmit(const RayDifferential<world_s> &ray, BSDF *bsdf, RNG &rng,
    const Intersection &isect, const Renderer *renderer, const Scene *scene,
    const Sample *sample, MemoryArena &arena);
Distribution1D<power_d> *ComputeLightSamplingCDF(const Scene *scene);

#endif // PBRT_CORE_INTEGRATOR_H
