
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

#ifndef PBRT_INTEGRATORS_DIRECTLIGHTING_H
#define PBRT_INTEGRATORS_DIRECTLIGHTING_H

// integrators/directlighting.h*
#include "pbrt.h"
#include "integrator.h"
#include "scene.h"

// DirectLightingIntegrator Declarations
enum LightStrategy { SAMPLE_ALL_UNIFORM, SAMPLE_ONE_UNIFORM };
class DirectLightingIntegrator : public SurfaceIntegrator {
public:
    // DirectLightingIntegrator Public Methods
    DirectLightingIntegrator(LightStrategy ls = SAMPLE_ALL_UNIFORM, int md = 5);
    ~DirectLightingIntegrator();
    Spectrum<radiance_d> Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential<world_s> &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
private:
    // DirectLightingIntegrator Private Data
    LightStrategy strategy;
    int maxDepth;

    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
    int lightNumOffset;
};


DirectLightingIntegrator *CreateDirectLightingIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_DIRECTLIGHTING_H
