
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

#ifndef PBRT_INTEGRATORS_PATH_H
#define PBRT_INTEGRATORS_PATH_H

// integrators/path.h*
#include "pbrt.h"
#include "integrator.h"
#include "light.h"

// PathIntegrator Declarations
class PathIntegrator : public SurfaceIntegrator {
public:
    // PathIntegrator Public Methods
    Spectrum<radiance_d> Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential<world_s> &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    PathIntegrator(int md) { maxDepth = md; }
private:
    // PathIntegrator Private Data
    int maxDepth;
#define SAMPLE_DEPTH 3
    LightSampleOffsets lightSampleOffsets[SAMPLE_DEPTH];
    int lightNumOffset[SAMPLE_DEPTH];
    BSDFSampleOffsets bsdfSampleOffsets[SAMPLE_DEPTH];
    BSDFSampleOffsets pathSampleOffsets[SAMPLE_DEPTH];
};


PathIntegrator *CreatePathSurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_PATH_H
