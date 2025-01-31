
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

#ifndef PBRT_INTEGRATORS_EMISSION_H
#define PBRT_INTEGRATORS_EMISSION_H

// integrators/emission.h*
#include "volume.h"
#include "integrator.h"
#include "scene.h"

// EmissionIntegrator Declarations
class EmissionIntegrator : public VolumeIntegrator {
public:
    // EmissionIntegrator Public Methods
    EmissionIntegrator(const mfloat<length_d> &ss) : stepSize(ss) { }
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    Spectrum<radiance_d> Li(const Scene *scene, const Renderer *renderer,
            const RayDifferential<world_s> &ray, const Sample *sample, RNG &rng,
            Spectrum<rho_d> *transmittance, MemoryArena &arena) const;
    Spectrum<rho_d> Transmittance(const Scene *scene, const Renderer *,
            const RayDifferential<world_s> &ray, const Sample *sample, RNG &rng,
            MemoryArena &arena) const;
private:
    // EmissionIntegrator Private Data
    mfloat<length_d> stepSize;
    int tauSampleOffset, scatterSampleOffset;
};


EmissionIntegrator *CreateEmissionVolumeIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_EMISSION_H
