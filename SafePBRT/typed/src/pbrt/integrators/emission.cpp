
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


// integrators/emission.cpp*
#include "integrators/emission.h"
#include "paramset.h"

// EmissionIntegrator Method Definitions
void EmissionIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
                                        const Scene *scene) {
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}


Spectrum<rho_d> EmissionIntegrator::Transmittance(const Scene *scene,
        const Renderer *, const RayDifferential<world_s> &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    if (!scene->volumeRegion) return Spectrum<rho_d>(1.0f);
    mfloat<length_d> step;
    float offset;
    if (sample) {
        step = stepSize;
        offset = sample->oneD[tauSampleOffset][0];
    }
    else {
        step = 4.0f * stepSize;
        offset = rng.RandomFloat();
    }
    Spectrum<rho_d> tau = scene->volumeRegion->tau(ray, step, offset);
    return Exp(-tau);
}


Spectrum<radiance_d> EmissionIntegrator::Li(const Scene *scene,
        const Renderer *renderer, const RayDifferential<world_s> &ray,
        const Sample *sample, RNG &rng, Spectrum<rho_d> *T,
        MemoryArena &arena) const {
    VolumeRegion *vr = scene->volumeRegion;
    Assert(sample != NULL);
    mfloat<length_d> t0, t1;
    if (!vr || !vr->IntersectP(ray, &t0, &t1)) {
        *T = Spectrum<rho_d>(1.f);
        return Spectrum<radiance_d>();
    }
    // Do emission-only volume integration in _vr_
    Spectrum<radiancediff_d> Lv;

    // Prepare for volume integration stepping
    int nSamples = Ceil2Int((t1-t0) / stepSize);
    mfloat<length_d> step = (t1 - t0) / nSamples;
    Spectrum<rho_d> Tr(1.0f);
    Point<world_s> p = ray(t0), pPrev;
    Direction<world_s> w = -ray.d;
    t0 += sample->oneD[scatterSampleOffset][0] * step;
    for (int i = 0; i < nSamples; ++i, t0 += step) {
        // Advance to sample at _t0_ and update _T_
        pPrev = p;
        p = ray(t0);
        mfloat<length_d> len = Distance(p, pPrev);
        Spectrum<rho_d> stepTau;
        if(len != mfloat<length_d>(0.0f))
        {
            Ray<world_s> tauRay(pPrev, Normalize(p - pPrev), 0.f*len, 1.f*len, ray.time, ray.depth);
            stepTau = vr->tau(tauRay,
                .5f * stepSize, rng.RandomFloat());
        }
        Tr *= Exp(-stepTau);

        // Possibly terminate ray marching if transmittance is small
        if (Tr.y() < 1e-3) {
            const float continueProb = .5f;
            if (rng.RandomFloat() > continueProb) break;
            Tr /= continueProb;
        }

        // Compute emission-only source term at _p_
        Lv += Tr * vr->Lve(p, w, ray.time);
    }
    *T = Tr;
    return Lv * step;
}


EmissionIntegrator *CreateEmissionVolumeIntegrator(const ParamSet &params) {
    mfloat<length_d> stepSize  = params.FindOneFloat<length_d>("stepsize", mfloat<length_d>(1.0f));
    return new EmissionIntegrator(stepSize);
}


