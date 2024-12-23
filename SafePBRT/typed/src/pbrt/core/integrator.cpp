
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


// core/integrator.cpp*
#include "integrator.h"
#include "scene.h"
#include "light.h"
#include "intersection.h"
#include "montecarlo.h"

// Integrator Method Definitions
Integrator::~Integrator() {
}



// Integrator Utility Functions
Spectrum<radiance_d> UniformSampleAllLights(const Scene *scene,
        const Renderer *renderer, MemoryArena &arena, const Point<world_s> &p,
        const Normal<world_s> &n, const Direction<world_s> &wo, const mfloat<length_d> &rayEpsilon,
        const mfloat<time_d> &time, BSDF *bsdf, const Sample *sample, RNG &rng,
        const LightSampleOffsets *lightSampleOffsets,
        const BSDFSampleOffsets *bsdfSampleOffsets) {
    Spectrum<radiance_d> L;
    for (uint32_t i = 0; i < scene->lights.size(); ++i) {
        Light *light = scene->lights[i];
        int nSamples = lightSampleOffsets ?
                       lightSampleOffsets[i].nSamples : 1;
        // Estimate direct lighting from _light_ samples
        Spectrum<radiance_d> Ld;
        for (int j = 0; j < nSamples; ++j) {
            // Find light and BSDF sample values for direct lighting estimate
            LightSample lightSample;
            BSDFSample bsdfSample;
            if (lightSampleOffsets != NULL && bsdfSampleOffsets != NULL) {
                lightSample = LightSample(sample, lightSampleOffsets[i], j);
                bsdfSample = BSDFSample(sample, bsdfSampleOffsets[i], j);
            }
            else {
                lightSample = LightSample(rng);
                bsdfSample = BSDFSample(rng);
            }
            Ld += EstimateDirect(scene, renderer, arena, light, p, n, wo,
                rayEpsilon, time, bsdf, rng, lightSample, bsdfSample);
        }
        L += Ld / nSamples;
    }
    return L;
}


Spectrum<radiance_d> UniformSampleOneLight(const Scene *scene,
        const Renderer *renderer, MemoryArena &arena, const Point<world_s> &p,
        const Normal<world_s> &n, const Direction<world_s> &wo, const mfloat<length_d> &rayEpsilon, const mfloat<time_d> &time,
        BSDF *bsdf, const Sample *sample, RNG &rng, int lightNumOffset,
        const LightSampleOffsets *lightSampleOffset,
        const BSDFSampleOffsets *bsdfSampleOffset) {
    // Randomly choose a single light to sample, _light_
    int nLights = int(scene->lights.size());
    if (nLights == 0.) return Spectrum<radiance_d>();
    int lightNum;
    if (lightNumOffset != -1)
        lightNum = Floor2Int(sample->oneD[lightNumOffset][0] * nLights);
    else
        lightNum = Floor2Int(rng.RandomFloat() * nLights);
    lightNum = min(lightNum, nLights-1);
    Light *light = scene->lights[lightNum];

    // Initialize light and bsdf samples for single light sample
    LightSample lightSample;
    BSDFSample bsdfSample;
    if (lightSampleOffset != NULL && bsdfSampleOffset != NULL) {
        lightSample = LightSample(sample, *lightSampleOffset, 0);
        bsdfSample = BSDFSample(sample, *bsdfSampleOffset, 0);
    }
    else {
        lightSample = LightSample(rng);
        bsdfSample = BSDFSample(rng);
    }
    return (float)nLights *
        EstimateDirect(scene, renderer, arena, light, p, n, wo,
                       rayEpsilon, time, bsdf, rng, lightSample,
                       bsdfSample);
}


Spectrum<radiance_d> EstimateDirect(const Scene *scene, const Renderer *renderer,
        MemoryArena &arena, const Light *light, const Point<world_s> &p,
        const Normal<world_s> &n, const Direction<world_s> &wo, const mfloat<length_d> &rayEpsilon, const mfloat<time_d> &time,
        BSDF *bsdf, RNG &rng, const LightSample &lightSample,
        const BSDFSample &bsdfSample) {
    Spectrum<radiance_d> Ld;
    // Sample light source with multiple importance sampling
    Direction<world_s> wi;
    mfloat<invsolidangle_d> lightPdf, bsdfPdf;
    VisibilityTester<world_s> visibility;
    Spectrum<radiance_d> Li = light->Sample_L(p, rayEpsilon, lightSample, time,
                                  &wi, &lightPdf, &visibility);
    if (lightPdf > mfloat<invsolidangle_d>(0.0f) && !Li.IsBlack()) {
        Spectrum<brdf_d> f = bsdf->f(wo, wi);
        if (!f.IsBlack() && visibility.Unoccluded(scene)) {
            // Add light's contribution to reflected radiance
            Li = Li * visibility.Transmittance(scene, renderer, NULL, rng, arena);
            if (light->IsDeltaLight())
                Ld += f * Li * AbsDot(wi, n) / lightPdf;
            else {
                bsdfPdf = bsdf->Pdf(wo, wi);
                float weight = PowerHeuristic(1, lightPdf, 1, bsdfPdf);
                Ld += f * Li * AbsDot(wi, n) * weight / lightPdf;
            }
        }
    }

    // Sample BSDF with multiple importance sampling
    if (!light->IsDeltaLight()) {
        BxDFType flags = BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
        Spectrum<brdf_d> f = bsdf->Sample_f(wo, &wi, bsdfSample, &bsdfPdf, flags);
        if (!f.IsBlack() && bsdfPdf > mfloat<invsolidangle_d>(0.0f)) {
            lightPdf = light->Pdf(p, wi);
            if (lightPdf > mfloat<invsolidangle_d>(0.0f)) {
                // Add light contribution from BSDF sampling
                float weight = PowerHeuristic(1, bsdfPdf, 1, lightPdf);
                Intersection lightIsect;
                Spectrum<radiance_d> Li;
                RayDifferential<world_s> ray(p, wi, rayEpsilon, INFINITY*meters, time);
                if (scene->Intersect(ray, &lightIsect)) {
                    if (lightIsect.primitive->GetAreaLight() == light)
                        Li = lightIsect.Le(-wi);
                }
                else
                    Li = light->Le(ray);
                if (!Li.IsBlack()) {
                    Li = Li * renderer->Transmittance(scene, ray, NULL, rng, arena);
                    Ld += f * Li * AbsDot(wi, n) * weight / bsdfPdf;
                }
            }
        }
    }
    return Ld;
}


Spectrum<radiance_d> SpecularReflect(const RayDifferential<world_s> &ray, BSDF *bsdf,
        RNG &rng, const Intersection &isect, const Renderer *renderer,
        const Scene *scene, const Sample *sample, MemoryArena &arena) {
    Direction<world_s> wo = -ray.d, wi;
    mfloat<invsolidangle_d> pdf;
    const Point<world_s> &p = bsdf->dgShading.p;
    const Normal<world_s> &n = bsdf->dgShading.nn;
    Spectrum<brdf_d> f = bsdf->Sample_f(wo, &wi, BSDFSample(rng), &pdf,
                                BxDFType(BSDF_REFLECTION | BSDF_SPECULAR));
    Spectrum<radiance_d> L;
    if (pdf > mfloat<invsolidangle_d>(0.0f) && !f.IsBlack() && AbsDot(wi, n) != mfloat<proj_d>(0.0f)) {
        // Compute ray differential _rd_ for specular reflection
        RayDifferential<world_s> rd(p, wi, ray, isect.rayEpsilon);
        if (ray.hasDifferentials) {
            rd.hasDifferentials = true;
            rd.rxOrigin = p + isect.dg.dpdx;
            rd.ryOrigin = p + isect.dg.dpdy;

            // Compute differential reflected directions
            auto dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx + 
                        bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
            auto dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy + 
                        bsdf->dgShading.dndv * bsdf->dgShading.dvdy;
            auto dwodx = (-ray.rxDirection) - wi, dwody = (-ray.ryDirection) - wo;
            auto dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
            auto dDNdy = Dot(dwody, n) + Dot(wo, dndy);
            rd.rxDirection = 
                Normalize<world_s>(wi - dwodx + 2 * (Dot(wo, n) * dndx + dDNdx * n));
            rd.ryDirection = 
                Normalize<world_s>(wi - dwody + 2 * (Dot(wo, n) * dndy + dDNdy * n));
        }
        PBRT_STARTED_SPECULAR_REFLECTION_RAY(const_cast<RayDifferential<world_s> *>(&rd));
        Spectrum<radiance_d> Li = renderer->Li(scene, rd, sample, rng, arena);
        L = f * Li * AbsDot(wi, n) / pdf;
        PBRT_FINISHED_SPECULAR_REFLECTION_RAY(const_cast<RayDifferential<world_s> *>(&rd));
    }
    return L;
}


Spectrum<radiance_d> SpecularTransmit(const RayDifferential<world_s> &ray, BSDF *bsdf,
        RNG &rng, const Intersection &isect, const Renderer *renderer,
        const Scene *scene, const Sample *sample, MemoryArena &arena) {
    Direction<world_s> wo = -ray.d, wi;
    mfloat<invsolidangle_d> pdf;
    const Point<world_s> &p = bsdf->dgShading.p;
    const Normal<world_s> &n = bsdf->dgShading.nn;
    Spectrum<brdf_d> f = bsdf->Sample_f(wo, &wi, BSDFSample(rng), &pdf,
                               BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));
    Spectrum<radiance_d> L;
    if (pdf > mfloat<invsolidangle_d>(0.0f) && !f.IsBlack() && AbsDot(wi, n) != mfloat<proj_d>(0.0f)) {
        // Compute ray differential _rd_ for specular transmission
        RayDifferential<world_s> rd(p, wi, ray, isect.rayEpsilon);
        if (ray.hasDifferentials) {
            rd.hasDifferentials = true;
            rd.rxOrigin = p + isect.dg.dpdx;
            rd.ryOrigin = p + isect.dg.dpdy;
        
            float eta = bsdf->eta;
            Direction<world_s> w = -wo;
            if (Dot(wo, n) < 0) eta = 1.f / eta;
        
            auto dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx + bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
            auto dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy + bsdf->dgShading.dndv * bsdf->dgShading.dvdy;
        
            auto dwodx = -ray.rxDirection - wo, dwody = -ray.ryDirection - wo;
            auto dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
            auto dDNdy = Dot(dwody, n) + Dot(wo, dndy);
        
            float mu = eta * Dot(w, n) - Dot(wi, n);
            float dmudx = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdx;
            float dmudy = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdy;
        
            rd.rxDirection = Normalize<world_s>(wi + eta * dwodx - (mu * dndx + dmudx * n));
            rd.ryDirection = Normalize<world_s>(wi + eta * dwody - (mu * dndy + dmudy * n));
        }
        PBRT_STARTED_SPECULAR_REFRACTION_RAY(const_cast<RayDifferential<world_s> *>(&rd));
        Spectrum<radiance_d> Li = renderer->Li(scene, rd, sample, rng, arena);
        L = f * Li * AbsDot(wi, n) / pdf;
        PBRT_FINISHED_SPECULAR_REFRACTION_RAY(const_cast<RayDifferential<world_s> *>(&rd));
    }
    return L;
}


Distribution1D<power_d> *ComputeLightSamplingCDF(const Scene *scene) {
    uint32_t nLights = int(scene->lights.size());
    Assert(nLights > 0);
    std::vector<mfloat<power_d> >lightPower(nLights, mfloat<power_d>(0.0f));
    for (uint32_t i = 0; i < nLights; ++i)
        lightPower[i] = scene->lights[i]->Power(scene).y();
    return new Distribution1D<power_d>(&lightPower[0], nLights);
}


