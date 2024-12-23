
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


// lights/infinite.cpp*
#include "lights/infinite.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

// InfiniteAreaLight Utility Classes
struct InfiniteAreaCube {
    // InfiniteAreaCube Public Methods
    InfiniteAreaCube(const InfiniteAreaLight *l, const Scene *s,
                     float t, bool cv, float pe)
        : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum<radiance_d> operator()(int, int, const Point<world_s> &p, const Direction<world_s> &w) {
        Ray<world_s> ray(p, w, pEpsilon, mfloat<length_d>(INFINITY), time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential<world_s>(ray));
        return Spectrum<radiance_d>();
    }
    const InfiniteAreaLight *light;
    const Scene *scene;
    mfloat<time_d> time;
    mfloat<length_d> pEpsilon;
    bool computeVis;
};



// InfiniteAreaLight Method Definitions
InfiniteAreaLight::~InfiniteAreaLight() {
    delete distribution;
    delete radianceMap;
}


InfiniteAreaLight::InfiniteAreaLight(const Transform<light_s, world_s> &light2world,
        const Spectrum<unit_d> &L, int ns, const string &texmap)
    : Light(light2world, ns) {
    int width = 0, height = 0;
    RGBSpectrum<radiance_d> *texels = NULL;
    // Read texel data from _texmap_ into _texels_
    if (texmap != "") {
        texels = ReadImage<radiance_d>(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum<radiance_d>[1];
        texels[0] = mfloat<radiance_d>(1.0f) * L.ToRGBSpectrum();
    }
    radianceMap = new MIPMap<RGBSpectrum<radiance_d> >(width, height, texels);
    delete[] texels;
    // Initialize sampling PDFs for infinite area light

    // Compute scalar-valued image _img_ from environment map
    float filter = 1.f / max(width, height);
    mfloat<radiance_d> *img = new mfloat<radiance_d>[width*height];
    for (int v = 0; v < height; ++v) {
        float vp = (float)v / (float)height;
        mfloat<unit_d> sinTheta = sin(halfCircleAngle * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            float up = (float)u / (float)width;
            img[u+v*width] = sinTheta * radianceMap->Lookup(up, vp, filter).y();
        }
    }

    // Compute sampling distributions for rows and columns of image
    distribution = new Distribution2D<radiance_d>(img, width, height);
    delete[] img;
}


Spectrum<power_d> InfiniteAreaLight::Power(const Scene *scene) const {
    Point<world_s> worldCenter;
    mfloat<length_d> worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return hemisphericalProjectedAngle * worldRadius * worldRadius *
        Spectrum<radiance_d>(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);

}


Spectrum<radiance_d> InfiniteAreaLight::Le(const RayDifferential<world_s> &r) const {
    Direction<light_s> wh = WorldToLight(r.d);
    float s = SphericalPhi(wh) / circleAngle;
    float t = SphericalTheta(wh) / halfCircleAngle;
    return Spectrum<radiance_d>(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


Spectrum<radiance_d> InfiniteAreaLight::Sample_L(const Point<world_s> &p, const mfloat<length_d> &pEpsilon,
        const LightSample &ls, const mfloat<time_d> &time, Direction<world_s> *wi, mfloat<invsolidangle_d> *pdf,
        VisibilityTester<world_s> *visibility) const {
    // Find $(u,v)$ sample coordinates in infinite light texture
    float uv[2];
    float mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);

    // Convert infinite light sample point to direction
    mfloat<angle_d> theta = uv[1] * halfCircleAngle;
    mfloat<angle_d> phi = uv[0] * circleAngle;
    mfloat<unit_d> costheta = cos(theta), sintheta = sin(theta);
    mfloat<unit_d> sinphi = sin(phi), cosphi = cos(phi);
    *wi = LightToWorld(Direction<light_s>(sintheta * cosphi, sintheta * sinphi, costheta));

    // Compute PDF for sampled infinite light direction
    *pdf = mfloat<invsolidangle_d>(mapPdf / (2 * M_PI * M_PI * sintheta));
    if (sintheta == 0.f) *pdf = mfloat<invsolidangle_d>(0.0f);

    // Return radiance value for infinite light direction
    visibility->SetRay(p, pEpsilon, *wi, time);
    return Spectrum<radiance_d>(radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
}


mfloat<invsolidangle_d> InfiniteAreaLight::Pdf(const Point<world_s> &, const Direction<world_s> &w) const {
    Direction<light_s> wi = WorldToLight(w);
    mfloat<angle_d> theta = SphericalTheta(wi);
    mfloat<angle_d> phi = SphericalPhi(wi);
    mfloat<unit_d> sintheta = sin(theta);
    if (sintheta == 0.f) 
        return mfloat<invsolidangle_d>(0.f);
    return distribution->Pdf(phi / circleAngle, theta / halfCircleAngle) / steradians /
           (2.0f * M_PI * M_PI * sintheta);
}


Spectrum<radiance_d> InfiniteAreaLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, const mfloat<time_d> &time,
        Ray<world_s> *ray, Normal<world_s> *Ns, mfloat<invsolidanglearea_d> *pdf) const {
    // Compute direction for infinite light sample ray

    // Find $(u,v)$ sample coordinates in infinite light texture
    float uv[2];
    float mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    mfloat<angle_d> theta = uv[1] * halfCircleAngle, phi = uv[0] * circleAngle;
    float costheta = cos(theta), sintheta = sin(theta);
    float sinphi = sin(phi), cosphi = cos(phi);
    Direction<world_s> d = -LightToWorld(Direction<light_s>(sintheta * cosphi, sintheta * sinphi, costheta));
    *Ns = d;

    // Compute origin for infinite light sample ray
    Point<world_s> worldCenter;
    mfloat<length_d> worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Direction<world_s> v1, v2;
    CoordinateSystem(-d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point<world_s> Pdisk = worldCenter + (d1 * worldRadius * v1  + d2 * worldRadius * v2);
    *ray = Ray<world_s>(Pdisk + worldRadius * -d, d, mfloat<length_d>(0.0f), mfloat<length_d>(INFINITY), time);

    // Compute _InfiniteAreaLight_ ray PDF
    mfloat<invsolidangle_d> directionPdf = mfloat<invsolidangle_d>(mapPdf / (2.0f * M_PI * M_PI * sintheta));
    mfloat<invarea_d> areaPdf = mfloat<unit_d>(1.0f) / (M_PI * worldRadius * worldRadius);
    *pdf = directionPdf * areaPdf;
    if (sintheta == 0.f) *pdf =  mfloat<invsolidanglearea_d>(0.0f);
    return Spectrum<radiance_d>(radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
}


InfiniteAreaLight *CreateInfiniteLight(const Transform<light_s, world_s> &light2world,
                                       const ParamSet &paramSet) 
{
    Spectrum<unit_d> L = paramSet.FindOneSpectrum("L", Spectrum<unit_d>(1.0));
    Spectrum<unit_d> sc = paramSet.FindOneSpectrum("scale", Spectrum<unit_d>(1.0));
    string texmap = paramSet.FindOneString("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new InfiniteAreaLight(light2world, L * sc, nSamples, texmap);
}


