
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


// lights/point.cpp*
#include "lights/point.h"
#include "scene.h"
#include "paramset.h"
#include "montecarlo.h"

// PointLight Method Definitions
PointLight::PointLight(const Transform<light_s, world_s> &light2world,
                       const Spectrum<intensity_d> &intensity)
    : Light(light2world) {
    lightPos = LightToWorld(Point<light_s>());
    Intensity = intensity;
}


Spectrum<radiance_d> PointLight::Sample_L(const Point<world_s> &p, const mfloat<length_d> &pEpsilon,
         const LightSample &ls, const mfloat<time_d> &time, Direction<world_s> *wi, mfloat<invsolidangle_d> *pdf,
         VisibilityTester<world_s> *visibility) const {
    *wi = Normalize(lightPos - p);
    *pdf = mfloat<invsolidangle_d>(1.0f);
    visibility->SetSegment(p, pEpsilon, lightPos, 0., time);
    return Intensity / (DistanceSquared(lightPos, p) * mfloat<proj_d>(1.0f));
}


Spectrum<power_d> PointLight::Power(const Scene *) const {
    return Intensity * sphericalAngle;
}


PointLight *CreatePointLight(const Transform<light_s, world_s> &light2world,
        const ParamSet &paramSet) {
    Spectrum<intensity_d> I = paramSet.FindOneSpectrum("I", Spectrum<intensity_d>(mfloat<intensity_d>(1.0f)));
    Spectrum<rho_d> sc = paramSet.FindOneSpectrum("scale", Spectrum<rho_d>(mfloat<rho_d>(1.0f)));
    Point<world_s> P = paramSet.FindOnePoint("from", Point<world_s>());
    Transform<light_s, world_s> l2w = Translate<world_s, world_s>(AsVector(P)) * light2world;
    return new PointLight(l2w, I * sc);
}


mfloat<invsolidangle_d> PointLight::Pdf(const Point<world_s> &, const Direction<world_s> &) const {
    return mfloat<invsolidangle_d>(0.0f);
}


Spectrum<radiance_d> PointLight::Sample_L(const Scene *scene, const LightSample &ls,
        float u1, float u2, const mfloat<time_d> &time, Ray<world_s> *ray, Normal<world_s> *Ns,
        mfloat<invsolidanglearea_d> *pdf) const {
    *ray = Ray<world_s>(lightPos, UniformSampleSphere<world_s>(ls.uPos[0], ls.uPos[1]),
               0.f*meters, INFINITY*meters, time);
    *Ns = ray->d;
    *pdf = UniformSpherePdf() / mfloat<area_d>(1.0f);
    return Intensity * (mfloat<invprojarea_d>(1.0f));
}
