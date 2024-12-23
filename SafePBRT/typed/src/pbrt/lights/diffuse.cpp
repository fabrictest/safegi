
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


// lights/diffuse.cpp*
#include "lights/diffuse.h"
#include "paramset.h"
#include "montecarlo.h"

// DiffuseAreaLight Method Definitions
DiffuseAreaLight::~DiffuseAreaLight() {
    delete shapeSet;
}

DiffuseAreaLight::DiffuseAreaLight(const Transform<light_s, world_s> &light2world,
    const Spectrum<radiance_d> &le, int ns, const Reference<Shape> &s)
    : AreaLight(light2world, ns) 
{
    Lemit = le;
    shapeSet = new ShapeSet(s);
    area = shapeSet->Area();
}


Spectrum<power_d> DiffuseAreaLight::Power(const Scene *) const 
{
    return Lemit * area * hemisphericalProjectedAngle;
}

Spectrum<radiance_d> DiffuseAreaLight::Sample_L(const Point<world_s> &p, const mfloat<length_d> &pEpsilon,
    const LightSample &ls, const mfloat<time_d> &time, Direction<world_s> *wi, 
    mfloat<invsolidangle_d> *pdf, VisibilityTester<world_s> *visibility) const 
{
    Normal<world_s> ns;
    Point<world_s> ps = shapeSet->Sample(p, ls, &ns);
    *wi = Normalize(ps - p);
    *pdf = shapeSet->Pdf(p, *wi);
    visibility->SetSegment(p, pEpsilon, ps, 1e-3f, time);
    return L(ps, ns, -*wi);
}


mfloat<invsolidangle_d> DiffuseAreaLight::Pdf(const Point<world_s> &p, const Direction<world_s> &wi) const {
    return shapeSet->Pdf(p, wi);
}


Spectrum<radiance_d> DiffuseAreaLight::Sample_L(const Scene *scene,
    const LightSample &ls, float u1, float u2, const mfloat<time_d> &time,
    Ray<world_s> *ray, Normal<world_s> *Ns, mfloat<invsolidanglearea_d> *pdf) const 
{
    Point<world_s> org = shapeSet->Sample(ls, Ns);
    Direction<world_s> dir = UniformSampleSphere<world_s>(u1, u2);
    if (Dot(dir, *Ns) < 0.0f) dir = -dir;
    *ray = Ray<world_s>(org, dir, mfloat<length_d>(1e-3f), mfloat<length_d>(INFINITY), time);
    *pdf = shapeSet->Pdf(org) / hemisphericalAngle;
    return L(org, *Ns, dir);
}


AreaLight *CreateDiffuseAreaLight(const Transform<light_s, world_s> &light2world, const ParamSet &paramSet,
    const Reference<Shape> &shape) 
{
    Spectrum<radiance_d> L = paramSet.FindOneSpectrum("L", Spectrum<radiance_d>(mfloat<radiance_d>(1.0f)));
    Spectrum<rho_d> sc = paramSet.FindOneSpectrum("scale", Spectrum<rho_d>(1.0f));
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new DiffuseAreaLight(light2world, L * sc, nSamples, shape);
}

