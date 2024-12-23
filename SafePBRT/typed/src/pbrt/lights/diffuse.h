
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

#ifndef PBRT_LIGHTS_DIFFUSE_H
#define PBRT_LIGHTS_DIFFUSE_H

// lights/diffuse.h*
#include "pbrt.h"
#include "light.h"
#include "primitive.h"

// DiffuseAreaLight Declarations
class DiffuseAreaLight : public AreaLight {
public:
    // DiffuseAreaLight Public Methods
    DiffuseAreaLight(const Transform<light_s, world_s> &light2world,
        const Spectrum<radiance_d> &Le, int ns, const Reference<Shape> &shape);
    ~DiffuseAreaLight();
    Spectrum<radiance_d> L(const Point<world_s> &p, const Normal<world_s> &n, const Direction<world_s> &w) const {
        return Dot(n, w) > 0.f ? Lemit : Spectrum<radiance_d>();
    }
    Spectrum<power_d> Power(const Scene *) const;
    bool IsDeltaLight() const { return false; }
    mfloat<invsolidangle_d> Pdf(const Point<world_s> &, const Direction<world_s> &) const;
    Spectrum<radiance_d> Sample_L(const Point<world_s> &P, const mfloat<length_d> &pEpsilon, const LightSample &ls, const mfloat<time_d> &time,
        Direction<world_s> *wo, mfloat<invsolidangle_d> *pdf, VisibilityTester<world_s> *visibility) const;
    Spectrum<radiance_d> Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
        const mfloat<time_d> &time, Ray<world_s> *ray, Normal<world_s> *Ns, mfloat<invsolidanglearea_d> *pdf) const;
protected:
    // DiffuseAreaLight Protected Data
    Spectrum<radiance_d> Lemit;
    ShapeSet *shapeSet;
    mfloat<area_d> area;
};


AreaLight *CreateDiffuseAreaLight(const Transform<light_s, world_s> &light2world, const ParamSet &paramSet,
        const Reference<Shape> &shape);

#endif // PBRT_LIGHTS_DIFFUSE_H
