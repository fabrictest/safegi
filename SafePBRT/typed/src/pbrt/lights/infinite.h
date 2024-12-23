
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

#ifndef PBRT_LIGHTS_INFINITE_H
#define PBRT_LIGHTS_INFINITE_H

// lights/infinite.h*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"

// InfiniteAreaLight Declarations
class InfiniteAreaLight : public Light {
public:
    // InfiniteAreaLight Public Methods
    InfiniteAreaLight(const Transform<light_s, world_s> &light2world, const Spectrum<unit_d> &power, int ns,
        const string &texmap);
    ~InfiniteAreaLight();
    Spectrum<power_d> Power(const Scene *) const;
    bool IsDeltaLight() const { return false; }
    Spectrum<radiance_d> Le(const RayDifferential<world_s> &r) const;
    Spectrum<radiance_d> Sample_L(const Point<world_s> &p, const mfloat<length_d> &pEpsilon, const LightSample &ls,
        const mfloat<time_d> &time, Direction<world_s> *wi, mfloat<invsolidangle_d> *pdf, VisibilityTester<world_s> *visibility) const;
    Spectrum<radiance_d> Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
        const mfloat<time_d> &time, Ray<world_s> *ray, Normal<world_s> *Ns, mfloat<invsolidanglearea_d> *pdf) const;
    mfloat<invsolidangle_d> Pdf(const Point<world_s> &, const Direction<world_s> &) const;
private:
    // InfiniteAreaLight Private Data
    MIPMap<RGBSpectrum<radiance_d> > *radianceMap;
    Distribution2D<radiance_d> *distribution;
};


InfiniteAreaLight *CreateInfiniteLight(const Transform<light_s, world_s> &light2world,
        const ParamSet &paramSet);

#endif // PBRT_LIGHTS_INFINITE_H
