
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

#ifndef PBRT_LIGHTS_POINT_H
#define PBRT_LIGHTS_POINT_H

// lights/point.h*
#include "pbrt.h"
#include "light.h"
#include "shape.h"

// PointLight Declarations
class PointLight : public Light {
public:
    // PointLight Public Methods
    PointLight(const Transform<light_s, world_s> &light2world, 
        const Spectrum<intensity_d> &intensity);
    Spectrum<radiance_d> Sample_L(const Point<world_s> &p, const mfloat<length_d> &pEpsilon, 
        const LightSample &ls, const mfloat<time_d> &time, Direction<world_s> *wi, 
        mfloat<invsolidangle_d> *pdf, VisibilityTester<world_s> *vis) const;
    Spectrum<power_d> Power(const Scene *) const;
    bool IsDeltaLight() const { return true; }
    Spectrum<radiance_d> Sample_L(const Scene *scene, const LightSample &ls, 
        float u1, float u2, const mfloat<time_d> &time,
        Ray<world_s> *ray, Normal<world_s> *Ns, mfloat<invsolidanglearea_d> *pdf) const;
    mfloat<invsolidangle_d> Pdf(const Point<world_s> &, const Direction<world_s> &) const;
private:
    // PointLight Private Data
    Point<world_s> lightPos;
    Spectrum<intensity_d> Intensity;
};


PointLight *CreatePointLight(const Transform<light_s, world_s> &light2world,
        const ParamSet &paramSet);

#endif // PBRT_LIGHTS_POINT_H
