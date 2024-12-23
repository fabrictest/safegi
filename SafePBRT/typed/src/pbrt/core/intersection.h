
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

#ifndef PBRT_CORE_INTERSECTION_H
#define PBRT_CORE_INTERSECTION_H

// core/intersection.h*
#include "pbrt.h"
#include "diffgeom.h"
#include "transform.h"

#include "shape.h"
#include "primitive.h"
#include "light.h"

// Intersection Declarations
class Intersection {
public:
    // Intersection Public Methods
    Intersection();
    BSDF *GetBSDF(const RayDifferential<world_s> &ray, MemoryArena &arena) const;
    BSSRDF *GetBSSRDF(const RayDifferential<world_s> &ray, MemoryArena &arena) const;
    Spectrum<radiance_d> Le(const Direction<world_s> &w) const;

    // Intersection Public Data
    DifferentialGeometry dg;
    const Primitive *primitive;
    Transform<world_s, object_s> WorldToObject;
    Transform<object_s, world_s> ObjectToWorld;
    uint32_t shapeId, primitiveId;
    mfloat<length_d> rayEpsilon;
private:
    // Intersection Private Methods
    Intersection(const Intersection &);
    Intersection &operator=(const Intersection &);
};



#endif // PBRT_CORE_INTERSECTION_H
