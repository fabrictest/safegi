
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


// core/intersection.cpp*
#include "intersection.h"

Intersection::Intersection() {
    primitive = NULL;
    shapeId = primitiveId = 0;
    rayEpsilon = mfloat<length_d>(0.0f);
}

BSDF* Intersection::GetBSDF(const RayDifferential<world_s> &ray, MemoryArena &arena) const
{
    PBRT_STARTED_BSDF_SHADING(const_cast<RayDifferential *>(&ray));
    dg.ComputeDifferentials(ray);
    BSDF *bsdf = primitive->GetBSDF(dg, ObjectToWorld, arena);
    PBRT_FINISHED_BSDF_SHADING(const_cast<RayDifferential *>(&ray), bsdf);
    return bsdf;
}

BSSRDF* Intersection::GetBSSRDF(const RayDifferential<world_s> &ray, MemoryArena &arena) const
{
    PBRT_STARTED_BSSRDF_SHADING(const_cast<RayDifferential *>(&ray));
    dg.ComputeDifferentials(ray);
    BSSRDF *bssrdf = primitive->GetBSSRDF(dg, ObjectToWorld, arena);
    PBRT_FINISHED_BSSRDF_SHADING(const_cast<RayDifferential *>(&ray), bssrdf);
    return bssrdf;
}

Spectrum<radiance_d> Intersection::Le(const Direction<world_s> &w) const
{
    const AreaLight *area = primitive->GetAreaLight();
    return area ? area->L(dg.p, dg.nn, w) : Spectrum<radiance_d>();
}