
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

#ifndef PBRT_CORE_DIFFGEOM_H
#define PBRT_CORE_DIFFGEOM_H

// core/diffgeom.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"

// DifferentialGeometry Declarations
class DifferentialGeometry {
public:
    DifferentialGeometry() { u = v = 0.0f; shape = NULL; }
    // DifferentialGeometry Public Methods
    DifferentialGeometry(const Point<world_s> &P, const Vector<world_s> &DPDU,
        const Vector<world_s> &DPDV, const Normal<world_s> &DNDU,
        const Normal<world_s> &DNDV, float uu, float vv,
        const Shape *sh);
    void ComputeDifferentials(const RayDifferential<world_s> &ray) const;

    // DifferentialGeometry Public Data
    Point<world_s> p;
    Normal<world_s> nn;
    float u, v;
    const Shape *shape;
    Vector<world_s> dpdu, dpdv;
    Normal<world_s> dndu, dndv;
    mutable Vector<world_s> dpdx, dpdy;
    mutable float dudx, dvdx, dudy, dvdy;
};



#endif // PBRT_CORE_DIFFGEOM_H
