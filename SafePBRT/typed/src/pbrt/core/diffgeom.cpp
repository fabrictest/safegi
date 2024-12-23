
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


// core/diffgeom.cpp*
#include "diffgeom.h"
#include "transform.h"
#include "shape.h"


DifferentialGeometry::DifferentialGeometry(const Point<world_s> &P, const Vector<world_s> &DPDU,
                     const Vector<world_s> &DPDV, const Normal<world_s> &DNDU,
                     const Normal<world_s> &DNDV, float uu, float vv,
                     const Shape *sh): p(P), dpdu(DPDU), dpdv(DPDV), dndu(DNDU), dndv(DNDV) 
{
    // Initialize _DifferentialGeometry_ from parameters
    nn = Normalize<world_s>(Cross(dpdu, dpdv));
    u = uu;
    v = vv;
    shape = sh;
    dudx = dvdx = dudy = dvdy = 0;

    // Adjust normal based on orientation and handedness
    if (shape && (shape->ReverseOrientation ^ shape->TransformSwapsHandedness))
        nn = -nn;
}

void DifferentialGeometry::ComputeDifferentials(const RayDifferential<world_s> &ray) const
{
    if (ray.hasDifferentials) {
        // Estimate screen space change in $\pt{}$ and $(u,v)$

        // Compute auxiliary intersection points with plane
        auto d = -Dot(nn, Vector<world_s>(p.x, p.y, p.z));
        Vector<world_s> rxv(ray.rxOrigin.x, ray.rxOrigin.y, ray.rxOrigin.z);
        auto tx = -(Dot(nn, rxv) + d) / Dot(nn, ray.rxDirection);
        auto px = ray.rxOrigin + tx * ray.rxDirection;
        Vector<world_s> ryv(ray.ryOrigin.x, ray.ryOrigin.y, ray.ryOrigin.z);
        auto ty = -(Dot(nn, ryv) + d) / Dot(nn, ray.ryDirection);
        auto py = ray.ryOrigin + ty * ray.ryDirection;
        dpdx = px - p;
        dpdy = py - p;

        // Compute $(u,v)$ offsets at auxiliary points

        // Initialize _A_, _Bx_, and _By_ matrices for offset computation
        mfloat<length_d> A[2][2], Bx[2], By[2];
        int axes[2];
        if (fabsf(nn.x) > fabsf(nn.y) && fabsf(nn.x) > fabsf(nn.z)) {
            axes[0] = 1; axes[1] = 2;
        }
        else if (fabsf(nn.y) > fabsf(nn.z)) {
            axes[0] = 0; axes[1] = 2;
        }
        else {
            axes[0] = 0; axes[1] = 1;
        }

        // Initialize matrices for chosen projection plane
        A[0][0] = dpdu[axes[0]];
        A[0][1] = dpdv[axes[0]];
        A[1][0] = dpdu[axes[1]];
        A[1][1] = dpdv[axes[1]];
        Bx[0] = px[axes[0]] - p[axes[0]];
        Bx[1] = px[axes[1]] - p[axes[1]];
        By[0] = py[axes[0]] - p[axes[0]];
        By[1] = py[axes[1]] - p[axes[1]];
        if (!SolveLinearSystem2x2(A, Bx, &dudx, &dvdx)) {
            dudx = 1.; dvdx = 0.;
        }
        if (!SolveLinearSystem2x2(A, By, &dudy, &dvdy)) {
            dudy = 0.; dvdy = 1.;
        }
    }
    else {
        dudx = dvdx = 0.;
        dudy = dvdy = 0.;
        dpdx = dpdy = Vector<world_s>();
    }
}