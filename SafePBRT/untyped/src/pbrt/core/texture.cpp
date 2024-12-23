
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


// core/texture.cpp*
#include "texture.h"
#include "shape.h"

// Texture Method Definitions
UVMapping2D::UVMapping2D(float ssu, float ssv, float ddu, float ddv)
    : su(ssu), sv(ssv), du(ddu), dv(ddv) { }
void UVMapping2D::Map(const DifferentialGeometry &dg,
                      float *s, float *t, float *dsdx, float *dtdx,
                      float *dsdy, float *dtdy) const {
    *s = su * dg.u + du;
    *t = sv * dg.v + dv;
    // Compute texture differentials for 2D identity mapping
    *dsdx = su * dg.dudx;
    *dtdx = sv * dg.dvdx;
    *dsdy = su * dg.dudy;
    *dtdy = sv * dg.dvdy;
}


void SphericalMapping2D::Map(const DifferentialGeometry &dg,
        float *s, float *t, float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const {
    sphere(dg.p, s, t);
    // Compute texture coordinate differentials for sphere $(u,v)$ mapping
    float sx, tx, sy, ty;
    const float delta = .1f;
    sphere(dg.p + delta * dg.dpdx, &sx, &tx);
    *dsdx = (sx - *s) / delta;
    *dtdx = (tx - *t) / delta;
    if (*dtdx > .5) *dtdx = 1.f - *dtdx;
    else if (*dtdx < -.5f) *dtdx = -(*dtdx + 1);
    sphere(dg.p + delta * dg.dpdy, &sy, &ty);
    *dsdy = (sy - *s) / delta;
    *dtdy = (ty - *t) / delta;
    if (*dtdy > .5) *dtdy = 1.f - *dtdy;
    else if (*dtdy < -.5f) *dtdy = -(*dtdy + 1);
}


void SphericalMapping2D::sphere(const Point &p, float *s, float *t) const {
    Vector vec = Normalize(WorldToTexture(p) - Point(0,0,0));
    float theta = SphericalTheta(vec);
    float phi = SphericalPhi(vec);
    *s = theta * INV_PI;
    *t = phi * INV_TWOPI;
}


void CylindricalMapping2D::Map(const DifferentialGeometry &dg,
        float *s, float *t, float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const {
    cylinder(dg.p, s, t);
    // Compute texture coordinate differentials for cylinder $(u,v)$ mapping
    float sx, tx, sy, ty;
    const float delta = .01f;
    cylinder(dg.p + delta * dg.dpdx, &sx, &tx);
    *dsdx = (sx - *s) / delta;
    *dtdx = (tx - *t) / delta;
    if (*dtdx > .5) *dtdx = 1.f - *dtdx;
    else if (*dtdx < -.5f) *dtdx = -(*dtdx + 1);
    cylinder(dg.p + delta * dg.dpdy, &sy, &ty);
    *dsdy = (sy - *s) / delta;
    *dtdy = (ty - *t) / delta;
    if (*dtdy > .5) *dtdy = 1.f - *dtdy;
    else if (*dtdy < -.5f) *dtdy = -(*dtdy + 1);
}


void PlanarMapping2D::Map(const DifferentialGeometry &dg,
        float *s, float *t, float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const {
    Vector vec = dg.p - Point(0,0,0);
    *s = ds + Dot(vec, vs);
    *t = dt + Dot(vec, vt);
    *dsdx = Dot(dg.dpdx, vs);
    *dtdx = Dot(dg.dpdx, vt);
    *dsdy = Dot(dg.dpdy, vs);
    *dtdy = Dot(dg.dpdy, vt);
}


Point IdentityMapping3D::Map(const DifferentialGeometry &dg,
                             Vector *dpdx, Vector *dpdy) const {
    *dpdx = WorldToTexture(dg.dpdx);
    *dpdy = WorldToTexture(dg.dpdy);
    return WorldToTexture(dg.p);
}



// Texture Function Definitions
float Lanczos(float x, float tau) {
    x = fabsf(x);
    if (x < 1e-5) return 1;
    if (x > 1.)    return 0;
    x *= M_PI;
    float s = sinf(x * tau) / (x * tau);
    float lanczos = sinf(x) / x;
    return s * lanczos;
}


