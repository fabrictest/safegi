
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


// shapes/disk.cpp*
#include "shapes/disk.h"
#include "paramset.h"
#include "montecarlo.h"

// Disk Method Definitions
Disk::Disk(const Transform<object_s, world_s> *o2w, const Transform<world_s, object_s> *w2o, 
    bool ro, const mfloat<length_d> &ht, const mfloat<length_d> &r, 
    const mfloat<length_d> &ri, float tmax) : Shape(o2w, w2o, ro) {
        height = ht;
        radius = r;
        innerRadius = ri;
        phiMax = Radians(Clamp(tmax, 0.0f, 360.0f));
}


BBox<object_s> Disk::ObjectBound() const {
    return BBox<object_s>(Point<object_s>(-radius, -radius, height),
        Point<object_s>( radius,  radius, height));
}


bool Disk::Intersect(const Ray<world_s> &r, mfloat<length_d> *tHit, mfloat<length_d> *rayEpsilon,
    DifferentialGeometry *dg) const {
        // Transform _Ray_ to object space
        Ray<object_s> ray;
        (*WorldToObject)(r, &ray);

        // Compute plane intersection for disk
        if (fabs(ray.d.z) < (1e-7)) return false;
        auto thit = (height - ray.o.z) / ray.d.z;
        if (thit < ray.mint || thit > ray.maxt)
            return false;

        // See if hit point is inside disk radii and $\phimax$
        auto phit = ray(thit);
        auto dist2 = phit.x * phit.x + phit.y * phit.y;
        if (dist2 > radius * radius || dist2 < innerRadius * innerRadius)
            return false;

        // Test disk $\phi$ value against $\phimax$
        auto phi = _atan2(phit.y/radius, phit.x/radius);
        if (phi < mfloat<angle_d>(0.0f)) phi += circleAngle;
        if (phi > phiMax)
            return false;

        // Find parametric representation of disk hit
        float u = phi / phiMax;
        float v = 1.f - ((sqrt(dist2)-innerRadius) / (radius-innerRadius));
        _Vector<length_d> dpdu(__asfloat(-phiMax) * phit.y, __asfloat(phiMax) * phit.x, mfloat<length_d>(0.0f));
        _Vector<length_d> dpdv(-phit.x / (1-v), -phit.y / (1-v), mfloat<length_d>(0.0f));
        dpdu *= phiMax / circleAngle;
        dpdv *= (radius - innerRadius) / radius;
        Normal<object_s> dndu(0,0,0), dndv(0,0,0);

        // Initialize _DifferentialGeometry_ from parametric information
        const Transform<object_s, world_s> &o2w = *ObjectToWorld;
        *dg = DifferentialGeometry(o2w(phit), o2w(AsVector<object_s>(dpdu)), o2w(AsVector<object_s>(dpdv)),
            o2w(dndu), o2w(dndv), u, v, this);

        // Update _tHit_ for quadric intersection
        *tHit = thit;

        // Compute _rayEpsilon_ for quadric intersection
        *rayEpsilon = 5e-4f * *tHit;
        return true;
}


bool Disk::IntersectP(const Ray<world_s> &r) const {
    // Transform _Ray_ to object space
    Ray<object_s> ray;
    (*WorldToObject)(r, &ray);

    // Compute plane intersection for disk
    if (fabsf(ray.d.z) < 1e-7) return false;
    mfloat<length_d> thit = (height - ray.o.z) / ray.d.z;
    if (thit < ray.mint || thit > ray.maxt)
        return false;

    // See if hit point is inside disk radii and $\phimax$
    Point<object_s> phit = ray(thit);
    mfloat<area_d> dist2 = phit.x * phit.x + phit.y * phit.y;
    if (dist2 > radius * radius || dist2 < innerRadius * innerRadius)
        return false;

    // Test disk $\phi$ value against $\phimax$
    mfloat<angle_d> phi = _atan2(phit.y/radius, phit.x/radius);
    if (phi < mfloat<angle_d>(0.0f)) phi += circleAngle;
    if (phi > phiMax)
        return false;
    return true;
}


mfloat<area_d> Disk::Area() const {
    return phiMax/circleAngle * M_PI * (radius * radius - innerRadius * innerRadius);
}


Disk *CreateDiskShape(const Transform<object_s, world_s> *o2w, const Transform<world_s, object_s> *w2o,
    bool reverseOrientation, const ParamSet &params) {
        mfloat<length_d> height = params.FindOneFloat<length_d>("height", mfloat<length_d>(0.0f));
        mfloat<length_d> radius = params.FindOneFloat<length_d>("radius", mfloat<length_d>(1.0f));
        mfloat<length_d> inner_radius = params.FindOneFloat<length_d>("innerradius", mfloat<length_d>(0.0f));
        float phimax = params.FindOneFloat("phimax", 360);
        return new Disk(o2w, w2o, reverseOrientation, height, radius, inner_radius, phimax);
}


Point<world_s> Disk::Sample(float u1, float u2, Normal<world_s> *Ns) const {
    Point<object_s> p;
    float x; 
    float y;
    ConcentricSampleDisk(u1, u2, &x, &y);
    p.x = x * radius;
    p.y = y * radius;
    p.z = height;
    *Ns = (*ObjectToWorld)(Normal<object_s>(0,0,1));
    if (ReverseOrientation) *Ns = -*Ns ;
    return (*ObjectToWorld)(p);
}


