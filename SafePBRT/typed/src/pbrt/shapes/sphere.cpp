
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


// shapes/sphere.cpp*
#include "shapes/sphere.h"
#include "montecarlo.h"
#include "paramset.h"

// Sphere Method Definitions
Sphere::Sphere(const Transform<object_s, world_s> *o2w, 
    const Transform<world_s, object_s>  *w2o, bool ro,
    mfloat<length_d> rad, mfloat<length_d> z0, 
    mfloat<length_d> z1, float pm)
    : Shape(o2w, w2o, ro) {
        radius = rad;
        zmin = Clamp(min(z0, z1), -radius, radius);
        zmax = Clamp(max(z0, z1), -radius, radius);
        thetaMin = _acos(Clamp(zmin/radius, -1.f, 1.f));
        thetaMax = _acos(Clamp(zmax/radius, -1.f, 1.f));
        phiMax = Radians(Clamp(pm, 0.0f, 360.0f));
}


BBox<object_s> Sphere::ObjectBound() const {
    return BBox<object_s>(Point<object_s>(-radius, -radius, zmin),
        Point<object_s>( radius,  radius, zmax));
}


bool Sphere::Intersect(const Ray<world_s> &r, mfloat<length_d> *tHit, mfloat<length_d> *rayEpsilon,
    DifferentialGeometry *dg) const {
        mfloat<angle_d> phi;
        Point<object_s> phit;
        // Transform _Ray_ to object space
        Ray<object_s> ray;
        (*WorldToObject)(r, &ray);

        // Compute quadratic sphere coefficients
        auto A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
        auto B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z);
        auto C = ray.o.x*ray.o.x + ray.o.y*ray.o.y +
            ray.o.z*ray.o.z - radius*radius;

        // Solve quadratic equation for _t_ values
        mfloat<length_d> t0, t1;
        if (!Quadratic(A, B, C, &t0, &t1))
            return false;

        // Compute intersection distance along ray
        if (t0 > ray.maxt || t1 < ray.mint)
            return false;
        auto thit = t0;
        if (t0 < ray.mint) {
            thit = t1;
            if (thit > ray.maxt) return false;
        }

        // Compute sphere hit position and $\phi$
        phit = ray(thit);
        if (phit.x == mfloat<length_d>(0.0f) && phit.y == mfloat<length_d>(0.0f)) 
            phit.x = 1e-5 * radius;
        phi = _atan2(phit.y/radius, phit.x/radius);
        if (phi < mfloat<angle_d>(0.0f)) 
            phi += circleAngle;

        // Test sphere intersection against clipping parameters
        if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
            if (thit == t1) return false;
            if (t1 > ray.maxt) return false;
            thit = t1;
            // Compute sphere hit position and $\phi$
            phit = ray(thit);
            if (phit.x == mfloat<length_d>(0.0f) && phit.y == mfloat<length_d>(0.0f)) 
                phit.x = 1e-5 * radius;
            phi = _atan2(phit.y/radius, phit.x/radius);
            if (phi < mfloat<angle_d>(0.0f)) 
                phi += circleAngle;
            if (phit.z < zmin || phit.z > zmax || phi > phiMax)
                return false;
        }

        // Find parametric representation of sphere hit
        float u = phi / phiMax;
        mfloat<angle_d> theta = _acos(phit.z / radius);
        float v = (theta - thetaMin) / (thetaMax - thetaMin);

        // Compute sphere $\dpdu$ and $\dpdv$
        mfloat<length_d> zradius = sqrt(phit.x*phit.x + phit.y*phit.y);
        float cosphi = phit.x / zradius;
        float sinphi = phit.y / zradius;

        Vector<object_s> dpdu(__asfloat(-phiMax) * phit.y, __asfloat(phiMax) * phit.x, 0.0f*meters);
        Vector<object_s> dpdv = Vector<object_s> (phit.z * cosphi, phit.z * sinphi,
            -radius * sin(theta)) * __asfloat(thetaMax-thetaMin);

        // Compute sphere $\dndu$ and $\dndv$
        auto d2Pduu = 
            Vector<object_s>(phit.x, phit.y, 0.0f*meters) * __asfloat(-phiMax) * __asfloat(phiMax);
        auto d2Pduv = 
            Vector<object_s>(-sinphi*phit.z, cosphi*phit.z, 0.0f*phit.z) * __asfloat(thetaMax - thetaMin) * __asfloat(phiMax);
        auto d2Pdvv = 
            _Vector<length_d>(phit.x, phit.y, phit.z) * __asfloat(-(thetaMax - thetaMin)) * __asfloat(thetaMax - thetaMin);

        // Compute coefficients for fundamental forms
        auto E = Dot(dpdu, dpdu);
        auto F = Dot(dpdu, dpdv);
        auto G = Dot(dpdv, dpdv);
        auto N = Normalize<object_s>(Cross(dpdu, dpdv));
        auto e = Dot(N, d2Pduu);
        auto f = Dot(N, d2Pduv);
        auto g = Dot(N, d2Pdvv);

        // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
        auto invEGF2 = mfloat<unit_d>(1.0f) / (E*G - F*F);
        auto dndu = (f*F - e*G) * invEGF2 * dpdu + (e*F - f*E) * invEGF2 * dpdv;
        auto dndv = (g*F - f*G) * invEGF2 * dpdu + (f*F - g*E) * invEGF2 * dpdv;

        // Initialize _DifferentialGeometry_ from parametric information
        const auto &o2w = *ObjectToWorld;
        *dg = DifferentialGeometry(o2w(phit), o2w(dpdu), 
            o2w(dpdv), o2w(AsNormalDiff<object_s>(dndu)), 
            o2w(AsNormalDiff<object_s>(dndv)), u, v, this);

        // Update _tHit_ for quadric intersection
        *tHit = thit;

        // Compute _rayEpsilon_ for quadric intersection
        *rayEpsilon = 5e-4f * *tHit;
        return true;
}


bool Sphere::IntersectP(const Ray<world_s> &r) const {
    mfloat<angle_d> phi;
    Point<object_s> phit;
    // Transform _Ray_ to object space
    Ray<object_s> ray;
    (*WorldToObject)(r, &ray);

    // Compute quadratic sphere coefficients
    auto A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
    auto B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z);
    auto C = ray.o.x*ray.o.x + ray.o.y*ray.o.y +
        ray.o.z*ray.o.z - radius*radius;

    // Solve quadratic equation for _t_ values
    mfloat<length_d> t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    // Compute intersection distance along ray
    if (t0 > ray.maxt || t1 < ray.mint)
        return false;
    auto thit = t0;
    if (t0 < ray.mint) {
        thit = t1;
        if (thit > ray.maxt) return false;
    }

    // Compute sphere hit position and $\phi$
    phit = ray(thit);
    if (phit.x == mfloat<length_d>(0.0f) && phit.y == mfloat<length_d>(0.0f)) 
        phit.x = 1e-5 * radius;
    phi = _atan2(phit.y/radius, phit.x/radius);
    if (phi < mfloat<angle_d>(0.0f) ) 
        phi += circleAngle;

    // Test sphere intersection against clipping parameters
    if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
        if (thit == t1) return false;
        if (t1 > ray.maxt) return false;
        thit = t1;
        // Compute sphere hit position and $\phi$
        phit = ray(thit);
        if (phit.x == mfloat<length_d>(0.0f) && phit.y == mfloat<length_d>(0.0f) ) phit.x = 1e-5 * radius;
        phi = _atan2(phit.y/radius, phit.x/radius);
        if (phi < mfloat<angle_d>(0.0f)) phi += circleAngle;
        if (phit.z < zmin || phit.z > zmax || phi > phiMax)
            return false;
    }
    return true;
}


mfloat<area_d> Sphere::Area() const {
    return (phiMax/circleAngle) * (zmax-zmin) * (2.0f * M_PI * radius);
}


Sphere *CreateSphereShape(const Transform<object_s, world_s> *o2w, const Transform<world_s, object_s> *w2o,
    bool reverseOrientation, const ParamSet &params) {
        mfloat<length_d> radius = params.FindOneFloat("radius", 1.0f*meters);
        mfloat<length_d> zmin = params.FindOneFloat("zmin", -radius);
        mfloat<length_d> zmax = params.FindOneFloat("zmax", radius);
        float phimax = params.FindOneFloat("phimax", 360.f);
        return new Sphere(o2w, w2o, reverseOrientation, radius,
            zmin, zmax, phimax);
}


Point<world_s> Sphere::Sample(float u1, float u2, Normal<world_s> *ns) const {
    Point<object_s> p = Point<object_s>() + radius * UniformSampleSphere<object_s>(u1, u2);
    *ns = (*ObjectToWorld)(Normalize(AsVector(p)));
    if (ReverseOrientation) *ns = -*ns;
    return (*ObjectToWorld)(p);
}


Point<world_s> Sphere::Sample(const Point<world_s> &p, float u1, float u2,
    Normal<world_s> *ns) const {
        // Compute coordinate system for sphere sampling
        Point<world_s> Pcenter = (*ObjectToWorld)(Point<object_s>());
        Direction<world_s> wc = Normalize(Pcenter - p);
        Direction<world_s> wcX, wcY;
        CoordinateSystem(wc, &wcX, &wcY);

        // Sample uniformly on sphere if $\pt{}$ is inside it
        if (DistanceSquared(p, Pcenter) - radius*radius < 1e-4f*sqr_meters)
            return Sample(u1, u2, ns);

        // Sample sphere uniformly inside subtended cone
        float sinThetaMax2 = radius*radius / DistanceSquared(p, Pcenter);
        float cosThetaMax = sqrtf(max(0.f, 1.f - sinThetaMax2));
        DifferentialGeometry dgSphere;
        mfloat<length_d> thit, rayEpsilon;
        Point<world_s> ps;
        Ray<world_s> r(p, UniformSampleCone(u1, u2, cosThetaMax, wcX, wcY, wc), 1e-3f*meters);
        if (!Intersect(r, &thit, &rayEpsilon, &dgSphere))
            ps = Pcenter - radius * wc;
        else
            ps = r(thit);
        *ns = Normal<world_s>(Normalize(ps - Pcenter));
        if (ReverseOrientation) *ns = -*ns;
        return ps;
}


mfloat<invsolidangle_d> Sphere::Pdf(const Point<world_s> &p, 
    const Direction<world_s> &wi) const {
    Point<world_s> Pcenter = (*ObjectToWorld)(Point<object_s>());
    // Return uniform weight if point inside sphere
    if (DistanceSquared(p, Pcenter) - radius*radius < 1e-4f*sqr_meters)
        return Shape::Pdf(p, wi);

    // Compute general sphere weight
    float sinThetaMax2 = radius*radius / DistanceSquared(p, Pcenter);
    float cosThetaMax = sqrt(max(0.f, 1.f - sinThetaMax2));
    return UniformConePdf(cosThetaMax);
}


