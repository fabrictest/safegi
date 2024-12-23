
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

#ifndef PBRT_CORE_MGEOMETRY_H
#define PBRT_CORE_MGEOMETRY_H

// core/geometry.h*
#include "pbrt.h"
#include "uvector.h"

// Geometry Declarations
template<typename S>
class Vector : public mtuple3<length_d>{
public:
    // Vector Public Methods
    Vector() { x = y = z = mfloat<length_d>(); }
    Vector(const mfloat<length_d> &xx, const mfloat<length_d> &yy, const mfloat<length_d> &zz)
        : mtuple3<length_d>(xx, yy, zz) {
        Assert(!HasNaNs());
    }
    bool HasNaNs() const { return isnan(x.v) || isnan(y.v) || isnan(z.v); }
    explicit Vector(const Point<S> &p);
#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    Vector(const Vector &v) {
        Assert(!v.HasNaNs());
        x = v.x; y = v.y; z = v.z;
    }
    
    Vector &operator=(const Vector &v) {
        Assert(!v.HasNaNs());
        x = v.x; y = v.y; z = v.z;
        return *this;
    }
#endif // !NDEBUG
    inline Vector operator+(const Vector &v) const {
        Assert(!v.HasNaNs());
        return Vector(x + v.x, y + v.y, z + v.z);
    }
    
    inline Vector& operator+=(const Vector &v) {
        Assert(!v.HasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    inline Vector operator-(const Vector &v) const {
        Assert(!v.HasNaNs());
        return Vector(x - v.x, y - v.y, z - v.z);
    }
    
    inline Vector& operator-=(const Vector &v) {
        Assert(!v.HasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    inline bool operator==(const Vector &v) const {
        return x == v.x && y == v.y && z == v.z;
    }
    inline Vector operator*(float f) const { return Vector(f*x, f*y, f*z); }
    
    inline Vector &operator*=(float f) {
        Assert(!isnan(f));
        x *= f; y *= f; z *= f;
        return *this;
    }
    inline Vector operator/(float f) const {
        Assert(f != 0);
        float inv = 1.f / f;
        return Vector(x * inv, y * inv, z * inv);
    }
    
    inline Vector &operator/=(float f) {
        Assert(f != 0);
        float inv = 1.f / f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    inline Vector operator-() const { return Vector(-x, -y, -z); }
    inline mfloat<length_d> operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    
    inline mfloat<length_d> &operator[](int i) {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    inline mfloat<area_d> LengthSquared() const { return x*x + y*y + z*z; }
    inline mfloat<length_d> Length() const { return sqrt(LengthSquared()); }
};

template<typename S>
class Point : public mtuple3<length_d>{
public:
    // Point Public Methods
    Point() { x = y = z = mfloat<length_d>(); }
    Point(const mfloat<length_d> &xx, const mfloat<length_d> &yy, const mfloat<length_d> &zz)
        : mtuple3<length_d>(xx, yy, zz) {
        Assert(!HasNaNs());
    }
#ifndef NDEBUG
    Point(const Point &p) {
        Assert(!p.HasNaNs());
        x = p.x; y = p.y; z = p.z;
    }
    
    Point &operator=(const Point &p) {
        Assert(!p.HasNaNs());
        x = p.x; y = p.y; z = p.z;
        return *this;
    }
#endif // !NDEBUG
    inline Point<S> operator+(const Vector<S> &v) const {
        Assert(!v.HasNaNs());
        return Point(x + v.x, y + v.y, z + v.z);
    }
    
    inline Point<S> &operator+=(const Vector<S> &v) {
        Assert(!v.HasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    inline Vector<S> operator-(const Point<S> &p) const {
        Assert(!p.HasNaNs());
        return Vector<S>(x - p.x, y - p.y, z - p.z);
    }
    
    inline Point<S> operator-(const Vector<S> &v) const {
        Assert(!v.HasNaNs());
        return Point<S>(x - v.x, y - v.y, z - v.z);
    }
    
    inline Point<S> &operator-=(const Vector<S> &v) {
        Assert(!v.HasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    inline mfloat<length_d> operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    
    inline mfloat<length_d> &operator[](int i) {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    inline bool HasNaNs() const {
        return isnan(x.v) || isnan(y.v) || isnan(z.v);
    }
};

template <typename S> 
class Direction : public tuple3 {
public:
    // Normal Public Methods
    Direction() { x = 0.0f; y = 1.0f; z = 0.0f; }
    Direction(float xx, float yy, float zz)
        : tuple3(xx,yy ,zz) {
        Assert(!HasNaNs());
    }

    inline Direction operator-() const {
        return Direction(-x, -y, -z);
    }

    inline bool HasNaNs() const {
        return isnan(x) || isnan(y) || isnan(z);
    }

    template<typename _S>
    friend Vector<_S> operator*(const Direction<_S> &d, const mfloat<length_d> &f);
    template<typename _S>
    friend Vector<_S> operator*(const mfloat<length_d> &f, const Direction<_S> &d);
    
#ifndef NDEBUG
    Direction(const Direction &n) {
        Assert(!n.HasNaNs());
        x = n.x; y = n.y; z = n.z;
    }
    
    Direction &operator=(const Direction &n) {
        Assert(!n.HasNaNs());
        x = n.x; y = n.y; z = n.z;
        return *this;
    }
#endif // !NDEBUG
    inline float operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    
    inline float &operator[](int i) {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
};

template<typename S>
class Normal : public Direction<S>
{
public:
    Normal() : Direction<S>(){}
    Normal(float xx, float yy, float zz)
        : Direction<S>(xx, yy, zz) {}
    Normal(const Direction<S> &d) : Direction<S>(d) {}
};

template<typename S>
class Ray {
public:
    // Ray Public Methods
    Ray() : mint(0.f), maxt(INFINITY), time(0.f), depth(0) { }
    Ray(const Point<S> &origin, const Direction<S> &direction,
        mfloat<length_d> start, mfloat<length_d> end = mfloat<length_d>(INFINITY), mfloat<time_d> t = mfloat<time_d>(0.0f), int d = 0)
        : o(origin), d(direction), mint(start), maxt(end), time(t), depth(d) { }
    template<typename S_>
    Ray(const Point<S> &origin, const Direction<S> &direction, const Ray<S_> &parent,
        mfloat<length_d> start, mfloat<length_d> end = INFINITY)
        : o(origin), d(direction), mint(start), maxt(end),
          time(parent.time), depth(parent.depth+1) { }
    inline Point<S> operator()(const mfloat<length_d> &t) const { return o + d * t; }
    inline bool HasNaNs() const {
        return (o.HasNaNs() || d.HasNaNs() ||
                isnan(mint.v) || isnan(maxt.v));
    }

    // Ray Public Data
    Point<S> o;
    Direction<S> d;
    mutable mfloat<length_d> mint, maxt;
    mfloat<time_d> time;
    int depth;
};

template<typename S>
class RayDifferential : public Ray<S> {
public:
    // RayDifferential Public Methods
    RayDifferential() { hasDifferentials = false; }
    RayDifferential(const Point<S> &org, const Direction<S> &dir, mfloat<length_d> start,
        mfloat<length_d> end = mfloat<length_d>(INFINITY), mfloat<time_d> t = mfloat<time_d>(0.f), int d = 0)
            : Ray<S>(org, dir, start, end, t, d) {
        hasDifferentials = false;
    }
    template<typename S_>
    RayDifferential(const Point<S> &org, const Direction<S> &dir, const Ray<S_> &parent,
        mfloat<length_d> start, mfloat<length_d> end = mfloat<length_d>(INFINITY))
            : Ray<S>(org, dir, start, end, parent.time, parent.depth+1) {
        hasDifferentials = false;
    }
    inline explicit RayDifferential(const Ray<S> &ray) : Ray<S>(ray) {
        hasDifferentials = false;
    }
    inline bool HasNaNs() const {
        return Ray<S>::HasNaNs() ||
           (hasDifferentials && (rxOrigin.HasNaNs() || ryOrigin.HasNaNs() ||
                                 rxDirection.HasNaNs() || ryDirection.HasNaNs()));
    }

    // RayDifferential Public Data
    bool hasDifferentials;
    Point<S> rxOrigin, ryOrigin;
    Direction<S> rxDirection, ryDirection;
};

template<typename S>
class BBox {
public:
    // BBox Public Methods
    BBox() {
        pMin = Point<S>( mfloat<length_d>(INFINITY),  mfloat<length_d>(INFINITY),  mfloat<length_d>(INFINITY));
        pMax = Point<S>(mfloat<length_d>(-INFINITY), mfloat<length_d>(-INFINITY), mfloat<length_d>(-INFINITY));
    }
    BBox(const Point<S> &p) : pMin(p), pMax(p) { }
    BBox(const Point<S> &p1, const Point<S> &p2) {
        pMin = Point<S>(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z));
        pMax = Point<S>(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z));
    }

    template<typename _S>
    friend BBox<_S> Union(const BBox<_S> &b, const Point<_S> &p);
    template<typename _S>
    friend BBox<_S> Union(const BBox<_S> &b, const BBox<_S> &b2);

    inline bool Overlaps(const BBox &b) const {
        bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
        bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
        bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
        return (x && y && z);
    }
    inline bool Inside(const Point<S> &pt) const {
        return (pt.x >= pMin.x && pt.x <= pMax.x &&
                pt.y >= pMin.y && pt.y <= pMax.y &&
                pt.z >= pMin.z && pt.z <= pMax.z);
    }
    inline void Expand(float delta) {
        pMin -= Vector<S>(delta, delta, delta);
        pMax += Vector<S>(delta, delta, delta);
    }
    inline mfloat<area_d> SurfaceArea() const {
        Vector<S> d = pMax - pMin;
        return 2.f * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    inline mfloat<volume_d> Volume() const {
        Vector<S> d = pMax - pMin;
        return d.x * d.y * d.z;
    }
    inline int MaximumExtent() const {
        Vector<S> diag = pMax - pMin;
        if (diag.x > diag.y && diag.x > diag.z)
            return 0;
        else if (diag.y > diag.z)
            return 1;
        else
            return 2;
    }

    inline const Point<S> &operator[](int i) const {
        Assert(i == 0 || i == 1);
        return (&pMin)[i];
    }

    inline Point<S> &operator[](int i) {
        Assert(i == 0 || i == 1);
        return (&pMin)[i];
    }

    inline Point<S> Lerp(float tx, float ty, float tz) const {
        return Point<S>(::Lerp(tx, pMin.x, pMax.x), ::Lerp(ty, pMin.y, pMax.y),
                     ::Lerp(tz, pMin.z, pMax.z));
    }

    inline Vector<S> Offset(const Point<S> &p) const {
        return Vector<S>((p.x - pMin.x) / (pMax.x - pMin.x),
                      (p.y - pMin.y) / (pMax.y - pMin.y),
                      (p.z - pMin.z) / (pMax.z - pMin.z));
    }

    inline void BoundingSphere(Point<S> *c, mfloat<length_d> *rad) const 
    {
        *c = pMin + (pMax - pMin) * 0.5f;
        *rad = Inside(*c) ? Distance(*c, pMax) : mfloat<length_d>(0.0f);
    }

    bool IntersectP(const Ray<S> &ray, mfloat<length_d> *hitt0 = NULL, mfloat<length_d> *hitt1 = NULL) const {
            mfloat<length_d> t0 = ray.mint, t1 = ray.maxt;
            for (int i = 0; i < 3; ++i) {
                // Update interval for _i_th bounding box slab
                mfloat<length_d> tNear = (pMin[i] - ray.o[i]) / ray.d[i];
                mfloat<length_d> tFar  = (pMax[i] - ray.o[i]) / ray.d[i];

                // Update parametric interval from slab intersection $t$s
                if (tNear > tFar) swap(tNear, tFar);
                t0 = tNear > t0 ? tNear : t0;
                t1 = tFar  < t1 ? tFar  : t1;
                if (t0 > t1) return false;
            }
            if (hitt0) *hitt0 = t0;
            if (hitt1) *hitt1 = t1;
            return true;
    }

    // BBox Public Data
    Point<S> pMin, pMax;
};

template<typename S>
BBox<S> Union(const BBox<S> &b, const Point<S> &p) {
    BBox<S> ret = b;
    ret.pMin.x = min(b.pMin.x, p.x);
    ret.pMin.y = min(b.pMin.y, p.y);
    ret.pMin.z = min(b.pMin.z, p.z);
    ret.pMax.x = max(b.pMax.x, p.x);
    ret.pMax.y = max(b.pMax.y, p.y);
    ret.pMax.z = max(b.pMax.z, p.z);
    return ret;
}

template<typename S>
BBox<S> Union(const BBox<S> &b, const BBox<S> &b2) {
    BBox<S> ret;
    ret.pMin.x = min(b.pMin.x, b2.pMin.x);
    ret.pMin.y = min(b.pMin.y, b2.pMin.y);
    ret.pMin.z = min(b.pMin.z, b2.pMin.z);
    ret.pMax.x = max(b.pMax.x, b2.pMax.x);
    ret.pMax.y = max(b.pMax.y, b2.pMax.y);
    ret.pMax.z = max(b.pMax.z, b2.pMax.z);
    return ret;
}


// Geometry Inline Functions
template<typename S>
inline Vector<S>::Vector(const Point<S> &p)
    : mtuple3<length_d>(p.x, p.y, p.z) {
    Assert(!HasNaNs());
}

template<typename S>
inline Vector<S> AsVector(const Point<S> &p)
{
    Vector<S> v(p.x, p.y, p.z);
    Assert(!v.HasNaNs());
    return v;
}

template<typename S>
inline Vector<S> operator*(const Direction<S> &d, const mfloat<length_d> &f) 
{
    return Vector<S>(f*d.x, f*d.y, f*d.z);
}
template<typename S>
inline Vector<S> operator*(const mfloat<length_d> &f, const Direction<S> &d) 
{
    return Vector<S>(f*d.x, f*d.y, f*d.z);
}
template<typename S>
inline Vector<S> operator*(float f, const Vector<S> &v) { return v*f; }

template<typename S>
inline mfloat<length_d> Dot(const Vector<S> &n1, const Direction<S> &n2) {
    Assert(!n1.HasNaNs() && !n2.HasNaNs());
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

template<typename S>
inline mfloat<length_d> Dot(const Direction<S> &n1, const Vector<S> &n2) {
    Assert(!n1.HasNaNs() && !n2.HasNaNs());
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

template<typename S>
inline mfloat<area_d> Dot(const Vector<S> &v1, const Vector<S> &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template<typename S>
inline float Dot(const Direction<S> &n1, const Direction<S> &n2) {
    Assert(!n1.HasNaNs() && !n2.HasNaNs());
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

template<typename S>
inline mfloat<area_d> AbsDot(const Vector<S> &v1, const Vector<S> &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    return mfloat<area_d>(fabs(Dot(v1, v2).v));
}

template<typename S>
inline Vector<S> Cross(const Vector<S> &v1, const Direction<S> &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    return Vector<S>((v1.y * v2.z) - (v1.z * v2.y),
                  (v1.z * v2.x) - (v1.x * v2.z),
                  (v1.x * v2.y) - (v1.y * v2.x));
}

template<typename S>
inline Vector<S> Cross(const Direction<S> &v1, const Vector<S> &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    return Vector<S>((v1.y * v2.z) - (v1.z * v2.y),
                  (v1.z * v2.x) - (v1.x * v2.z),
                  (v1.x * v2.y) - (v1.y * v2.x));
}

template<typename S>
inline Vector<S> Cross(const Direction<S> &v1, const Direction<S> &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    return Vector<S>((v1.y * v2.z) - (v1.z * v2.y),
        (v1.z * v2.x) - (v1.x * v2.z),
        (v1.x * v2.y) - (v1.y * v2.x));
}

template<typename S>
inline Direction<S> Normalize(const Vector<S> &v) 
{
    mfloat<length_d> l = v.Length();
    return Direction<S>(v.x/l, v.y/l, v.z/l); 
}

template<typename S>
inline Direction<S> Bisector(const Direction<S>& a, const Direction<S>& b) {
    tuple3 r(a.x+b.x,a.y+b.y,a.z+b.z);
    float l = sqrt(r.x*r.x+r.y*r.y+r.z*r.z);
    return Direction<S>(r.x/l,r.y/l,r.z/l);
}

template<typename S>
inline Direction<S> Orthogonal(const Direction<S>& a, const Direction<S>& b) {
    tuple3 r(a.y*b.z-a.z*b.y, -a.x*b.z+a.z*b.x, a.x*b.y-a.y*b.x);
    float l = sqrt(r.x*r.x+r.y*r.y+r.z*r.z);
    return Direction<S>(r.x/l,r.y/l,r.z/l);
}

template<typename S>
inline Direction<S> Reflect(const Direction<S>& d, const Direction<S>& n) {
    float dot2 = 2.0f * Dot(d, n);
    tuple3 r(dot2 * n.x - d.x, 
        dot2 * n.y - d.y,
        dot2 * n.z - d.z);
    float l = sqrt(r.x*r.x+r.y*r.y+r.z*r.z);
    return Direction<S>(r.x/l,r.y/l,r.z/l);
}


template<typename S>
inline void CoordinateSystem(const Direction<S> &v1, Direction<S> *v2, Direction<S> *v3) {
    if (fabs(v1.x) > fabs(v1.y)) {
        float invLen = 1.f / sqrt(v1.x*v1.x + v1.z*v1.z);
        *v2 = Direction<S>(-v1.z * invLen, 0.f, v1.x * invLen);
    }
    else {
        float invLen = 1.f / sqrtf(v1.y*v1.y + v1.z*v1.z);
        *v2 = Direction<S>(0.f, v1.z * invLen, -v1.y * invLen);
    }
    *v3 = Orthogonal(v1, *v2);
}

template<typename S>
inline mfloat<length_d> Distance(const Point<S> &p1, const Point<S> &p2) {
    return (p1 - p2).Length();
}

template<typename S>
inline mfloat<area_d> DistanceSquared(const Point<S> &p1, const Point<S> &p2) {
    return (p1 - p2).LengthSquared();
}

template<typename S>
inline mfloat<proj_d> AbsDot(const Direction<S> &d1, const Direction<S> &d2) {
    Assert(!d1.HasNaNs() && !d2.HasNaNs());
    return mfloat<proj_d>(fabs(d1.x * d2.x + d1.y * d2.y + d1.z * d2.z));
}


template<typename S>
inline Normal<S> Faceforward(const Normal<S> &n, const Normal<S> &n2) {
    if(Dot(n, n2) < 0.f) 
        return -n;
    else
        return n;
}


template<typename S>
inline Vector<S> Faceforward(const Vector<S> &v, const Normal<S> &n2) {
    return (Dot(Normalize(v), n2) < 0.f) ? -v : v;
}

template<typename S>
inline Direction<S> Faceforward(const Direction<S> &v, const Direction<S> &n2) {
    return (Dot(v, n2) < 0.f) ? -v : v;
}

template<typename S>
inline Direction<S> SphericalDirection(const float sintheta,
                                 const float costheta, const mfloat<angle_d> &phi) {
    return Direction<S>(sintheta * cos(phi),
                  sintheta * sin(phi),
                  costheta);
}

template<typename S>
inline Direction<S> SphericalDirection(const float sintheta, const float costheta,
                                 const mfloat<angle_d> &phi, const Direction<S> &x,
                                 const Direction<S> &y, const Direction<S> &z) {
    return Normalize(sintheta * cos(phi) * x +
           sintheta * sin(phi) * y + costheta * z);
}

template<typename S>
inline mfloat<angle_d> SphericalTheta(const Direction<S> &v) {
    return _acos(Clamp(v.z, -1.f, 1.f));
}

template<typename S>
inline mfloat<angle_d> SphericalPhi(const Direction<S> &v) {
    mfloat<angle_d> p = _atan2(v.y, v.x);
    return (p < mfloat<angle_d>(0.0f)) ? p + circleAngle : p;
}



#endif // PBRT_CORE_MGEOMETRY_H
