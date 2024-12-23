
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

#ifndef PBRT_CORE_TRANSFORM_H
#define PBRT_CORE_TRANSFORM_H

// core/transform.h*
#include "pbrt.h"
#include "geometry.h"

// Matrix4x4 Declarations
struct Matrix4x4 {
    // Matrix4x4 Public Methods
    Matrix4x4() {
        m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.f;
        m[0][1] = m[0][2] = m[0][3] = m[1][0] =
             m[1][2] = m[1][3] = m[2][0] = m[2][1] = m[2][3] =
             m[3][0] = m[3][1] = m[3][2] = 0.f;
    }
    Matrix4x4(float mat[4][4]);
    Matrix4x4(float t00, float t01, float t02, float t03,
              float t10, float t11, float t12, float t13,
              float t20, float t21, float t22, float t23,
              float t30, float t31, float t32, float t33);
    bool operator==(const Matrix4x4 &m2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (m[i][j] != m2.m[i][j]) return false;
        return true;
    }
    bool operator!=(const Matrix4x4 &m2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (m[i][j] != m2.m[i][j]) return true;
        return false;
    }
    friend Matrix4x4 Transpose(const Matrix4x4 &);
    void Print(FILE *f) const {
        fprintf(f, "[ ");
        for (int i = 0; i < 4; ++i) {
            fprintf(f, "  [ ");
            for (int j = 0; j < 4; ++j)  {
                fprintf(f, "%f", m[i][j]);
                if (j != 3) fprintf(f, ", ");
            }
            fprintf(f, " ]\n");
        }
        fprintf(f, " ] ");
    }
    static Matrix4x4 Mul(const Matrix4x4 &m1, const Matrix4x4 &m2) {
        Matrix4x4 r;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                r.m[i][j] = m1.m[i][0] * m2.m[0][j] +
                            m1.m[i][1] * m2.m[1][j] +
                            m1.m[i][2] * m2.m[2][j] +
                            m1.m[i][3] * m2.m[3][j];
        return r;
    }
    friend Matrix4x4 Inverse(const Matrix4x4 &);
    float m[4][4];
};



// Transform Declarations
template<typename S1, typename S2>
class Transform {
public:
    // Transform Public Methods
    Transform() { }
    Transform(const float mat[4][4]) {
        m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3],
                      mat[1][0], mat[1][1], mat[1][2], mat[1][3],
                      mat[2][0], mat[2][1], mat[2][2], mat[2][3],
                      mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
        mInv = Inverse(m);
    }
    Transform(const Matrix4x4 &mat)
        : m(mat), mInv(Inverse(mat)) {
    }
    Transform(const Matrix4x4 &mat, const Matrix4x4 &minv)
       : m(mat), mInv(minv) {
    }

    template<typename S1_, typename S2_>    
    friend Transform<S2_, S1_> Inverse(const Transform<S1_, S2_> &t);

    bool operator==(const Transform &t) const {
        return t.m == m && t.mInv == mInv;
    }
    bool operator!=(const Transform &t) const {
        return t.m != m || t.mInv != mInv;
    }
    bool operator<(const Transform &t2) const {
        for (uint32_t i = 0; i < 4; ++i)
            for (uint32_t j = 0; j < 4; ++j) {
                if (m.m[i][j] < t2.m.m[i][j]) return true;
                if (m.m[i][j] > t2.m.m[i][j]) return false;
            }
        return false;
    }
    bool IsIdentity() const {
        return (m.m[0][0] == 1.f && m.m[0][1] == 0.f &&
                m.m[0][2] == 0.f && m.m[0][3] == 0.f &&
                m.m[1][0] == 0.f && m.m[1][1] == 1.f &&
                m.m[1][2] == 0.f && m.m[1][3] == 0.f &&
                m.m[2][0] == 0.f && m.m[2][1] == 0.f &&
                m.m[2][2] == 1.f && m.m[2][3] == 0.f &&
                m.m[3][0] == 0.f && m.m[3][1] == 0.f &&
                m.m[3][2] == 0.f && m.m[3][3] == 1.f);
    }
    const Matrix4x4 &GetMatrix() const { return m; }
    const Matrix4x4 &GetInverseMatrix() const { return mInv; }
    bool HasScale() const {
        mfloat<area_d> la2 = (*this)(Vector<S1>(mfloat<length_d>(1.0f),mfloat<length_d>(0.0f),mfloat<length_d>(0.0f))).LengthSquared();
        mfloat<area_d> lb2 = (*this)(Vector<S1>(mfloat<length_d>(0.0f),mfloat<length_d>(1.0f),mfloat<length_d>(0.0f))).LengthSquared();
        mfloat<area_d> lc2 = (*this)(Vector<S1>(mfloat<length_d>(0.0f),mfloat<length_d>(0.0f),mfloat<length_d>(1.0f))).LengthSquared();
#define NOT_ONE(x) ((x) < mfloat<area_d>(.999f) || (x) > mfloat<area_d>(1.001f))
        return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
    }

    inline Point<S2> operator()(const Point<S1> &pt) const {
        mfloat<length_d> x = pt.x, y = pt.y, z = pt.z;
        mfloat<length_d> xp = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z + mfloat<length_d>(m.m[0][3]);
        mfloat<length_d> yp = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z + mfloat<length_d>(m.m[1][3]);
        mfloat<length_d> zp = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z + mfloat<length_d>(m.m[2][3]);
        float wp = __asfloat(m.m[3][0]*x + m.m[3][1]*y + m.m[3][2]*z + mfloat<length_d>(m.m[3][3]));
        Assert(wp != 0);
        if (wp == 1.0) return Point<S2>(xp, yp, zp);
        else          return Point<S2>(xp/wp, yp/wp, zp/wp);
    }

    inline void operator()(const Point<S1> &pt, Point<S2> *ptrans) const {
            mfloat<length_d> x = pt.x, y = pt.y, z = pt.z;
            ptrans->x = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z + mfloat<length_d>(m.m[0][3]);
            ptrans->y = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z + mfloat<length_d>(m.m[1][3]);
            ptrans->z = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z + mfloat<length_d>(m.m[2][3]);
            float w   = __asfloat(m.m[3][0]*x + m.m[3][1]*y + m.m[3][2]*z + mfloat<length_d>(m.m[3][3]));
            Assert(w != 0);
            if (w != 1.0f) *ptrans = Point<S2>(ptrans->x/w, ptrans->y/w, ptrans->z/w);
    }

    inline Vector<S2> operator()(const Vector<S1> &v) const {
        mfloat<length_d> x = v.x, y = v.y, z = v.z;
        return Vector<S2>(m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z,
            m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z,
            m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z);
    }

    inline void operator()(const Vector<S1> &v,
        Vector<S2> *vt) const {
            mfloat<length_d> x = v.x, y = v.y, z = v.z;
            vt->x = m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z;
            vt->y = m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z;
            vt->z = m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z;
    }

    inline Direction<S2> operator()(const Direction<S1> &n) const {
        float x = n.x, y = n.y, z = n.z;
        return Direction<S2>(mInv.m[0][0]*x + mInv.m[1][0]*y + mInv.m[2][0]*z,
            mInv.m[0][1]*x + mInv.m[1][1]*y + mInv.m[2][1]*z,
            mInv.m[0][2]*x + mInv.m[1][2]*y + mInv.m[2][2]*z);
    }

    inline void operator()(const Direction<S1> &n, Direction<S2> *nt) const {
            float x = n.x, y = n.y, z = n.z;
            nt->x = mInv.m[0][0] * x + mInv.m[1][0] * y +
                mInv.m[2][0] * z;
            nt->y = mInv.m[0][1] * x + mInv.m[1][1] * y +
                mInv.m[2][1] * z;
            nt->z = mInv.m[0][2] * x + mInv.m[1][2] * y +
                mInv.m[2][2] * z;
    }

    inline Direction<S2> operator()(const DiffNormal<S1> &n) const {
        float x = n.x, y = n.y, z = n.z;
        return Direction<S2>(mInv.m[0][0]*x + mInv.m[1][0]*y + mInv.m[2][0]*z,
            mInv.m[0][1]*x + mInv.m[1][1]*y + mInv.m[2][1]*z,
            mInv.m[0][2]*x + mInv.m[1][2]*y + mInv.m[2][2]*z);
    }

    inline void operator()(const DiffNormal<S1> &n, DiffNormal<S2> *nt) const {
        float x = n.x, y = n.y, z = n.z;
        nt->x = mInv.m[0][0] * x + mInv.m[1][0] * y +
            mInv.m[2][0] * z;
        nt->y = mInv.m[0][1] * x + mInv.m[1][1] * y +
            mInv.m[2][1] * z;
        nt->z = mInv.m[0][2] * x + mInv.m[1][2] * y +
            mInv.m[2][2] * z;
    }

    inline Ray<S2> operator()(const Ray<S1> &r) const {
        Ray<S2> ret;
        (*this)(r.o, &ret.o);
        (*this)(r.d, &ret.d);
        ret.maxt = r.maxt;
        ret.mint = r.mint;
        ret.time = r.time;
        ret.depth = r.depth;
        return ret;
    }


    inline void operator()(const Ray<S1> &r, Ray<S2> *rt) const {
        (*this)(r.o, &rt->o);
        (*this)(r.d, &rt->d);
        //if (rt != reinterpret_cast<const Ray<S2>*>(&r)) {
        rt->mint = r.mint;
        rt->maxt = r.maxt;
        rt->time = r.time;
        rt->depth = r.depth;
        //}
    }

    inline void operator()(const RayDifferential<S1> &r, RayDifferential<S2> *rt) const {
        (*this)(Ray<S1>(r), rt);
        rt->hasDifferentials = r.hasDifferentials;
        (*this)(r.rxOrigin, &rt->rxOrigin);
        (*this)(r.ryOrigin, &rt->ryOrigin);
        (*this)(r.rxDirection, &rt->rxDirection);
        (*this)(r.ryDirection, &rt->ryDirection);
    }



    inline RayDifferential<S2> operator()(const RayDifferential<S1> &r) const {
        RayDifferential<S2> ret;
        (*this)(Ray<S1>(r), &ret);
        ret.hasDifferentials = r.hasDifferentials;
        (*this)(r.rxOrigin, &ret.rxOrigin);
        (*this)(r.ryOrigin, &ret.ryOrigin);
        (*this)(r.rxDirection, &ret.rxDirection);
        (*this)(r.ryDirection, &ret.ryDirection);
        return ret;
    }

    BBox<S2> operator()(const BBox<S1> &b) const {
        const Transform<S1, S2> &M = *this;
        BBox<S2> ret(    M(Point<S1>(b.pMin.x, b.pMin.y, b.pMin.z)));
        ret = Union(ret, M(Point<S1>(b.pMax.x, b.pMin.y, b.pMin.z)));
        ret = Union(ret, M(Point<S1>(b.pMin.x, b.pMax.y, b.pMin.z)));
        ret = Union(ret, M(Point<S1>(b.pMin.x, b.pMin.y, b.pMax.z)));
        ret = Union(ret, M(Point<S1>(b.pMin.x, b.pMax.y, b.pMax.z)));
        ret = Union(ret, M(Point<S1>(b.pMax.x, b.pMax.y, b.pMin.z)));
        ret = Union(ret, M(Point<S1>(b.pMax.x, b.pMin.y, b.pMax.z)));
        ret = Union(ret, M(Point<S1>(b.pMax.x, b.pMax.y, b.pMax.z)));
        return ret;
    }

    template<typename S1_, typename S2_, typename S3_>
    friend Transform<S1_, S3_> operator*(const Transform<S2_, S3_> &t1, const Transform<S1_, S2_> &t2);


    void Print(FILE *f) const {
        m.Print(f);
    }

    bool SwapsHandedness() const {
        float det = ((m.m[0][0] *
            (m.m[1][1] * m.m[2][2] -
            m.m[1][2] * m.m[2][1])) -
            (m.m[0][1] *
            (m.m[1][0] * m.m[2][2] -
            m.m[1][2] * m.m[2][0])) +
            (m.m[0][2] *
            (m.m[1][0] * m.m[2][1] -
            m.m[1][1] * m.m[2][0])));
        return det < 0.f;
    }

private:
    // Transform Private Data
    Matrix4x4 m, mInv;
    template<typename S1_, typename S2_> friend class AnimatedTransform;
    template<typename S_> friend struct Quaternion;
};

template<typename S1, typename S2>    
inline Transform<S2, S1> Inverse(const Transform<S1, S2> &t) {
    return Transform<S2, S1>(t.mInv, t.m);
}

template<typename S1, typename S2, typename S3>
inline Transform<S1, S3> operator*(const Transform<S2, S3> &t1, const Transform<S1, S2> &t2)
{
    Matrix4x4 m1 = Matrix4x4::Mul(t1.m, t2.m);
    Matrix4x4 m2 = Matrix4x4::Mul(t2.mInv, t1.mInv);
    return Transform<S1, S3>(m1, m2);
}


template<typename S1, typename S2>
Transform<S1, S2> Translate(const Vector<S2> &delta) {
    Matrix4x4 m(1, 0, 0, __asfloat(delta.x),
        0, 1, 0, __asfloat(delta.y),
        0, 0, 1, __asfloat(delta.z),
        0, 0, 0,       1);
    Matrix4x4 minv(1, 0, 0, __asfloat(-delta.x),
        0, 1, 0, __asfloat(-delta.y),
        0, 0, 1, __asfloat(-delta.z),
        0, 0, 0,        1);
    return Transform<S1, S2>(m, minv);
}

template<typename S1, typename S2>
Transform<S1, S2> Scale(const float x, const float y, const float z) {
    Matrix4x4 m(x, 0, 0, 0,
                0, y, 0, 0,
                0, 0, z, 0,
                0, 0, 0, 1);
    Matrix4x4 minv(1.f/x,     0,     0,     0,
                   0,     1.f/y,     0,     0,
                   0,         0,     1.f/z, 0,
                   0,         0,     0,     1);
    return Transform<S1, S2>(m, minv);
}

template<typename S1, typename S2>
Transform<S1, S2> RotateX(float angle) {
    float sin_t = sin(Radians(angle));
    float cos_t = cos(Radians(angle));
    Matrix4x4 m(1,     0,      0, 0,
                0, cos_t, -sin_t, 0,
                0, sin_t,  cos_t, 0,
                0,     0,      0, 1);
    return Transform<S1, S2>(m, Transpose(m));
}

template<typename S1, typename S2>
Transform<S1, S2> RotateY(float angle) {
    float sin_t = sin(Radians(angle));
    float cos_t = cos(Radians(angle));
    Matrix4x4 m( cos_t,   0,  sin_t, 0,
                 0,   1,      0, 0,
                -sin_t,   0,  cos_t, 0,
                 0,   0,   0,    1);
    return Transform<S1, S2>(m, Transpose(m));
}


template<typename S1, typename S2>
Transform<S1, S2> RotateZ(float angle) {
    float sin_t = sin(Radians(angle));
    float cos_t = cos(Radians(angle));
    Matrix4x4 m(cos_t, -sin_t, 0, 0,
                sin_t,  cos_t, 0, 0,
                0,      0, 1, 0,
                0,      0, 0, 1);
    return Transform<S1, S2>(m, Transpose(m));
}

template<typename S1, typename S2>
Transform<S1, S2> Rotate(float angle, const Direction<S2> &axis) {
    Direction<S2> a = axis;
    float s = sin(_Radians(angle));
    float c = cos(_Radians(angle));
    float m[4][4];

    m[0][0] = a.x * a.x + (1.f - a.x * a.x) * c;
    m[0][1] = a.x * a.y * (1.f - c) - a.z * s;
    m[0][2] = a.x * a.z * (1.f - c) + a.y * s;
    m[0][3] = 0;

    m[1][0] = a.x * a.y * (1.f - c) + a.z * s;
    m[1][1] = a.y * a.y + (1.f - a.y * a.y) * c;
    m[1][2] = a.y * a.z * (1.f - c) - a.x * s;
    m[1][3] = 0;

    m[2][0] = a.x * a.z * (1.f - c) - a.y * s;
    m[2][1] = a.y * a.z * (1.f - c) + a.x * s;
    m[2][2] = a.z * a.z + (1.f - a.z * a.z) * c;
    m[2][3] = 0;

    m[3][0] = 0;
    m[3][1] = 0;
    m[3][2] = 0;
    m[3][3] = 1;

    Matrix4x4 mat(m);
    return Transform<S1, S2>(mat, Transpose(mat));
}

template<typename S1, typename S2>
Transform<S1, S2> LookAt(const Point<S2> &pos, const Point<S2> &look, const Direction<S2> &up) {
    float m[4][4];
    // Initialize fourth column of viewing matrix
    m[0][3] = pos.x.v;
    m[1][3] = pos.y.v;
    m[2][3] = pos.z.v;
    m[3][3] = 1;

    // Initialize first three columns of viewing matrix
    Direction<S2> dir = Normalize<S2>(look - pos);
    Direction<S2> left = Orthogonal<S2>(up,  dir);
    Direction<S2> newUp = Orthogonal<S2>(dir, left);
    m[0][0] = left.x;
    m[1][0] = left.y;
    m[2][0] = left.z;
    m[3][0] = 0.;
    m[0][1] = newUp.x;
    m[1][1] = newUp.y;
    m[2][1] = newUp.z;
    m[3][1] = 0.;
    m[0][2] = dir.x;
    m[1][2] = dir.y;
    m[2][2] = dir.z;
    m[3][2] = 0.;
    Matrix4x4 camToWorld(m);
    return Transform<S1, S2>(Inverse(camToWorld), camToWorld);
}


bool SolveLinearSystem2x2(const mfloat<length_d> A[2][2],
                          const mfloat<length_d> B[2], float *x0, float *x1);
bool SolveLinearSystem2x2(const float A[2][2],
                          const float B[2], float *x0, float *x1);

inline Transform<camera_s, screen_s> Orthographic(const mfloat<length_d> &znear, const mfloat<length_d> &zfar) {
    return Scale<screen_s, screen_s>(1.f, 1.f, 1.f / __asfloat(zfar-znear)) *
        Translate<camera_s, screen_s>(Vector<screen_s>(mfloat<length_d>(0.0f), mfloat<length_d>(0.0f), mfloat<length_d>(znear)));
}

inline Transform<camera_s, screen_s> Perspective(float fov, const mfloat<length_d> &n, const mfloat<length_d> &f) {
    // Perform projective divide
    Matrix4x4 persp = Matrix4x4(1, 0,           0,              0,
        0, 1,           0,              0,
        0, 0, f.v / (f - n).v, (-f*n).v / (f - n).v,
        0, 0,           1,              0);

    // Scale to canonical viewing volume
    float invTanAng = 1.f / tanf(_Radians(fov) / 2.f);
    return Scale<screen_s, screen_s>(invTanAng, invTanAng, 1) * Transform<camera_s, screen_s>(persp);
}


#include "quaternion.h"
// AnimatedTransform Declarations
template<typename S1, typename S2>
class AnimatedTransform {
public:
    // AnimatedTransform Public Methods
    AnimatedTransform(const Transform<S1, S2> *transform1, const mfloat<time_d> &time1,
                      const Transform<S1, S2> *transform2, const mfloat<time_d> &time2)
        : startTime(time1), endTime(time2),
          startTransform(transform1), endTransform(transform2),
          actuallyAnimated(*startTransform != *endTransform) {
        Decompose(startTransform->m, &T[0], &R[0], &S[0]);
        Decompose(endTransform->m, &T[1], &R[1], &S[1]);
    }

    void Interpolate(const mfloat<time_d> &time, Transform<S1, S2> *t) const {
        // Handle boundary conditions for matrix interpolation
        if (!actuallyAnimated || time <= startTime) {
            *t = *startTransform;
            return;
        }
        if (time >= endTime) {
            *t = *endTransform;
            return;
        }
        float dt = (time - startTime) / (endTime - startTime);
        // Interpolate translation at _dt_
        Vector<S2> trans = (1.f - dt) * T[0] + dt * T[1];

        // Interpolate rotation at _dt_
        Quaternion<S1> rotate = Slerp(dt, R[0], R[1]);

        // Interpolate scale at _dt_
        Matrix4x4 scale;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                scale.m[i][j] = Lerp(dt, S[0].m[i][j], S[1].m[i][j]);

        // Compute interpolated matrix as product of interpolated components
        *t = Translate<S1, S2>(trans) * rotate.ToTransform() * Transform<S1, S1>(scale);
    }

    void InterpolateInverse(const mfloat<time_d> &time, Transform<S1, S2> *t) const;

    void operator()(const Ray<S1> &r, Ray<S2> *tr) const
    {
        if (!actuallyAnimated || r.time <= startTime)
            (*startTransform)(r, tr);
        else if (r.time >= endTime)
            (*endTransform)(r, tr);
        else {
            Transform<S1, S2> t;
            Interpolate(r.time, &t);
            t(r, tr);
        }
        tr->time = r.time;
    }

    void operator()(const RayDifferential<S1> &r,
        RayDifferential<S2> *tr) const {
            if (!actuallyAnimated || r.time <= startTime)
                (*startTransform)(r, tr);
            else if (r.time >= endTime)
                (*endTransform)(r, tr);
            else {
                Transform<S1, S2> t;
                Interpolate(r.time, &t);
                t(r, tr);
            }
            tr->time = r.time;
    }


    Point<S2> operator()(mfloat<time_d> time, const Point<S1> &p) const {
        if (!actuallyAnimated || time <= startTime)
            return (*startTransform)(p);
        else if (time >= endTime)
            return (*endTransform)(p);
        Transform<S1, S2> t;
        Interpolate(time, &t);
        return t(p);
    }


    Vector<S2> operator()(mfloat<time_d> time, const Vector<S1> &v) const {
        if (!actuallyAnimated || time <= startTime)
            return (*startTransform)(v);
        else if (time >= endTime)
            return (*endTransform)(v);
        Transform<S1, S2> t;
        Interpolate(time, &t);
        return t(v);
    }

    Ray<S2> operator()(const Ray<S1> &r) const {
        Ray<S2> ret;
        (*this)(r, &ret);
        return ret;
    }

    BBox<S2> MotionBounds(const BBox<S1> &b) const 
    {
            if (!actuallyAnimated) return Inverse(*startTransform)(b);
            BBox<S2> ret;
            const int nSteps = 128;
            for (int i = 0; i < nSteps; ++i) {
                Transform<S1, S2> t;
                mfloat<time_d> time = Lerp(float(i)/float(nSteps-1), startTime, endTime);
                Interpolate(time, &t);
                ret = Union(ret, t(b));
            }
            return ret;
    }

    BBox<S1> InverseMotionBounds(const BBox<S2> &b) const 
    {
            if (!actuallyAnimated) return Inverse(*startTransform)(b);
            BBox<S1> ret;
            const int nSteps = 128;
            for (int i = 0; i < nSteps; ++i) {
                Transform<S1, S2> t;
                mfloat<time_d> time = Lerp(float(i)/float(nSteps-1), startTime, endTime);
                Interpolate(time, &t);
                Transform<S2, S1> it = Inverse(t);
                ret = Union(ret, it(b));
            }
            return ret;
    }


    bool HasScale() const { return startTransform->HasScale() || endTransform->HasScale(); }
private:
    // AnimatedTransform Method Definitions
    static void Decompose(const Matrix4x4 &m, Vector<S2> *T,
        Quaternion<S1> *Rquat, Matrix4x4 *S) {
            // Extract transformation _T_ from transformation matrix
            T->x = mfloat<length_d>(m.m[0][3]);
            T->y = mfloat<length_d>(m.m[1][3]);
            T->z = mfloat<length_d>(m.m[2][3]);

            // Compute new transformation matrix _M_ without translation
            Matrix4x4 M = m;
            for (int i = 0; i < 3; ++i)
                M.m[i][3] = M.m[3][i] = 0.f;
            M.m[3][3] = 1.f;

            // Extract rotation _R_ from transformation matrix
            float norm;
            int count = 0;
            Matrix4x4 R = M;
            do {
                // Compute next matrix _Rnext_ in series
                Matrix4x4 Rnext;
                Matrix4x4 Rit = Inverse(Transpose(R));
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                        Rnext.m[i][j] = 0.5f * (R.m[i][j] + Rit.m[i][j]);

                // Compute norm of difference between _R_ and _Rnext_
                norm = 0.f;
                for (int i = 0; i < 3; ++i) {
                    float n = fabsf(R.m[i][0] - Rnext.m[i][0]) +
                        fabsf(R.m[i][1] - Rnext.m[i][1]) +
                        fabsf(R.m[i][2] - Rnext.m[i][2]);
                    norm = max(norm, n);
                }
                R = Rnext;
            } while (++count < 100 && norm > .0001f);
            // XXX TODO FIXME deal with flip...
            *Rquat = Quaternion<S1>(R);

            // Compute scale _S_ using rotation and original matrix
            *S = Matrix4x4::Mul(Inverse(R), M);
    }

    // AnimatedTransform Private Data
    const mfloat<time_d> startTime, endTime;
    const Transform<S1, S2> *startTransform, *endTransform;
    const bool actuallyAnimated;
    Vector<S2> T[2];
    Quaternion<S1> R[2];
    Matrix4x4 S[2];
};

#endif // PBRT_CORE_MTRANSFORM_H
