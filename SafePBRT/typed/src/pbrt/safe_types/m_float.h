#ifndef _M_REAL_H_
#define _M_REAL_H_

#include "dimensions.h"

#ifdef _WIN32
#include <memory>
#include <random>
#else
#include <tr1/memory>
#include <tr1/random>
#endif
using std::tr1::minstd_rand;
using std::tr1::uniform_real;

const float eps = 1e-6f;
const float inf = 1e12f;
const float pi = 3.1415926535897932384626433832795f;

struct tuple2 {
    float x, y;

    tuple2() : x(0), y(0) { }
    tuple2(const float& x, const float& y) : x(x), y(y) { }

    inline tuple2& operator+= (const float &a) { this->x += a; this->y += a; return *this;}
    inline tuple2& operator-= (const float &a) { this->x -= a; this->y -= a; return *this;}
    inline tuple2& operator*= (const float &a) { this->x *= a; this->y *= a; return *this;}
    inline tuple2& operator/= (const float &a) { this->x /= a; this->y /= a; return *this;}

    inline tuple2& operator+= (const tuple2 &a) { this->x *= a.x; this->y *= a.y; return *this;}
    inline tuple2& operator-= (const tuple2 &a) { this->x -= a.x; this->y -= a.y; return *this;}
    inline tuple2& operator*= (const tuple2 &a) { this->x *= a.x; this->y *= a.y; return *this;}
    inline tuple2& operator/= (const tuple2 &a) { this->x /= a.x; this->y /= a.y; return *this;}
};

struct tuple3 {
    float x, y, z;

    tuple3() : x(0), y(0), z(0) { }
    tuple3(const float& x, const float& y, const float& z) : x(x), y(y), z(z) { }
    inline tuple3& operator+= (const float &a) { this->x += a; this->y += a; this->z += a; return *this;}
    inline tuple3& operator-= (const float &a) { this->x -= a; this->y -= a; this->z -= a; return *this;}
    inline tuple3& operator*= (const float &a) { this->x *= a; this->y *= a; this->z *= a; return *this;}
    inline tuple3& operator/= (const float &a) { this->x /= a; this->y /= a; this->z /= a; return *this;}

    inline tuple3& operator+= (const tuple3 &a) { this->x += a.x; this->y += a.y; this->z += a.z; return *this;}
    inline tuple3& operator-= (const tuple3 &a) { this->x -= a.x; this->y -= a.y; this->z -= a.z; return *this;}
    inline tuple3& operator*= (const tuple3 &a) { this->x *= a.x; this->y *= a.y; this->z *= a.z; return *this;}
    inline tuple3& operator/= (const tuple3 &a) { this->x /= a.x; this->y /= a.y; this->z /= a.z; return *this;}
};

inline bool iszero(const float& v) { return v == 0; }

inline float urandom() {
    static minstd_rand eng;
    static uniform_real<float> dist;
    return dist(eng);
}

template<typename D>
struct mfloat {
    float v;

    mfloat<D>() : v(0) { }
    inline mfloat& operator+= (const mfloat& a) { this->v += a.v; return *this;}
    inline mfloat& operator-= (const mfloat& a) { this->v -= a.v; return *this;}
    inline mfloat& operator*= (const float &a) { this->v *= a; return *this;}
    inline mfloat& operator/= (const float &a) { this->v /= a; return *this;}
    explicit mfloat<D>(const float& v) : v(v) { }
};

template<>
struct mfloat<unit_d> {
    float v;

    mfloat() : v(0) { }
    mfloat(const float& v) : v(v) { }
    inline mfloat& operator+= (const mfloat& a) { this->v += a.v; return *this;}
    inline mfloat& operator-= (const mfloat& a) { this->v -= a.v; return *this;}
    inline mfloat& operator*= (const float &a) { this->v *= a; return *this;}
    inline mfloat& operator/= (const float &a) { this->v /= a; return *this;}
    operator float() { return v; }
};

template<typename D>
inline mfloat<D> operator+ (const mfloat<D>& a, const mfloat<D>& b) {
    return mfloat<D>(a.v+b.v);
}

template<typename D>
inline mfloat<D> operator- (const mfloat<D>& a, const mfloat<D>& b) {
    return mfloat<D>(a.v-b.v);
}

template<typename D>
inline mfloat<D> operator- (const mfloat<D>& a) {
    return mfloat<D>(-a.v);
}

template<typename D>
inline mfloat<D> operator* (const mfloat<D>& a, const float& b) {
    return mfloat<D>(a.v*b);
}

template<typename D>
inline mfloat<D> operator* (const float& a, const mfloat<D>& b) {
    return mfloat<D>(a*b.v);
}

template<typename D>
inline mfloat<D> operator/ (const mfloat<D>& a, const float& b) {
    return mfloat<D>(a.v/b);
}

template<typename D>
inline bool operator == (const mfloat<D>& a, const mfloat<D>& b) {
    return a.v == b.v;
}

template<typename D>
inline bool operator != (const mfloat<D>& a, const mfloat<D>& b) {
    return a.v != b.v;
}

template<typename D>
inline bool iszero(const mfloat<D>& a) {
    return a.v == 0;
}

template<typename D>
inline bool operator> (const mfloat<D>& a, const mfloat<D>& b) {
    return a.v > b.v;
}

template<typename D>
inline bool operator< (const mfloat<D>& a, const mfloat<D>& b) {
    return a.v < b.v;
}

template<typename D>
inline bool operator>= (const mfloat<D>& a, const mfloat<D>& b) {
    return a.v >= b.v;
}

template<typename D>
inline bool operator<= (const mfloat<D>& a, const mfloat<D>& b) {
    return a.v <= b.v;
}

template<typename D1, typename D2>
inline mfloat<typename dimension_product<D1,D2>::type> operator* (const mfloat<D1>& a, const mfloat<D2>& b) {
    return mfloat<typename dimension_product<D1,D2>::type>(a.v*b.v);
}

template<typename D1, typename D2>
inline mfloat<typename dimension_division<D1,D2>::type> operator/ (const mfloat<D1>& a, const mfloat<D2>& b) {
    return mfloat<typename dimension_division<D1,D2>::type>(a.v/b.v);
}

template<typename D>
inline float operator/ (const mfloat<D>& a, const mfloat<D>& b) {
    return a.v/b.v;
}

template<typename D>
inline mfloat<typename dimension_division<unit_d,D>::type> operator/ (const float& a, const mfloat<D>& b) {
    return mfloat<typename dimension_division<unit_d,D>::type>(a/b.v);
}

template<typename D>
inline mfloat<typename dimension_sqrt<D>::type> sqrt(const mfloat<D>& a) {
    return mfloat<typename dimension_sqrt<D>::type>(sqrt(a.v));
}

template<typename D>
inline mfloat<D> abs(const mfloat<D> &a)
{
    return mfloat<D>(abs(a.v));
}

inline mfloat<unit_d> tan(const mfloat<angle_d> &a)
{
    return tan(a.v);
}

inline mfloat<unit_d> sin(const mfloat<angle_d> &a)
{
    return sin(a.v);
}

inline mfloat<unit_d> cos(const mfloat<angle_d> &a)
{
    return cos(a.v);
}

inline mfloat<angle_d> _acos(const mfloat<unit_d> &a)
{
    return mfloat<angle_d>(acos(a.v));
}

inline mfloat<angle_d> _atan2(const mfloat<unit_d> &x, const mfloat<unit_d> &y)
{
    return mfloat<angle_d>(atan2f(x.v, y.v));
}

inline mfloat<angle_d> _atan(const mfloat<unit_d> &x)
{
    return mfloat<angle_d>(atanf(x.v));
}


template<typename D>
class mtuple3 {
public:
	mfloat<D> x, y, z;

    mtuple3() : x(), y(), z() { }
    mtuple3(const mfloat<D>& x, const mfloat<D>& y, const mfloat<D>& z) : x(x), y(y), z(z) { }
};

template<typename D>
inline float __asfloat(mfloat<D>& a) { return a.v; }
template<typename D>
inline const float __asfloat(const mfloat<D>& a) { return a.v; }

const mfloat<length_d> eps_length = mfloat<length_d>(eps);
const mfloat<length_d> inf_length = mfloat<length_d>(inf);

const mfloat<projsolidangle_d> hemisphericalProjectedAngle = mfloat<projsolidangle_d>(pi);
const mfloat<solidangle_d> sphericalAngle = mfloat<solidangle_d>(4.0f*pi);
const mfloat<solidangle_d> hemisphericalAngle = mfloat<solidangle_d>(2.0f*pi);
const mfloat<angle_d> circleAngle = mfloat<angle_d>(2.0f*pi);
const mfloat<angle_d> halfCircleAngle = mfloat<angle_d>(pi);

#endif
