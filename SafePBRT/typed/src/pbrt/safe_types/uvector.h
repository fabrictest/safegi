#ifndef _SAFE_INTERNAL_VECTOR
#define _SAFE_INTERNAL_VECTOR

#include<m_float.h>

template<typename S> class DiffNormal;
template<typename S> class Direction;
template<typename S> class Normal;
template<typename S> class Point;
template<typename S> class Vector;

template<typename D>
class _Vector{ //internal used vector has only dimension
public:
    mfloat<D> x, y, z;
    // Vector Public Methods
    _Vector() : x(), y(), z() {}
    _Vector(const mfloat<D> &xx, const mfloat<D> &yy, const mfloat<D> &zz)
    {
        x = xx;
        y = yy;
        z = zz;
        Assert(!HasNaNs());
    }
    inline bool HasNaNs() const { return isnan(x.v) || isnan(y.v) || isnan(z.v); }
    inline _Vector operator+(const _Vector &v) const { return _Vector(x + v.x, y + v.y, z + v.z); }
    inline _Vector operator-(const _Vector &v) const { return _Vector(x - v.x, y - v.y, z - v.z); }
    inline _Vector& operator+=(const _Vector &v) { x += v.x; y += v.y; z += v.z;  return *this; }
    inline _Vector& operator-=(const _Vector &v) { x -= v.x; y -= v.y; z -= v.z;  return *this; }
    inline bool operator==(const _Vector &v) const { return x == v.x && y == v.y && z == v.z; }
    inline _Vector operator*(float f) const { return _Vector(f*x, f*y, f*z); }
    inline _Vector &operator*=(float f) { x *= f; y *= f; z *= f; return *this; }
    inline _Vector operator/(float f) const {  return _Vector(x / f, y / f, z / f); }
    inline _Vector &operator/=(float f) {  x /= f; y /= f; z /= f; return *this;  }  
    inline _Vector operator-() const { return _Vector(-x, -y, -z); }
    template<typename S>
    inline _Vector operator*(const Direction<S> &d) const { return _Vector(x*d.x, y*d.y, z*d.z); }
    template<typename S>
    inline _Vector operator/(const Direction<S> &d) const { return _Vector(x/d.x, y/d.y, z/d.z); }
    template<typename S>
    inline _Vector<typename dimension_product<D, length_d>::type> operator*(const Vector<S> &d) const { 
        return _Vector(x-d.x, y-d.y, z-d.z); 
    }
    template<typename S>
    inline _Vector<typename dimension_division<D, length_d>::type> operator/(const Vector<S> &d) const { 
        return _Vector(x-d.x, y-d.y, z-d.z); 
    }
};

template<typename D1>
inline _Vector<D1> operator*(float f, const _Vector<D1> &v)
{
    return v * f;
}
template<typename D1, typename S>
inline _Vector<D1> operator*(Direction<S> &d, const _Vector<D1> &v)
{
    return v * d;
}

template<typename D1, typename D2>
inline _Vector<typename dimension_product<D1, D2>::type> operator*(const _Vector<D1> &v1, const _Vector<D2> &v2)
{
    return _Vector<typename dimension_product<D1, D2>::type>(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}
template<typename D1, typename D2>
inline _Vector<typename dimension_division<D2, D1>::type> operator/(const mfloat<D2> &f, const _Vector<D1> &v)
{
    return _Vector<typename dimension_division<D2, D1>::type>(f / v.x, f / v.y, f / v.z);
}
template<typename D1, typename D2>
inline _Vector<typename dimension_product<D1, D2>::type> operator*(const mfloat<D2> &f, const _Vector<D1> &v)
{
    return v * f.v;
}
template<typename D1, typename D2>
inline _Vector<typename dimension_division<D1, D2>::type> operator/(const _Vector<D1> &v, const mfloat<D2> &f)
{
    return _Vector<typename dimension_division<D1, D2>::type>(v.x / f, v.y / f, v.z / f);
}
template<typename D1, typename D2>
inline _Vector<typename dimension_division<D1, D2>::type> operator/(const _Vector<D1> &v1, const _Vector<D2> &v2)
{
    return _Vector<typename dimension_division<D1, D2>::type>(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
}
template<typename D1>
inline _Vector<typename dimension_sqrt<D1>::type> Sqrt(const _Vector<D1> &v1)
{
    return _Vector<typename dimension_sqrt<D1>::type>(sqrt(v1.x), sqrt(v1.y), sqrt(v1.z));
}
template<typename D1, typename D2>
inline mfloat<typename dimension_product<D1, D2>::type> Dot(const _Vector<D1> &v1, const _Vector<D2> &v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
template<typename S, typename D2>
inline mfloat<typename dimension_product<length_d, D2>::type> Dot(const Vector<S> &v1, const _Vector<D2> &v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
template<typename S>
inline _Vector<area_d> Cross(const Vector<S> &v1, const Vector<S> &v2)
{
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    return _Vector<area_d>((v1.y * v2.z) - (v1.z * v2.y),
        (v1.z * v2.x) - (v1.x * v2.z),
        (v1.x * v2.y) - (v1.y * v2.x));
}
template<typename S, typename D>
inline _Vector<typename dimension_product<length_d, D>::type> operator*(const Vector<S> &v1, const mfloat<D> &f)
{
    return _Vector<typename dimension_product<length_d, D>::type>
        (v1.x * f, v1.y * f, v1.z * f);
}
template<typename S, typename D>
inline _Vector<typename dimension_product<length_d, D>::type> operator*(const mfloat<D> &f, const Vector<S> &v1)
{
    return v1 * f;
}
template<typename S>
inline _Vector<unit_d> operator+(const Direction<S> &v1, const _Vector<unit_d> &v2)
{
    return _Vector<unit_d>(v1.x - v2.x.v, v1.y - v2.y.v, v1.z - v2.z.v);
}
template<typename S>
inline _Vector<unit_d> operator-(const Direction<S> &v1, const _Vector<unit_d> &v2)
{
    return _Vector<unit_d>(v1.x - v2.x.v, v1.y - v2.y.v, v1.z - v2.z.v);
}
template<typename S>
inline _Vector<unit_d> operator-(const _Vector<unit_d> &v1, const Direction<S> &v2)
{
    return _Vector<unit_d>(v1.x.v - v2.x, v1.y.v - v2.y, v1.z.v - v2.z);
}
template<typename S>
inline _Vector<unit_d> operator-(const Direction<S> &v1, const Direction<S> &v2)
{
    return _Vector<unit_d>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}
template<typename S>
inline _Vector<unit_d> operator*(const mfloat<unit_d> &f, const Direction<S> &d) 
{
    return d * f;
}
template<typename S>
inline _Vector<unit_d> operator*(const Direction<S> &d, const mfloat<unit_d> &f) 
{
    return _Vector<unit_d>(f.v*d.x, f.v*d.y, f.v*d.z);
}
template<typename S, typename D>
mfloat<D> Dot(const Direction<S> &v1, const _Vector<D> &v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
template<typename S, typename D>
mfloat<D> Dot(const _Vector<D> &v2, const Direction<S> &v1) {
    return Dot(v1, v2);
}
template<typename S, typename D>
inline Direction<S> Normalize(const _Vector<D> &v) 
{
    mfloat<D> l = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    return Direction<S>(v.x/l, v.y/l, v.z/l); 
}
template<typename S>
inline Vector<S> AsVector(const _Vector<length_d> &v) 
{
    return Vector<S>(v.x, v.y, v.z); 
}
template<typename S>
inline Normal<S> AsNormalDiff(const _Vector<unit_d> &v)
{
    return Normal<S>(v.x.v, v.y.v, v.z.v); 
}
#endif