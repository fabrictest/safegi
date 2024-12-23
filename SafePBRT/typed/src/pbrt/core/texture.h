
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

#ifndef PBRT_CORE_MTEXTURE_H
#define PBRT_CORE_MTEXTURE_H

// core/texture.h*
#include "pbrt.h"
#include "spectrum.h"
#include "geometry.h"
#include "transform.h"
#include "memory.h"

template<typename T>
struct texval{}; 
// This template object is use for template deduction in 
// texture template, since we do not allow implicit convertion
// between mfloat<> and Spectrum<>

template<>
struct texval<float>
{
    static inline float black() { return 0.0f; }
    static inline float white() { return 1.0f; }
    static inline float convert_float(float a) { return a; }
};

template<typename D>
struct texval<mfloat<D> >
{
    static inline mfloat<D> black() { return mfloat<D>(0.0f); }
    static inline mfloat<D> white() { return mfloat<D>(1.0f); }
    static inline mfloat<D> convert_float(float a) { return mfloat<D>(a); }
};

template<typename D>
struct texval<RGBSpectrum<D> >
{
    static inline RGBSpectrum<D> black() { return RGBSpectrum<D>(mfloat<D>(0.0f)); }
    static inline RGBSpectrum<D> white() { return RGBSpectrum<D>(mfloat<D>(1.0f)); }
    static inline RGBSpectrum<D> convert_float(float a) { return RGBSpectrum<D>(mfloat<D>(a)); }
};


// Texture Declarations
class TextureMapping2D {
public:
    // TextureMapping2D Interface
    virtual ~TextureMapping2D() { }
    virtual void Map(const DifferentialGeometry &dg,
                     float *s, float *t, float *dsdx, float *dtdx,
                     float *dsdy, float *dtdy) const = 0;
};


class UVMapping2D : public TextureMapping2D {
public:
    // UVMapping2D Public Methods
    UVMapping2D(float su = 1, float sv = 1, float du = 0, float dv = 0);
    void Map(const DifferentialGeometry &dg, float *s, float *t,
        float *dsdx, float *dtdx, float *dsdy, float *dtdy) const;
private:
    float su, sv, du, dv;
};


class SphericalMapping2D : public TextureMapping2D {
public:
    // SphericalMapping2D Public Methods
    SphericalMapping2D(const Transform<world_s, texture_s> &toSph)
        : WorldToTexture(toSph) {
    }
    void Map(const DifferentialGeometry &dg, float *s, float *t,
        float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const;
private:
    void sphere(const Point<world_s> &P, float *s, float *t) const;
    Transform<world_s, texture_s> WorldToTexture;
};


class CylindricalMapping2D : public TextureMapping2D {
public:
    // CylindricalMapping2D Public Methods
    CylindricalMapping2D(const Transform<world_s, texture_s> &toCyl)
        : WorldToTexture(toCyl) {
    }
    void Map(const DifferentialGeometry &dg, float *s, float *t,
        float *dsdx, float *dtdx, float *dsdy, float *dtdy) const;
private:
    // CylindricalMapping2D Private Methods
    void cylinder(const Point<world_s> &p, float *s, float *t) const {
        Direction<texture_s> vec = 
            Normalize(WorldToTexture(p) - Point<texture_s>());
        *s = (halfCircleAngle + _atan2(vec.y, vec.x)) / (circleAngle);
        *t = vec.z;
    }
    Transform<world_s, texture_s> WorldToTexture;
};


class PlanarMapping2D : public TextureMapping2D {
public:
    // PlanarMapping2D Public Methods
    void Map(const DifferentialGeometry &dg, float *s, float *t,
             float *dsdx, float *dtdx, float *dsdy, float *dtdy) const;
    PlanarMapping2D(const Direction<world_s> &vv1, const Direction<world_s> &vv2,
                    float dds = 0, float ddt = 0)
        : vs(vv1), vt(vv2), ds(dds), dt(ddt) { }
private:
    const Direction<world_s> vs, vt;
    const float ds, dt;
};


class TextureMapping3D {
public:
    // TextureMapping3D Interface
    virtual ~TextureMapping3D() { }
    virtual Point<texture_s> Map(const DifferentialGeometry &dg,
                      Vector<texture_s> *dpdx, Vector<texture_s> *dpdy) const = 0;
};


class IdentityMapping3D : public TextureMapping3D {
public:
    IdentityMapping3D(const Transform<world_s, texture_s> &x)
        : WorldToTexture(x) { }
    Point<texture_s> Map(const DifferentialGeometry &dg, Vector<texture_s> *dpdx,
              Vector<texture_s> *dpdy) const;
private:
    Transform<world_s, texture_s> WorldToTexture;
};

template <typename T> class Texture : public ReferenceCounted {
public:
    // Texture Interface
    virtual T Evaluate(const DifferentialGeometry &) const = 0;
    virtual ~Texture() { }
};

float Lanczos(float, float tau=2);
#endif // PBRT_CORE_MTEXTURE_H
