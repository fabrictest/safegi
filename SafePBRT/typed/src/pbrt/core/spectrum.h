
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

#ifndef PBRT_CORE_MSPECTRUM_H
#define PBRT_CORE_MSPECTRUM_H

// core/spectrum.h*
#include "pbrt.h"
#include "parallel.h"

// Spectrum Utility Declarations
static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectralSamples = 30;
extern bool SpectrumSamplesSorted(const float *lambda, const float *vals, int n);
extern void SortSpectrumSamples(float *lambda, float *vals, int n);
extern float AverageSpectrumSamples(const float *lambda, const float *vals,
    int n, float lambdaStart, float lambdaEnd);
inline void XYZToRGB(const float xyz[3], float rgb[3]) {
    rgb[0] =  3.240479f*xyz[0] - 1.537150f*xyz[1] - 0.498535f*xyz[2];
    rgb[1] = -0.969256f*xyz[0] + 1.875991f*xyz[1] + 0.041556f*xyz[2];
    rgb[2] =  0.055648f*xyz[0] - 0.204043f*xyz[1] + 1.057311f*xyz[2];
}


inline void RGBToXYZ(const float rgb[3], float xyz[3]) {
    xyz[0] = 0.412453f*rgb[0] + 0.357580f*rgb[1] + 0.180423f*rgb[2];
    xyz[1] = 0.212671f*rgb[0] + 0.715160f*rgb[1] + 0.072169f*rgb[2];
    xyz[2] = 0.019334f*rgb[0] + 0.119193f*rgb[1] + 0.950227f*rgb[2];
}


enum SpectrumType { SPECTRUM_REFLECTANCE, SPECTRUM_ILLUMINANT };
extern void Blackbody(const float *wl, int n, float temp, float *vals);
extern float InterpolateSpectrumSamples(const float *lambda, const float *vals,
                                        int n, float l);

// Spectral Data Declarations
static const int nCIESamples = 471;
extern const float CIE_X[nCIESamples];
extern const float CIE_Y[nCIESamples];
extern const float CIE_Z[nCIESamples];
extern const float CIE_lambda[nCIESamples];
static const int nRGB2SpectSamples = 32;
extern const float RGB2SpectLambda[nRGB2SpectSamples];
extern const float RGBRefl2SpectWhite[nRGB2SpectSamples];
extern const float RGBRefl2SpectCyan[nRGB2SpectSamples];
extern const float RGBRefl2SpectMagenta[nRGB2SpectSamples];
extern const float RGBRefl2SpectYellow[nRGB2SpectSamples];
extern const float RGBRefl2SpectRed[nRGB2SpectSamples];
extern const float RGBRefl2SpectGreen[nRGB2SpectSamples];
extern const float RGBRefl2SpectBlue[nRGB2SpectSamples];
extern const float RGBIllum2SpectWhite[nRGB2SpectSamples];
extern const float RGBIllum2SpectCyan[nRGB2SpectSamples];
extern const float RGBIllum2SpectMagenta[nRGB2SpectSamples];
extern const float RGBIllum2SpectYellow[nRGB2SpectSamples];
extern const float RGBIllum2SpectRed[nRGB2SpectSamples];
extern const float RGBIllum2SpectGreen[nRGB2SpectSamples];
extern const float RGBIllum2SpectBlue[nRGB2SpectSamples];

// Spectrum Declarations
template <typename D, int nSamples>
class CoefficientSpectrum {
public:
    CoefficientSpectrum() {
        for (int i = 0; i < nSamples; ++i)
            c[i] = mfloat<D>(0.0f);
    }
    // CoefficientSpectrum Public Methods
    explicit CoefficientSpectrum(mfloat<D> v) {
        for (int i = 0; i < nSamples; ++i)
            c[i] = v;
        Assert(!HasNaNs());
    }
#ifdef DEBUG
    CoefficientSpectrum(const CoefficientSpectrum<D> &s) {
        Assert(!s.HasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] = s.c[i];
    }
    
    CoefficientSpectrum<D> &operator=(const CoefficientSpectrum<D> &s) {
        Assert(!s.HasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] = s.c[i];
        return *this;
    }
#endif // DEBUG
    void Print(FILE *f) const {
        fprintf(f, "[ ");
        for (int i = 0; i < nSamples; ++i) {
            fprintf(f, "%f", c[i]);
            if (i != nSamples-1) fprintf(f, ", ");
        }
        fprintf(f, "]");
    }
    inline CoefficientSpectrum &operator+=(const CoefficientSpectrum &s2) {
        Assert(!s2.HasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i].v += s2.c[i].v;
        return *this;
    }
    inline CoefficientSpectrum operator+(const CoefficientSpectrum &s2) const {
        Assert(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i].v += s2.c[i].v;
        return ret;
    }
    inline CoefficientSpectrum operator-(const CoefficientSpectrum &s2) const {
        Assert(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i].v -= s2.c[i].v;
        return ret;
    }

    template<typename D1, typename D2, int nSamples_>
    friend CoefficientSpectrum<typename dimension_division<D1,D2>::type, nSamples_> operator/(const CoefficientSpectrum<D1, nSamples_> &s1, const CoefficientSpectrum<D2, nSamples_> &s2);
    template<typename D1, typename D2, int nSamples_>
    friend CoefficientSpectrum<typename dimension_product<D1,D2>::type, nSamples_> operator*(const CoefficientSpectrum<D1, nSamples_> &s1, const CoefficientSpectrum<D2, nSamples_> &s2);
    template<typename D1, typename D2, int nSamples_>
    friend CoefficientSpectrum<typename dimension_product<D1,D2>::type, nSamples_> operator*(const CoefficientSpectrum<D1, nSamples_> &s1, const mfloat<D2> &a);


    inline CoefficientSpectrum operator*(const float a) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i].v *= a;
        Assert(!ret.HasNaNs());
        return ret;
    }

    inline CoefficientSpectrum &operator*=(const float a) {
        for (int i = 0; i < nSamples; ++i)
            c[i].v *= a;
        Assert(!HasNaNs());
        return *this;
    }

    template<typename D_, int n> friend class CoefficientSpectrum;

    inline CoefficientSpectrum &operator*=(const CoefficientSpectrum<rho_d, nSamples> &a) {
        for (int i = 0; i < nSamples; ++i)
            c[i].v *= a.c[i].v;
        Assert(!HasNaNs());
        return *this;
    }

    friend 
    inline CoefficientSpectrum operator*(const float a, const CoefficientSpectrum &s) {
        Assert(!isnan(a) && !s.HasNaNs());
        return s * a;
    }
    template<typename D_> 
    friend inline CoefficientSpectrum<typename dimension_product<D,D_>::type, nSamples> operator*(const mfloat<D_> &a, const CoefficientSpectrum &s) {
            Assert(!isnan(a.v) && !s.HasNaNs());
            return s * a;
    }
    inline CoefficientSpectrum operator/(const float a) const {
        Assert(!isnan(a));
        return *this * (1.f / a);
    }
    template<typename D_> 
    inline CoefficientSpectrum<typename dimension_division<D,D_>::type, nSamples> operator/(const mfloat<D_> &a) const {
        Assert(!isnan(a.v));
        return *this * (1.f / a);
    }
    inline CoefficientSpectrum &operator/=(const float a) {
        Assert(!isnan(a));
        float inv = 1.f / a;
        for (int i = 0; i < nSamples; ++i)
            c[i].v *= inv;
        return *this;
    }

    inline CoefficientSpectrum &operator/=(const CoefficientSpectrum<rho_d, nSamples> &a) {
        for (int i = 0; i < nSamples; ++i)
            c[i] /= a.c[i];
        Assert(!HasNaNs());
        return *this;
    }

    inline bool operator==(const CoefficientSpectrum &sp) const {
        for (int i = 0; i < nSamples; ++i)
            if (c[i] != sp.c[i]) return false;
        return true;
    }
    inline bool operator!=(const CoefficientSpectrum &sp) const {
        return !(*this == sp);
    }
    inline bool IsBlack() const {
        for (int i = 0; i < nSamples; ++i)
            if (c[i] != mfloat<D>(0.0f)) return false;
        return true;
    }
    inline friend CoefficientSpectrum<typename dimension_sqrt<D>::type, nSamples> Sqrt(const CoefficientSpectrum &s) {
        CoefficientSpectrum<typename dimension_sqrt<D>::type, nSamples> ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i].v = sqrt(s.c[i].v);
        Assert(!ret.HasNaNs());
        return ret;
    }

    inline CoefficientSpectrum operator-() const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = -c[i];
        return ret;
    }

    inline CoefficientSpectrum Clamp(const mfloat<D> &low = mfloat<D>(0.0f), const mfloat<D> &high = mfloat<D>(INFINITY)) const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = ::Clamp(c[i], low, high);
        Assert(!ret.HasNaNs());
        return ret;
    }
    inline bool HasNaNs() const {
        for (int i = 0; i < nSamples; ++i)
            if (isnan(c[i].v)) return true;
        return false;
    }
    inline bool Write(FILE *f) const {
        for (int i = 0; i < nSamples; ++i)
            if (fprintf(f, "%f ", c[i]) < 0) return false;
        return true;
    }
    inline bool Read(FILE *f) {
        for (int i = 0; i < nSamples; ++i)
            if (fscanf(f, "%f ", &c[i]) != 1) return false;
        return true;
    }

    template<int nSamples_>
    friend inline CoefficientSpectrum<rho_d, nSamples_> Pow(const CoefficientSpectrum<rho_d, nSamples_> &s, float e);
    template<int nSamples_>
    friend inline CoefficientSpectrum<rho_d, nSamples_> Exp(const CoefficientSpectrum<rho_d, nSamples_> &s);
protected:
    // CoefficientSpectrum Protected Data
    mfloat<D> c[nSamples];
};

template<int nSamples>
inline CoefficientSpectrum<rho_d, nSamples> Pow(const CoefficientSpectrum<rho_d, nSamples> &s, float e)
{
    CoefficientSpectrum<rho_d, nSamples> ret;
    for (int i = 0; i < nSamples; ++i)
        ret.c[i] = pow(s.c[i].v, e);
    Assert(!ret.HasNaNs());
    return ret;
}

template<int nSamples>
inline CoefficientSpectrum<rho_d, nSamples> Exp(const CoefficientSpectrum<rho_d, nSamples> &s)
{
    CoefficientSpectrum<rho_d, nSamples> ret;
    for (int i = 0; i < nSamples; ++i)
        ret.c[i] = expf(s.c[i].v);
    Assert(!ret.HasNaNs());
    return ret;
}

template<typename D1, typename D2, int nSamples>
inline CoefficientSpectrum<typename dimension_division<D1,D2>::type, nSamples> operator/(const CoefficientSpectrum<D1, nSamples> &s1, const CoefficientSpectrum<D2, nSamples> &s2) 
{
    Assert(!s2.HasNaNs() || !s1.HasNaNs());
    CoefficientSpectrum<typename dimension_division<D1,D2>::type, nSamples> ret;
    for (int i = 0; i < nSamples; ++i)
        ret.c[i].v = s1.c[i].v / s2.c[i].v;
    return ret;
}
template<typename D1, typename D2, int nSamples>
inline CoefficientSpectrum<typename dimension_product<D1,D2>::type, nSamples> operator*(const CoefficientSpectrum<D1, nSamples> &s1, const CoefficientSpectrum<D2, nSamples> &s2) 
{
    Assert(!s2.HasNaNs() || !s1.HasNaNs());
    CoefficientSpectrum<typename dimension_product<D1,D2>::type, nSamples> ret;
    for (int i = 0; i < nSamples; ++i)
        ret.c[i].v = s1.c[i].v * s2.c[i].v;
    return ret;
}
template<typename D1, typename D2, int nSamples>
inline CoefficientSpectrum<typename dimension_product<D1,D2>::type, nSamples> operator*(const CoefficientSpectrum<D1, nSamples> &s1, const mfloat<D2> &a)
{
    CoefficientSpectrum<typename dimension_product<D1,D2>::type, nSamples> ret;
    for (int i = 0; i < nSamples; ++i)
        ret.c[i].v = s1.c[i].v * a.v;
    Assert(!ret.HasNaNs());
    return ret;
}

template<typename D>
class RGBSpectrum : public CoefficientSpectrum<D, 3> {
    using CoefficientSpectrum<D, 3>::c;
public:
    // RGBSpectrum Public Methods
    RGBSpectrum() : CoefficientSpectrum<D, 3>(mfloat<D>(0.0f)) { }
    explicit RGBSpectrum(mfloat<D> v) : CoefficientSpectrum<D, 3>(v) { }
    RGBSpectrum(const CoefficientSpectrum<D, 3> &v)
        : CoefficientSpectrum<D, 3>(v) { }
    RGBSpectrum(const RGBSpectrum &s, SpectrumType type = SPECTRUM_REFLECTANCE) {
        *this = s;
    }
    inline static RGBSpectrum FromRGB(const float rgb[3],
            SpectrumType type = SPECTRUM_REFLECTANCE) {
        RGBSpectrum s;
        s.c[0] = mfloat<D>(rgb[0]);
        s.c[1] = mfloat<D>(rgb[1]);
        s.c[2] = mfloat<D>(rgb[2]);
        Assert(!s.HasNaNs());
        return s;
    }
    inline void ToRGB(float *rgb) const {
        rgb[0] = c[0].v;
        rgb[1] = c[1].v;
        rgb[2] = c[2].v;
    }
    inline const RGBSpectrum &ToRGBSpectrum() const {
        return *this;
    }
    inline void ToXYZ(float xyz[3]) const {
        float rgb[3];
        ToRGB(rgb);
        RGBToXYZ(rgb, xyz);
    }
    inline static RGBSpectrum<D> FromXYZ(const float xyz[3],
            SpectrumType type = SPECTRUM_REFLECTANCE) {
        RGBSpectrum r;
        float rgb[3];
        XYZToRGB(xyz, rgb);
        r.FromRGB(rgb, type);
        return r;
    }
    inline mfloat<D> y() const {
        const float YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
        return YWeight[0] * c[0] + YWeight[1] * c[1] + YWeight[2] * c[2];
    }
    static RGBSpectrum<D> FromSampled(const float *lambda, const float *v,
                                   int n) {
        // Sort samples if unordered, use sorted for returned spectrum
        if (!SpectrumSamplesSorted(lambda, v, n)) {
            std::vector<float> slambda(&lambda[0], &lambda[n]);
            std::vector<float> sv(&v[0], &v[n]);
            SortSpectrumSamples(&slambda[0], &sv[0], n);
            return FromSampled(&slambda[0], &sv[0], n);
        }
        float xyz[3] = { 0, 0, 0 };
        float yint = float();
        for (int i = 0; i < nCIESamples; ++i) {
            yint += CIE_Y[i];
            float val = InterpolateSpectrumSamples(lambda, v, n,
                                                   CIE_lambda[i]);
            xyz[0] += val * CIE_X[i];
            xyz[1] += val * CIE_Y[i];
            xyz[2] += val * CIE_Z[i];
        }
        xyz[0] /= yint;
        xyz[1] /= yint;
        xyz[2] /= yint;
        return FromXYZ(xyz);
    }
};

template<typename D>
inline Spectrum<D> Lerp(float t, const Spectrum<D> &s1, const Spectrum<D> &s2) {
    return (float(1.0) - t) * s1 + t * s2;
}


#endif // PBRT_CORE_MSPECTRUM_H
