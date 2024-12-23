
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

#ifndef PBRT_CORE_SPECTRUM_H
#define PBRT_CORE_SPECTRUM_H

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

// Spectrum Declarations
template <int nSamples> class CoefficientSpectrum {
public:
    // CoefficientSpectrum Public Methods
    CoefficientSpectrum(float v = 0.f) {
        for (int i = 0; i < nSamples; ++i)
            c[i] = v;
        Assert(!HasNaNs());
    }
#ifdef DEBUG
    CoefficientSpectrum(const CoefficientSpectrum &s) {
        Assert(!s.HasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] = s.c[i];
    }
    
    CoefficientSpectrum &operator=(const CoefficientSpectrum &s) {
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
    CoefficientSpectrum &operator+=(const CoefficientSpectrum &s2) {
        Assert(!s2.HasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] += s2.c[i];
        return *this;
    }
    CoefficientSpectrum operator+(const CoefficientSpectrum &s2) const {
        Assert(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] += s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator-(const CoefficientSpectrum &s2) const {
        Assert(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] -= s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator/(const CoefficientSpectrum &s2) const {
        Assert(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] /= s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator*(const CoefficientSpectrum &sp) const {
        Assert(!sp.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] *= sp.c[i];
        return ret;
    }
    CoefficientSpectrum &operator*=(const CoefficientSpectrum &sp) {
        Assert(!sp.HasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] *= sp.c[i];
        return *this;
    }
    CoefficientSpectrum operator*(float a) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] *= a;
        Assert(!ret.HasNaNs());
        return ret;
    }
    CoefficientSpectrum &operator*=(float a) {
        for (int i = 0; i < nSamples; ++i)
            c[i] *= a;
        Assert(!HasNaNs());
        return *this;
    }
    friend inline
    CoefficientSpectrum operator*(float a, const CoefficientSpectrum &s) {
        Assert(!isnan(a) && !s.HasNaNs());
        return s * a;
    }
    CoefficientSpectrum operator/(float a) const {
        Assert(!isnan(a));
        return *this * (1.f / a);
    }
    CoefficientSpectrum &operator/=(float a) {
        Assert(!isnan(a));
        float inv = 1.f / a;
        for (int i = 0; i < nSamples; ++i)
            c[i] *= inv;
        return *this;
    }
    bool operator==(const CoefficientSpectrum &sp) const {
        for (int i = 0; i < nSamples; ++i)
            if (c[i] != sp.c[i]) return false;
        return true;
    }
    bool operator!=(const CoefficientSpectrum &sp) const {
        return !(*this == sp);
    }
    bool IsBlack() const {
        for (int i = 0; i < nSamples; ++i)
            if (c[i] != 0.) return false;
        return true;
    }
    friend CoefficientSpectrum Sqrt(const CoefficientSpectrum &s) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = sqrtf(s.c[i]);
        Assert(!ret.HasNaNs());
        return ret;
    }
    template <int n> friend inline CoefficientSpectrum<n> Pow(const CoefficientSpectrum<n> &s, float e);
    CoefficientSpectrum operator-() const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = -c[i];
        return ret;
    }
    friend CoefficientSpectrum Exp(const CoefficientSpectrum &s) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = expf(s.c[i]);
        Assert(!ret.HasNaNs());
        return ret;
    }
    CoefficientSpectrum Clamp(float low = 0, float high = INFINITY) const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = ::Clamp(c[i], low, high);
        Assert(!ret.HasNaNs());
        return ret;
    }
    bool HasNaNs() const {
        for (int i = 0; i < nSamples; ++i)
            if (isnan(c[i])) return true;
        return false;
    }
    bool Write(FILE *f) const {
        for (int i = 0; i < nSamples; ++i)
            if (fprintf(f, "%f ", c[i]) < 0) return false;
        return true;
    }
    bool Read(FILE *f) {
        for (int i = 0; i < nSamples; ++i)
            if (fscanf(f, "%f ", &c[i]) != 1) return false;
        return true;
    }
protected:
    // CoefficientSpectrum Protected Data
    float c[nSamples];
};

class RGBSpectrum : public CoefficientSpectrum<3> {
    using CoefficientSpectrum<3>::c;
public:
    // RGBSpectrum Public Methods
    RGBSpectrum(float v = 0.f) : CoefficientSpectrum<3>(v) { }
    RGBSpectrum(const CoefficientSpectrum<3> &v)
        : CoefficientSpectrum<3>(v) { }
    RGBSpectrum(const RGBSpectrum &s, SpectrumType type = SPECTRUM_REFLECTANCE) {
        *this = s;
    }
    static RGBSpectrum FromRGB(const float rgb[3],
            SpectrumType type = SPECTRUM_REFLECTANCE) {
        RGBSpectrum s;
        s.c[0] = rgb[0];
        s.c[1] = rgb[1];
        s.c[2] = rgb[2];
        Assert(!s.HasNaNs());
        return s;
    }
    void ToRGB(float *rgb) const {
        rgb[0] = c[0];
        rgb[1] = c[1];
        rgb[2] = c[2];
    }
    const RGBSpectrum &ToRGBSpectrum() const {
        return *this;
    }
    void ToXYZ(float xyz[3]) const {
        RGBToXYZ(c, xyz);
    }
    static RGBSpectrum FromXYZ(const float xyz[3],
            SpectrumType type = SPECTRUM_REFLECTANCE) {
        RGBSpectrum r;
        XYZToRGB(xyz, r.c);
        return r;
    }
    float y() const {
        const float YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
        return YWeight[0] * c[0] + YWeight[1] * c[1] + YWeight[2] * c[2];
    }
    static RGBSpectrum FromSampled(const float *lambda, const float *v,
                                   int n) {
        // Sort samples if unordered, use sorted for returned spectrum
        if (!SpectrumSamplesSorted(lambda, v, n)) {
            vector<float> slambda(&lambda[0], &lambda[n]);
            vector<float> sv(&v[0], &v[n]);
            SortSpectrumSamples(&slambda[0], &sv[0], n);
            return FromSampled(&slambda[0], &sv[0], n);
        }
        float xyz[3] = { 0, 0, 0 };
        float yint = 0.f;
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



// Spectrum Inline Functions
template <int nSamples> inline CoefficientSpectrum<nSamples>
Pow(const CoefficientSpectrum<nSamples> &s, float e) {
    CoefficientSpectrum<nSamples> ret;
    for (int i = 0; i < nSamples; ++i)
        ret.c[i] = powf(s.c[i], e);
    Assert(!ret.HasNaNs());
    return ret;
}


inline Spectrum Lerp(float t, const Spectrum &s1, const Spectrum &s2) {
    return (1.f - t) * s1 + t * s2;
}



#endif // PBRT_CORE_SPECTRUM_H
