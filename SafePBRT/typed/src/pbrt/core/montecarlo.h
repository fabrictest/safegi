
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

#ifndef PBRT_CORE_MONTECARLO_H
#define PBRT_CORE_MONTECARLO_H

// core/montecarlo.h*
#include "pbrt.h"
#include "geometry.h"
#include "rng.h"

int LDPixelSampleFloatsNeeded(const Sample *sample, int nPixelSamples);

// Monte Carlo Utility Declarations
template<typename D>
struct Distribution1D {
    // Distribution1D Public Methods
    Distribution1D(const mfloat<D> *f, int n) {
        count = n;
        func = new mfloat<D>[n];
        memcpy(func, f, n*sizeof(float));
        cdf = new float[n+1];
        std::vector<mfloat<D> > tmp(n+1, mfloat<D>(0.0f));
        // Compute integral of step function at $x_i$
        tmp[0] = mfloat<D>(0.0f);
        for (int i = 1; i < n+1; ++i)
            tmp[i] = tmp[i-1] + f[i-1] / float(n);

        // Transform step function integral into CDF
        funcInt = tmp[n];
        for (int i = 1; i < n+1; ++i)
            cdf[i] = tmp[i] / funcInt;
    }

    ~Distribution1D() {
        delete[] func;
        delete[] cdf;
    }

    float SampleContinuous(float u, float *pdf) const {
        // Find surrounding CDF segments and _offset_
        float *ptr = std::lower_bound(cdf, cdf+count+1, u);
        int offset = max(0, int(ptr-cdf-1));

        // Compute offset along CDF segment
        float du = (u - cdf[offset]) / (cdf[offset+1] - cdf[offset]);

        // Compute PDF for sampled offset
        if (pdf) *pdf = func[offset] / funcInt;

        // Return $x \in [0,1)$ corresponding to sample
        return (offset + du) / count;
    }

    int SampleDiscrete(float u, float *pdf) const {
        // Find surrounding CDF segments and _offset_
        float *ptr = std::lower_bound(cdf, cdf+count+1, u);
        int offset = max(0, int(ptr-cdf-1));

        // Compute PDF for sampled offset
        if (pdf) *pdf = func[offset] / funcInt;
        return offset;
    }
private:
    template<typename D_>
    friend struct Distribution2D;
    // Distribution1D Private Data
    mfloat<D> *func;
    float *cdf;
    mfloat<D> funcInt;
    int count;
};


void RejectionSampleDisk(float u1, float u2, float *x, float *y);

template<typename S>
Direction<S> UniformSampleHemisphere(const float &u1, const float &u2) {
    float z = u1;
    float r = sqrt(max(0.0f, 1.0f - z*z));
    float phi = M_PI * 2.0f * u2;
    float x = r * cos(phi);
    float y = r * sin(phi);
    return Direction<S>(x, y, z);
}
template<typename S>
Direction<S> UniformSampleSphere(const float &u1, const float &u2) {
    float z = 1.f - 2.f * u1;
    float r = sqrt(max(0.f, 1.f - z*z));
    float phi = M_PI * 2.0f * u2;
    float x = r * cos(phi);
    float y = r * sin(phi);
    return Direction<S>(x, y, z);
}
mfloat<invsolidangle_d>  UniformHemispherePdf();
mfloat<invsolidangle_d>  UniformSpherePdf();
template<typename S>
Direction<S> UniformSampleCone(float u1, float u2, float costhetamax)
{
    float costheta = (1.f - u1) + u1 * costhetamax;
    float sintheta = sqrt(1.f - costheta*costheta);
    float phi = u2 * circleAngle;
    return Direction<S>(cos(phi) * sintheta, sin(phi) * sintheta, costheta);
}
template<typename S>
Direction<S> UniformSampleCone(float u1, float u2, float costhetamax,
                            const Direction<S> &x, const Direction<S> &y, const Direction<S> &z)
{
    float costheta = Lerp(u1, costhetamax, 1.f);
    float sintheta = sqrtf(1.f - costheta*costheta);
    mfloat<angle_d> phi = u2 * circleAngle;
    return Normalize<S>((cos(phi) * sintheta) * x + (sin(phi) * sintheta) * y +
        (costheta) * z);
}
mfloat<invsolidangle_d> UniformConePdf(float thetamax);
void UniformSampleDisk(float u1, float u2, float *x, float *y);
void ConcentricSampleDisk(float u1, float u2, float *dx, float *dy);
template<typename S>
inline Direction<S> CosineSampleHemisphere(float u1, float u2) {
    float x, y, z;
    ConcentricSampleDisk(u1, u2, &x, &y);
    z = sqrtf(max(0.f, 1.f - x*x - y*y));
    return Direction<S>(x, y, z);
}

inline mfloat<invsolidangle_d> CosineHemispherePdf(const mfloat<proj_d> &costheta, float phi) {
    return costheta / hemisphericalProjectedAngle;
}

void UniformSampleTriangle(float ud1, float ud2, float *u, float *v);
template<typename D>
struct Distribution2D {
    // Distribution2D Public Methods
    Distribution2D(const mfloat<D> *data, int nu, int nv)
    {
        pConditionalV.reserve(nv);
        for (int v = 0; v < nv; ++v) {
            // Compute conditional sampling distribution for $\tilde{v}$
            pConditionalV.push_back(new Distribution1D<D>(&data[v*nu], nu));
        }
        // Compute marginal sampling distribution $p[\tilde{v}]$
        std::vector<mfloat<D> > marginalFunc;
        marginalFunc.reserve(nv);
        for (int v = 0; v < nv; ++v)
            marginalFunc.push_back(pConditionalV[v]->funcInt);
        pMarginal = new Distribution1D<D>(&marginalFunc[0], nv);
    }
    ~Distribution2D() {
            delete pMarginal;
            for (uint32_t i = 0; i < pConditionalV.size(); ++i)
                delete pConditionalV[i];
    }

    void SampleContinuous(float u0, float u1, float uv[2],
                          float *pdf) const {
        float pdfs[2];
        uv[1] = pMarginal->SampleContinuous(u1, &pdfs[1]);
        int v = Clamp(Float2Int(uv[1] * pMarginal->count), 0,
                      pMarginal->count-1);
        uv[0] = pConditionalV[v]->SampleContinuous(u0, &pdfs[0]);
        *pdf = pdfs[0] * pdfs[1];
    }
    float Pdf(float u, float v) const {
        int iu = Clamp(Float2Int(u * pConditionalV[0]->count), 0,
                       pConditionalV[0]->count-1);
        int iv = Clamp(Float2Int(v * pMarginal->count), 0,
                       pMarginal->count-1);
        return (pConditionalV[iv]->func[iu] * pMarginal->func[iv]) /
               (pConditionalV[iv]->funcInt * pMarginal->funcInt);
    }
private:
    // Distribution2D Private Data
    std::vector<Distribution1D<D> *> pConditionalV;
    Distribution1D<D> *pMarginal;
};


void StratifiedSample1D(float *samples, int nsamples, RNG &rng,
                        bool jitter = true);
void StratifiedSample2D(float *samples, int nx, int ny, RNG &rng,
                        bool jitter = true);
template <typename T>
void Shuffle(T *samp, uint32_t count, uint32_t dims, RNG &rng) {
    for (uint32_t i = 0; i < count; ++i) {
        uint32_t other = i + (rng.RandomUInt() % (count - i));
        for (uint32_t j = 0; j < dims; ++j)
            swap(samp[dims*i + j], samp[dims*other + j]);
    }
}


void LatinHypercube(float *samples, uint32_t nSamples, uint32_t nDim, RNG &rng);
inline double RadicalInverse(int n, int base) {
    double val = 0;
    double invBase = 1. / base, invBi = invBase;
    while (n > 0) {
        // Compute next digit of radical inverse
        int d_i = (n % base);
        val += d_i * invBi;
        n *= invBase;
        invBi *= invBase;
    }
    return val;
}


inline void GeneratePermutation(uint32_t *buf, uint32_t b, RNG &rng) {
    for (uint32_t i = 0; i < b; ++i)
        buf[i] = i;
    Shuffle(buf, b, 1, rng);
}


inline double PermutedRadicalInverse(uint32_t n, uint32_t base,
                                     const uint32_t *p) {
    double val = 0;
    double invBase = 1. / base, invBi = invBase;
    while (n > 0) {
        uint32_t d_i = p[n % base];
        val += d_i * invBi;
        n *= invBase;
        invBi *= invBase;
    }
    return val;
}


class PermutedHalton {
public:
    // PermutedHalton Public Methods
    PermutedHalton(uint32_t d, RNG &rng);
    ~PermutedHalton() {
        delete[] b;
        delete[] permute;
    }
    void Sample(uint32_t n, float *out) const {
        uint32_t *p = permute;
        for (uint32_t i = 0; i < dims; ++i) {
            out[i] = PermutedRadicalInverse(n, b[i], p);
            p += b[i];
        }
    }
private:
    // PermutedHalton Private Data
    uint32_t dims;
    uint32_t *b, *permute;
    PermutedHalton(const PermutedHalton &);
    PermutedHalton &operator=(const PermutedHalton &);
};


inline float VanDerCorput(uint32_t n, uint32_t scramble = 0);
inline float Sobol2(uint32_t n, uint32_t scramble = 0);
inline float LarcherPillichshammer2(uint32_t n, uint32_t scramble = 0);
inline void Sample02(uint32_t n, uint32_t scramble[2], float sample[2]);
int LDPixelSamplerealsNeeded(const Sample *sample, int nPixelSamples);
void LDPixelSample(int xPos, int yPos, const mfloat<time_d> &shutterOpen,
    const mfloat<time_d> &shutterClose, int nPixelSamples, Sample *samples, float *buf, RNG &rng);

template<typename S>
Direction<S> SampleHG(const Direction<S> &w, float g, float u1, float u2)
{
    float costheta;
    if (fabsf(g) < 1e-3)
        costheta = 1.f - 2.f * u1;
    else {
        float sqrTerm = (1.f - g * g) /
            (1.f - g + 2.f * g * u1);
        costheta = (1.f + g * g - sqrTerm * sqrTerm) / (2.f * g);
    }
    float sintheta = sqrtf(max(0.f, 1.f-costheta*costheta));
    mfloat<angle_d> phi = circleAngle * u2;
    Direction<S> v1, v2;
    CoordinateSystem(w, &v1, &v2);
    return SphericalDirection(sintheta, costheta, phi, v1, v2, w);
}
template<typename S>
float HGPdf(const Direction<S> &w, const Direction<S> &wp, float g)
{
    return PhaseHG(w, wp, g);
}

// Monte Carlo Inline Functions
inline float BalanceHeuristic(int nf, const mfloat<invsolidangle_d> &fPdf, int ng, const mfloat<invsolidangle_d> &gPdf) {
    return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}


inline float PowerHeuristic(int nf, const mfloat<invsolidangle_d> &fPdf, int ng, const mfloat<invsolidangle_d> &gPdf) {
    mfloat<invsolidangle_d> f = nf * fPdf, g = ng * gPdf;
    return (f*f) / (f*f + g*g);
}



// Sampling Inline Functions
inline void Sample02(uint32_t n, uint32_t scramble[2], float sample[2]) {
    sample[0] = VanDerCorput(n, scramble[0]);
    sample[1] = Sobol2(n, scramble[1]);
}


inline float VanDerCorput(uint32_t n, uint32_t scramble) {
    // Reverse bits of _n_
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    n ^= scramble;
    return (float)n / (float)0x100000000LL;
}


inline float Sobol2(uint32_t n, uint32_t scramble) {
    for (uint32_t v = 1 << 31; n != 0; n >>= 1, v ^= v >> 1)
        if (n & 0x1) scramble ^= v;
    return (float)scramble / (float)0x100000000LL;
}


inline float
LarcherPillichshammer2(uint32_t n, uint32_t scramble) {
    for (uint32_t v = 1 << 31; n != 0; n >>= 1, v |= v >> 1)
        if (n & 0x1) scramble ^= v;
    return (float)scramble / (float)0x100000000LL;
}


inline void LDShuffleScrambled1D(int nSamples, int nPixel,
                                 float *samples, RNG &rng) {
    uint32_t scramble = rng.RandomUInt();
    for (int i = 0; i < nSamples * nPixel; ++i)
        samples[i] = VanDerCorput(i, scramble);
    for (int i = 0; i < nPixel; ++i)
        Shuffle(samples + i * nSamples, nSamples, 1, rng);
    Shuffle(samples, nPixel, nSamples, rng);
}


inline void LDShuffleScrambled2D(int nSamples, int nPixel,
                                 float *samples, RNG &rng) {
    uint32_t scramble[2] = { static_cast<uint32_t>(rng.RandomUInt()), 
        static_cast<uint32_t>(rng.RandomUInt()) };
    for (int i = 0; i < nSamples * nPixel; ++i)
        Sample02(i, scramble, &samples[2*i]);
    for (int i = 0; i < nPixel; ++i)
        Shuffle(samples + 2 * i * nSamples, nSamples, 2, rng);
    Shuffle(samples, nPixel, 2 * nSamples, rng);
}



#endif // PBRT_CORE_MONTECARLO_H
