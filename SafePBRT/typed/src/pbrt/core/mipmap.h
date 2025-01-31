
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

#ifndef PBRT_CORE_MIPMAP_H
#define PBRT_CORE_MIPMAP_H

// core/mipmap.h*
#include "pbrt.h"
#include "spectrum.h"
#include "texture.h"

// MIPMap Declarations
typedef enum {
    TEXTURE_REPEAT,
    TEXTURE_BLACK,
    TEXTURE_CLAMP
} ImageWrap;
template <typename T> class MIPMap {
public:
    // MIPMap Public Methods
    MIPMap() { pyramid = NULL; width = height = nLevels = 0; }
    MIPMap(uint32_t xres, uint32_t yres, const T *data, bool doTri = false,
           float maxAniso = 8.f, ImageWrap wrapMode = TEXTURE_REPEAT);
    ~MIPMap();
    uint32_t Width() const { return width; }
    uint32_t Height() const { return height; }
    uint32_t Levels() const { return nLevels; }
    const T &Texel(uint32_t level, int s, int t) const;
    T Lookup(float s, float t, float width = 0.f) const;
    T Lookup(float s, float t, float ds0, float dt0,
        float ds1, float dt1) const;
private:
    // MIPMap Private Methods
    struct ResampleWeight;
    ResampleWeight *resampleWeights(uint32_t oldres, uint32_t newres) {
        Assert(newres >= oldres);
        ResampleWeight *wt = new ResampleWeight[newres];
        float filterwidth = 2.f;
        for (uint32_t i = 0; i < newres; ++i) {
            // Compute image resampling weights for _i_th texel
            float center = (i + .5f) * oldres / newres;
            wt[i].firstTexel = Floor2Int((center - filterwidth) + 0.5f);
            for (int j = 0; j < 4; ++j) {
                float pos = wt[i].firstTexel + j + .5f;
                wt[i].weight[j] = Lanczos((pos - center) / filterwidth);
            }

            // Normalize filter weights for texel resampling
            float invSumWts = 1.f / (wt[i].weight[0] + wt[i].weight[1] +
                                     wt[i].weight[2] + wt[i].weight[3]);
            for (uint32_t j = 0; j < 4; ++j)
                wt[i].weight[j] *= invSumWts;
        }
        return wt;
    }
    float clamp(float v) { return Clamp(v, 0.f, INFINITY); }
    template<typename D>
    RGBSpectrum<D> clamp(const RGBSpectrum<D> &v) { return v.Clamp(mfloat<D>(0.0f), mfloat<D>(INFINITY)); }
    //template<typename D>
    //SampledSpectrum<D> clamp(const SampledSpectrum<D> &v) { return v.Clamp(0.f, INFINITY); }
    T triangle(uint32_t level, float s, float t) const;
    T EWA(uint32_t level, float s, float t, float ds0, float dt0, float ds1, float dt1) const;

    // MIPMap Private Data
    bool doTrilinear;
    float maxAnisotropy;
    ImageWrap wrapMode;
    struct ResampleWeight {
        int firstTexel;
        float weight[4];
    };
    BlockedArray<T> **pyramid;
    uint32_t width, height, nLevels;
#define WEIGHT_LUT_SIZE 128
    static float *weightLut;
};



// MIPMap Method Definitions
template <typename T>
MIPMap<T>::MIPMap(uint32_t sres, uint32_t tres, const T *img, bool doTri,
                  float maxAniso, ImageWrap wm) {
    doTrilinear = doTri;
    maxAnisotropy = maxAniso;
    wrapMode = wm;
    T *resampledImage = NULL;
    if (!IsPowerOf2(sres) || !IsPowerOf2(tres)) {
        // Resample image to power-of-two resolution
        uint32_t sPow2 = RoundUpPow2(sres), tPow2 = RoundUpPow2(tres);

        // Resample image in $s$ direction
        ResampleWeight *sWeights = resampleWeights(sres, sPow2);
        resampledImage = new T[sPow2 * tPow2];

        // Apply _sWeights_ to zoom in $s$ direction
        for (uint32_t t = 0; t < tres; ++t) {
            for (uint32_t s = 0; s < sPow2; ++s) {
                // Compute texel $(s,t)$ in $s$-zoomed image
                resampledImage[t*sPow2+s] = texval<T>::black();
                for (int j = 0; j < 4; ++j) {
                    int origS = sWeights[s].firstTexel + j;
                    if (wrapMode == TEXTURE_REPEAT)
                        origS = Mod(origS, sres);
                    else if (wrapMode == TEXTURE_CLAMP)
                        origS = Clamp(origS, 0, sres-1);
                    if (origS >= 0 && origS < (int)sres)
                        resampledImage[t*sPow2+s] += sWeights[s].weight[j] *
                                                     img[t*sres + origS];
                }
            }
        }
        delete[] sWeights;

        // Resample image in $t$ direction
        ResampleWeight *tWeights = resampleWeights(tres, tPow2);
        T *workData = new T[tPow2];
        for (uint32_t s = 0; s < sPow2; ++s) {
            for (uint32_t t = 0; t < tPow2; ++t) {
                workData[t] = texval<T>::black();
                for (uint32_t j = 0; j < 4; ++j) {
                    int offset = tWeights[t].firstTexel + j;
                    if (wrapMode == TEXTURE_REPEAT) offset = Mod(offset, tres);
                    else if (wrapMode == TEXTURE_CLAMP) offset = Clamp(offset, 0, tres-1);
                    if (offset >= 0 && offset < (int)tres)
                        workData[t] += tWeights[t].weight[j] *
                            resampledImage[offset*sPow2 + s];
                }
            }
            for (uint32_t t = 0; t < tPow2; ++t)
                resampledImage[t*sPow2 + s] = clamp(workData[t]);
        }
        delete[] workData;
        delete[] tWeights;
        img = resampledImage;
        sres = sPow2;
        tres = tPow2;
    }
    width = sres;
    height = tres;
    // Initialize levels of MIPMap from image
    nLevels = 1 + Log2Int(float(max(sres, tres)));
    pyramid = new BlockedArray<T> *[nLevels];

    // Initialize most detailed level of MIPMap
    pyramid[0] = new BlockedArray<T>(sres, tres, img);
    for (uint32_t i = 1; i < nLevels; ++i) {
        // Initialize $i$th MIPMap level from $i-1$st level
        uint32_t sRes = max(1u, pyramid[i-1]->uSize()/2);
        uint32_t tRes = max(1u, pyramid[i-1]->vSize()/2);
        pyramid[i] = new BlockedArray<T>(sRes, tRes);

        // Filter four texels from finer level of pyramid
        for (uint32_t t = 0; t < tRes; ++t)
            for (uint32_t s = 0; s < sRes; ++s)
                (*pyramid[i])(s, t) = .25f *
                   (Texel(i-1, 2*s, 2*t)   + Texel(i-1, 2*s+1, 2*t) +
                    Texel(i-1, 2*s, 2*t+1) + Texel(i-1, 2*s+1, 2*t+1));
    }
    if (resampledImage) delete[] resampledImage;
    // Initialize EWA filter weights if needed
    if (!weightLut) {
        weightLut = AllocAligned<float>(WEIGHT_LUT_SIZE);
        for (int i = 0; i < WEIGHT_LUT_SIZE; ++i) {
            float alpha = 2;
            float r2 = float(i) / float(WEIGHT_LUT_SIZE - 1);
            weightLut[i] = expf(-alpha * r2) - expf(-alpha);
        }
    }
}


template <typename T>
const T &MIPMap<T>::Texel(uint32_t level, int s, int t) const {
    Assert(level < nLevels);
    const BlockedArray<T> &l = *pyramid[level];
    // Compute texel $(s,t)$ accounting for boundary conditions
    switch (wrapMode) {
        case TEXTURE_REPEAT:
            s = Mod(s, l.uSize());
            t = Mod(t, l.vSize());
            break;
        case TEXTURE_CLAMP:
            s = Clamp(s, 0, l.uSize() - 1);
            t = Clamp(t, 0, l.vSize() - 1);
            break;
        case TEXTURE_BLACK: {
            static const T black = texval<T>::black();
            if (s < 0 || s >= (int)l.uSize() ||
                t < 0 || t >= (int)l.vSize())
                return black;
            break;
        }
    }
    PBRT_ACCESSED_TEXEL(const_cast<MIPMap<T> *>(this), level, s, t);
    return l(s, t);
}


template <typename T>
MIPMap<T>::~MIPMap() {
    for (uint32_t i = 0; i < nLevels; ++i)
        delete pyramid[i];
    delete[] pyramid;
}


template <typename T>
T MIPMap<T>::Lookup(float s, float t, float width) const {
    // Compute MIPMap level for trilinear filtering
    float level = nLevels - 1 + Log2(max(width, 1e-8f));

    // Perform trilinear interpolation at appropriate MIPMap level
    PBRT_MIPMAP_TRILINEAR_FILTER(const_cast<MIPMap<T> *>(this), s, t, width, level, nLevels);
    if (level < 0)
        return triangle(0, s, t);
    else if (level >= nLevels - 1)
        return Texel(nLevels-1, 0, 0);
    else {
        uint32_t iLevel = Floor2Int(level);
        float delta = level - iLevel;
        return (1.f-delta) * triangle(iLevel, s, t) +
               delta * triangle(iLevel+1, s, t);
    }
}


template <typename T>
T MIPMap<T>::triangle(uint32_t level, float s, float t) const {
    level = Clamp(level, 0, nLevels-1);
    s = s * pyramid[level]->uSize() - 0.5f;
    t = t * pyramid[level]->vSize() - 0.5f;
    int s0 = Floor2Int(s), t0 = Floor2Int(t);
    float ds = s - s0, dt = t - t0;
    return (1.f-ds) * (1.f-dt) * Texel(level, s0, t0) +
           (1.f-ds) * dt       * Texel(level, s0, t0+1) +
           ds       * (1.f-dt) * Texel(level, s0+1, t0) +
           ds       * dt       * Texel(level, s0+1, t0+1);
}


template <typename T>
T MIPMap<T>::Lookup(float s, float t, float ds0, float dt0,
                    float ds1, float dt1) const {
    if (doTrilinear) {
        PBRT_STARTED_TRILINEAR_TEXTURE_LOOKUP(s, t);
        T val = Lookup(s, t,
                       2.f * max(max(fabsf(ds0), fabsf(dt0)),
                                 max(fabsf(ds1), fabsf(dt1))));
        PBRT_FINISHED_TRILINEAR_TEXTURE_LOOKUP();
        return val;
    }
    PBRT_STARTED_EWA_TEXTURE_LOOKUP(s, t);
    // Compute ellipse minor and major axes
    if (ds0*ds0 + dt0*dt0 < ds1*ds1 + dt1*dt1) {
        swap(ds0, ds1);
        swap(dt0, dt1);
    }
    float majorLength = sqrtf(ds0*ds0 + dt0*dt0);
    float minorLength = sqrtf(ds1*ds1 + dt1*dt1);

    // Clamp ellipse eccentricity if too large
    if (minorLength * maxAnisotropy < majorLength) {
        float scale = majorLength / (minorLength * maxAnisotropy);
        ds1 *= scale;
        dt1 *= scale;
        minorLength *= scale;
    }
    if (minorLength == 0.f)
       return triangle(0, s, t);

    // Choose level of detail for EWA lookup and perform EWA filtering
    float lod = max(0.f, nLevels - 1.f + Log2(minorLength));
    uint32_t ilod = Floor2Int(lod);
    PBRT_MIPMAP_EWA_FILTER(const_cast<MIPMap<T> *>(this), s, t, ds0, ds1, dt0, dt1, minorLength, majorLength, lod, nLevels);
    float d = lod - ilod;
    T val = (1.f - d) * EWA(ilod,   s, t, ds0, dt0, ds1, dt1) +
                   d  * EWA(ilod+1, s, t, ds0, dt0, ds1, dt1);
    PBRT_FINISHED_EWA_TEXTURE_LOOKUP();
    return val;
}


template <typename T>
T MIPMap<T>::EWA(uint32_t level, float s, float t, float ds0, float dt0,
                 float ds1, float dt1) const {
    if (level >= nLevels) return Texel(nLevels-1, 0, 0);
    // Convert EWA coordinates to appropriate scale for level
    s = s * pyramid[level]->uSize() - 0.5f;
    t = t * pyramid[level]->vSize() - 0.5f;
    ds0 *= pyramid[level]->uSize();
    dt0 *= pyramid[level]->vSize();
    ds1 *= pyramid[level]->uSize();
    dt1 *= pyramid[level]->vSize();

    // Compute ellipse coefficients to bound EWA filter region
    auto A = dt0*dt0 + dt1*dt1 + 1;
    auto B = -2.f * (ds0*dt0 + ds1*dt1);
    auto C = ds0*ds0 + ds1*ds1 + 1;
    auto invF = 1.f / (A*C - B*B*0.25f);
    A *= invF;
    B *= invF;
    C *= invF;

    // Compute the ellipse's $(s,t)$ bounding box in texture space
    auto det = -B*B + 4.f*A*C;
    auto invDet = 1.f / det;
    auto uSqrt = sqrtf(det * C), vSqrt = sqrtf(A * det);
    int s0 = Ceil2Int (s - 2.f * invDet * uSqrt);
    int s1 = Floor2Int(s + 2.f * invDet * uSqrt);
    int t0 = Ceil2Int (t - 2.f * invDet * vSqrt);
    int t1 = Floor2Int(t + 2.f * invDet * vSqrt);

    // Scan over ellipse bound and compute quadratic equation
    T sum(0.);
    float sumWts = 0.f;
    for (int it = t0; it <= t1; ++it) {
        float tt = it - t;
        for (int is = s0; is <= s1; ++is) {
            float ss = is - s;
            // Compute squared radius and filter texel if inside ellipse
            float r2 = A*ss*ss + B*ss*tt + C*tt*tt;
            if (r2 < 1.) {
                float weight = weightLut[min(Float2Int(r2 * WEIGHT_LUT_SIZE),
                                             WEIGHT_LUT_SIZE-1)];
                sum += Texel(level, is, it) * weight;
                sumWts += weight;
            }
        }
    }
    return sum / sumWts;
}


template <typename T> float *MIPMap<T>::weightLut = NULL;

#endif // PBRT_CORE_MIPMAP_H
