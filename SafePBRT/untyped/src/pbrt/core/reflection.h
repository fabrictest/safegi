
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

#ifndef PBRT_CORE_REFLECTION_H
#define PBRT_CORE_REFLECTION_H

// core/reflection.h*
#include "pbrt.h"
#include "geometry.h"
#include "shape.h"
#include "rng.h"
#include "spectrum.h"
#include "kdtree.h"

// Reflection Declarations
Spectrum FrDiel(float cosi, float cost, const Spectrum &etai,
                const Spectrum &etat);
Spectrum FrCond(float cosi, const Spectrum &n, const Spectrum &k);

// BSDF Inline Functions
inline float CosTheta(const Vector &w) { return w.z; }
inline float AbsCosTheta(const Vector &w) { return fabsf(w.z); }
inline float SinTheta2(const Vector &w) {
    return 1.f - CosTheta(w)*CosTheta(w);
}


inline float SinTheta(const Vector &w) {
    return sqrtf(SinTheta2(w));
}


inline float CosPhi(const Vector &w) {
    return w.x / SinTheta(w);
}


inline float SinPhi(const Vector &w) {
    return w.y / SinTheta(w);
}


inline bool SameHemisphere(const Vector &w, const Vector &wp) {
    return w.z * wp.z > 0.f;
}



// BSDF Declarations
enum BxDFType {
    BSDF_REFLECTION   = 1<<0,
    BSDF_TRANSMISSION = 1<<1,
    BSDF_DIFFUSE      = 1<<2,
    BSDF_GLOSSY       = 1<<3,
    BSDF_SPECULAR     = 1<<4,
    BSDF_ALL_TYPES        = BSDF_DIFFUSE |
                            BSDF_GLOSSY |
                            BSDF_SPECULAR,
    BSDF_ALL_REFLECTION   = BSDF_REFLECTION |
                            BSDF_ALL_TYPES,
    BSDF_ALL_TRANSMISSION = BSDF_TRANSMISSION |
                            BSDF_ALL_TYPES,
    BSDF_ALL              = BSDF_ALL_REFLECTION |
                            BSDF_ALL_TRANSMISSION
};


struct BSDFSample {
   // BSDFSample Public Methods
   BSDFSample(float up0, float up1, float ucomp) {
       uDir[0] = up0;
       uDir[1] = up1;
       uComponent = ucomp;
   }
   BSDFSample(RNG &rng) {
      uDir[0] = rng.RandomFloat();
      uDir[1] = rng.RandomFloat();
      uComponent = rng.RandomFloat();
   }
   BSDFSample(const Sample *sample, const BSDFSampleOffsets &offsets, uint32_t num);
   BSDFSample() { }
   float uDir[2], uComponent;
};


struct BSDFSampleOffsets {
    BSDFSampleOffsets() { }
    BSDFSampleOffsets(int count, Sample *sample);
    int nSamples, componentOffset, dirOffset;
};


class BSDF {
public:
    // BSDF Public Methods
    Spectrum Sample_f(const Vector &wo, Vector *wi, const BSDFSample &bsdfSample,
                      float *pdf, BxDFType flags = BSDF_ALL,
                      BxDFType *sampledType = NULL) const;
    float Pdf(const Vector &wo, const Vector &wi,
              BxDFType flags = BSDF_ALL) const;
    BSDF(const DifferentialGeometry &dgs, const Normal &ngeom,
         float eta = 1.f);
    inline void Add(BxDF *bxdf);
    int NumComponents() const { return nBxDFs; }
    int NumComponents(BxDFType flags) const;
    Vector WorldToLocal(const Vector &v) const {
        return Vector(Dot(v, sn), Dot(v, tn), Dot(v, nn));
    }
    Vector LocalToWorld(const Vector &v) const {
        return Vector(sn.x * v.x + tn.x * v.y + nn.x * v.z,
                      sn.y * v.x + tn.y * v.y + nn.y * v.z,
                      sn.z * v.x + tn.z * v.y + nn.z * v.z);
    }
    Spectrum f(const Vector &woW, const Vector &wiW, BxDFType flags = BSDF_ALL) const;
    Spectrum rho(RNG &rng, BxDFType flags = BSDF_ALL,
                 int sqrtSamples = 6) const;
    Spectrum rho(const Vector &wo, RNG &rng, BxDFType flags = BSDF_ALL,
                 int sqrtSamples = 6) const;

    // BSDF Public Data
    const DifferentialGeometry dgShading;
    const float eta;
private:
    // BSDF Private Methods
    ~BSDF() { }

    // BSDF Private Data
    Normal nn, ng;
    Vector sn, tn;
    int nBxDFs;
#define MAX_BxDFS 8
    BxDF *bxdfs[MAX_BxDFS];
    friend class MixMaterial;
};


#define BSDF_ALLOC(arena, Type) new (arena.Alloc(sizeof(Type))) Type

// BxDF Declarations
class BxDF {
public:
    // BxDF Interface
    virtual ~BxDF() { }
    BxDF(BxDFType t) : type(t) { }
    bool MatchesFlags(BxDFType flags) const {
        return (type & flags) == type;
    }
    virtual Spectrum f(const Vector &wo, const Vector &wi) const = 0;
    virtual Spectrum Sample_f(const Vector &wo, Vector *wi,
                              float u1, float u2, float *pdf) const;
    virtual Spectrum rho(const Vector &wo, int nSamples,
                         const float *samples) const;
    virtual Spectrum rho(int nSamples, const float *samples1,
                         const float *samples2) const;
    virtual float Pdf(const Vector &wi, const Vector &wo) const;

    // BxDF Public Data
    const BxDFType type;
};


class BRDFToBTDF : public BxDF {
public:
    // BRDFToBTDF Public Methods
    BRDFToBTDF(BxDF *b)
        : BxDF(BxDFType(b->type ^ (BSDF_REFLECTION | BSDF_TRANSMISSION))) {
        brdf = b;
    }
    static Vector otherHemisphere(const Vector &w) {
        return Vector(w.x, w.y, -w.z);
    }
    Spectrum f(const Vector &wo, const Vector &wi) const;
    Spectrum Sample_f(const Vector &wo, Vector *wi, float u1, float u2,
                      float *pdf) const;
    Spectrum rho(const Vector &w, int nSamples, const float *samples) const {
        return brdf->rho(otherHemisphere(w), nSamples, samples);
    }
    Spectrum rho(int nSamples, const float *samples1, const float *samples2) const {
        return brdf->rho(nSamples, samples1, samples2);
    }
    float Pdf(const Vector &wo, const Vector &wi) const;
private:
    BxDF *brdf;
};


class Fresnel {
public:
    // Fresnel Interface
    virtual ~Fresnel();
    virtual Spectrum Evaluate(float cosi) const = 0;
};


class FresnelConductor : public Fresnel {
public:
    // FresnelConductor Public Methods
    Spectrum Evaluate(float cosi) const;
    FresnelConductor(const Spectrum &e, const Spectrum &kk)
        : eta(e), k(kk) {
    }
private:
    Spectrum eta, k;
};


class FresnelDielectric : public Fresnel {
public:
    // FresnelDielectric Public Methods
    Spectrum Evaluate(float cosi) const;
    FresnelDielectric(float ei, float et) : eta_i(ei), eta_t(et) { }
private:
    float eta_i, eta_t;
};


class FresnelNoOp : public Fresnel {
public:
    Spectrum Evaluate(float) const { return Spectrum(1.); }
};


class SpecularReflection : public BxDF {
public:
    // SpecularReflection Public Methods
    SpecularReflection(const Spectrum &r, Fresnel *f)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)),
          R(r), fresnel(f) {
    }
    Spectrum f(const Vector &, const Vector &) const {
        return Spectrum(0.);
    }
    Spectrum Sample_f(const Vector &wo, Vector *wi,
                      float u1, float u2, float *pdf) const;
    float Pdf(const Vector &wo, const Vector &wi) const {
        return 0.;
    }
private:
    // SpecularReflection Private Data
    Spectrum R;
    Fresnel *fresnel;
};


class SpecularTransmission : public BxDF {
public:
    // SpecularTransmission Public Methods
    SpecularTransmission(const Spectrum &t, float ei, float et)
        : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)),
          fresnel(ei, et) {
        T = t;
        etai = ei;
        etat = et;
    }
    Spectrum f(const Vector &, const Vector &) const {
        return Spectrum(0.);
    }
    Spectrum Sample_f(const Vector &wo, Vector *wi, float u1, float u2, float *pdf) const;
    float Pdf(const Vector &wo, const Vector &wi) const {
        return 0.;
    }
private:
    // SpecularTransmission Private Data
    Spectrum T;
    float etai, etat;
    FresnelDielectric fresnel;
};


class Lambertian : public BxDF {
public:
    // Lambertian Public Methods
    Lambertian(const Spectrum &reflectance)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(reflectance) { }
    Spectrum f(const Vector &wo, const Vector &wi) const;
    Spectrum rho(const Vector &, int, const float *) const { return R; }
    Spectrum rho(int, const float *, const float *) const { return R; }
private:
    // Lambertian Private Data
    Spectrum R;
};


class OrenNayar : public BxDF {
public:
    // OrenNayar Public Methods
    Spectrum f(const Vector &wo, const Vector &wi) const;
    OrenNayar(const Spectrum &reflectance, float sig)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
          R(reflectance) {
        float sigma = Radians(sig);
        float sigma2 = sigma*sigma;
        A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
        B = 0.45f * sigma2 / (sigma2 + 0.09f);
    }
private:
    // OrenNayar Private Data
    Spectrum R;
    float A, B;
};


class MicrofacetDistribution {
public:
    // MicrofacetDistribution Interface
    virtual ~MicrofacetDistribution() { }
    virtual float D(const Vector &wh) const = 0;
    virtual void Sample_f(const Vector &wo, Vector *wi,
                          float u1, float u2, float *pdf) const = 0;
    virtual float Pdf(const Vector &wo, const Vector &wi) const = 0;
};


class Microfacet : public BxDF {
public:
    // Microfacet Public Methods
    Microfacet(const Spectrum &reflectance, Fresnel *f,
        MicrofacetDistribution *d);
    Spectrum f(const Vector &wo, const Vector &wi) const;
    float G(const Vector &wo, const Vector &wi, const Vector &wh) const {
        float NdotWh = AbsCosTheta(wh);
        float NdotWo = AbsCosTheta(wo);
        float NdotWi = AbsCosTheta(wi);
        float WOdotWh = AbsDot(wo, wh);
        return min(1.f, min((2.f * NdotWh * NdotWo / WOdotWh),
                            (2.f * NdotWh * NdotWi / WOdotWh)));
    }
    Spectrum Sample_f(const Vector &wo, Vector *wi,
                              float u1, float u2, float *pdf) const;
    float Pdf(const Vector &wo, const Vector &wi) const;
private:
    // Microfacet Private Data
    Spectrum R;
    MicrofacetDistribution *distribution;
    Fresnel *fresnel;
};


class Blinn : public MicrofacetDistribution {
public:
    Blinn(float e) { if (e > 10000.f || isnan(e)) e = 10000.f;
                     exponent = e; }
    // Blinn Public Methods
    float D(const Vector &wh) const {
        float costhetah = AbsCosTheta(wh);
        return (exponent+2) * INV_TWOPI * powf(costhetah, exponent);
    }
    virtual void Sample_f(const Vector &wi, Vector *sampled_f, float u1, float u2, float *pdf) const;
    virtual float Pdf(const Vector &wi, const Vector &wo) const;
private:
    float exponent;
};


class Anisotropic : public MicrofacetDistribution {
public:
    // Anisotropic Public Methods
    Anisotropic(float x, float y) {
        ex = x; ey = y;
        if (ex > 10000.f || isnan(ex)) ex = 10000.f;
        if (ey > 10000.f || isnan(ey)) ey = 10000.f;
    }
    float D(const Vector &wh) const {
        float costhetah = AbsCosTheta(wh);
        float d = 1.f - costhetah * costhetah;
        if (d == 0.f) return 0.f;
        float e = (ex * wh.x * wh.x + ey * wh.y * wh.y) / d;
        return sqrtf((ex+2.f) * (ey+2.f)) * INV_TWOPI * powf(costhetah, e);
    }
    void Sample_f(const Vector &wo, Vector *wi, float u1, float u2, float *pdf) const;
    float Pdf(const Vector &wo, const Vector &wi) const;
    void sampleFirstQuadrant(float u1, float u2, float *phi, float *costheta) const;
private:
    float ex, ey;
};

// BSDF Inline Method Definitions
inline void BSDF::Add(BxDF *b) {
    Assert(nBxDFs < MAX_BxDFS);
    bxdfs[nBxDFs++] = b;
}


inline int BSDF::NumComponents(BxDFType flags) const {
    int num = 0;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags)) ++num;
    return num;
}



#endif // PBRT_CORE_REFLECTION_H
