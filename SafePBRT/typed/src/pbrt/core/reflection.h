
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
Spectrum<rho_d> FrDiel(float cosi, float cost, const float &etai,
                const float &etat);
Spectrum<rho_d> FrCond(float cosi, const Spectrum<rho_d> &n, const Spectrum<rho_d> &k);


// BSDF Inline Functions
inline float AbsCosTheta(const Direction<local_s> &w) 
{ 
    return (fabsf(w.z)); 
}

inline float CosTheta(const Direction<local_s> &w) 
{ 
    return (w.z); 
}

inline mfloat<proj_d> ProjectTheta(const Direction<local_s> &w) 
{ 
    return mfloat<proj_d>(CosTheta(w)); 
}

inline mfloat<proj_d> AbsProjectTheta(const Direction<local_s> &w) 
{ 
    return mfloat<proj_d>(AbsCosTheta(w)); 
}

inline float SinTheta2(const Direction<local_s> &w) 
{
    return  1.0f - CosTheta(w) * CosTheta(w);
}

inline float SinTheta(const Direction<local_s> &w) {
    return sqrt(SinTheta2(w));
}

inline float CosPhi(const Direction<local_s> &w) {
    return w.x / SinTheta(w);
}

inline float SinPhi(const Direction<local_s> &w) {
    return w.y / SinTheta(w);
}

inline bool SameHemisphere(const Direction<local_s> &w, const Direction<local_s> &wp) {
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
    Spectrum<brdf_d> Sample_f(const Direction<world_s> &wo, Direction<world_s> *wi, const BSDFSample &bsdfSample,
                      mfloat<invsolidangle_d> *pdf, BxDFType flags = BSDF_ALL,
                      BxDFType *sampledType = NULL) const;
    mfloat<invsolidangle_d> Pdf(const Direction<world_s> &wo, const Direction<world_s> &wi,
              BxDFType flags = BSDF_ALL) const;
    BSDF(const DifferentialGeometry &dgs, const Normal<world_s> &ngeom,
         float eta = 1.f);
    inline void Add(BxDF *bxdf);
    int NumComponents() const { return nBxDFs; }
    int NumComponents(BxDFType flags) const;
    Direction<local_s> WorldToLocal(const Direction<world_s> &v) const {
        return Direction<local_s>(Dot(v, sn), Dot(v, tn), Dot(v, nn));
    }
    Direction<world_s> LocalToWorld(const Direction<local_s> &v) const {
        return Direction<world_s>(sn.x * v.x + tn.x * v.y + nn.x * v.z,
                      sn.y * v.x + tn.y * v.y + nn.y * v.z,
                      sn.z * v.x + tn.z * v.y + nn.z * v.z);
    }
    Spectrum<brdf_d> f(const Direction<world_s> &woW, const Direction<world_s> &wiW, BxDFType flags = BSDF_ALL) const;
    Spectrum<rho_d> rho(RNG &rng, BxDFType flags = BSDF_ALL,
                 int sqrtSamples = 6) const;
    Spectrum<rho_d> rho(const Direction<local_s> &wo, RNG &rng, BxDFType flags = BSDF_ALL,
                 int sqrtSamples = 6) const;

    // BSDF Public Data
    const DifferentialGeometry dgShading;
    const float eta;
private:
    // BSDF Private Methods
    ~BSDF() { }

    // BSDF Private Data
    Normal<world_s> nn, ng;
    Direction<world_s> sn, tn;
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
    virtual Spectrum<brdf_d> f(const Direction<local_s> &wo, const Direction<local_s> &wi) const = 0;
    virtual Spectrum<brdf_d> Sample_f(const Direction<local_s> &wo, Direction<local_s> *wi,
                              float u1, float u2, mfloat<invsolidangle_d> *pdf) const;
    virtual Spectrum<rho_d> rho(const Direction<local_s> &wo, int nSamples,
                         const float *samples) const;
    virtual Spectrum<rho_d> rho(int nSamples, const float *samples1,
                         const float *samples2) const;
    virtual mfloat<invsolidangle_d> Pdf(const Direction<local_s> &wi, const Direction<local_s> &wo) const;

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
    template<typename S>
    static Direction<S> otherHemisphere(const Direction<S> &w) {
        return Direction<S>(w.x, w.y, -w.z);
    }
    Spectrum<brdf_d> f(const Direction<local_s> &wo, const Direction<local_s> &wi) const;
    Spectrum<brdf_d> Sample_f(const Direction<local_s> &wo, Direction<local_s> *wi, 
        float u1, float u2, mfloat<invsolidangle_d> *pdf) const;
    Spectrum<rho_d> rho(const Direction<local_s> &w, int nSamples, 
        const float *samples) const {
        return brdf->rho(otherHemisphere(w), nSamples, samples);
    }
    Spectrum<rho_d> rho(int nSamples, const float *samples1, const float *samples2) const {
        return brdf->rho(nSamples, samples1, samples2);
    }
    mfloat<invsolidangle_d> Pdf(const Direction<local_s> &wo, 
        const Direction<local_s> &wi) const;
private:
    BxDF *brdf;
};


class Fresnel {
public:
    // Fresnel Interface
    virtual ~Fresnel();
    virtual Spectrum<rho_d> Evaluate(float cosi) const = 0;
};


class FresnelConductor : public Fresnel {
public:
    // FresnelConductor Public Methods
    Spectrum<rho_d> Evaluate(float cosi) const;
    FresnelConductor(const Spectrum<rho_d> &e, const Spectrum<rho_d> &kk)
        : eta(e), k(kk) {
    }
private:
    Spectrum<rho_d> eta, k;
};


class FresnelDielectric : public Fresnel {
public:
    // FresnelDielectric Public Methods
    Spectrum<rho_d> Evaluate(float cosi) const;
    FresnelDielectric(const float ei, const float et) : eta_i(ei), eta_t(et) { }
private:
    float eta_i, eta_t;
};


class FresnelNoOp : public Fresnel {
public:
    Spectrum<rho_d> Evaluate(float cosi) const { 
        return Spectrum<rho_d>(mfloat<rho_d>(1.0f)); 
    }
};


class SpecularReflection : public BxDF {
public:
    // SpecularReflection Public Methods
    SpecularReflection(const Spectrum<rho_d> &r, Fresnel *f)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)),
          R(r), fresnel(f) {
    }
    Spectrum<brdf_d> f(const Direction<local_s> &, const Direction<local_s> &) const {
        return Spectrum<brdf_d>();
    }
    Spectrum<brdf_d> Sample_f(const Direction<local_s> &wo, Direction<local_s> *wi,
                      float u1, float u2, mfloat<invsolidangle_d> *pdf) const;
    mfloat<invsolidangle_d> Pdf(const Direction<local_s> &wo, 
        const Direction<local_s> &wi) const {
        return mfloat<invsolidangle_d>(0.0f);
    }
private:
    // SpecularReflection Private Data
    Spectrum<rho_d> R;
    Fresnel *fresnel;
};


class SpecularTransmission : public BxDF {
public:
    // SpecularTransmission Public Methods
    SpecularTransmission(const Spectrum<rho_d> &t, float ei, float et)
        : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)),
          fresnel(ei, et) {
        T = t;
        etai = ei;
        etat = et;
    }
    Spectrum<brdf_d> f(const Direction<local_s> &, const Direction<local_s> &) const {
        return Spectrum<brdf_d>();
    }
    Spectrum<brdf_d> Sample_f(const Direction<local_s> &wo, Direction<local_s> *wi, 
        float u1, float u2, mfloat<invsolidangle_d> *pdf) const;
    mfloat<invsolidangle_d>  Pdf(const Direction<local_s>  
        &wo, const Direction<local_s>  &wi) const {
        return mfloat<invsolidangle_d> (0.0f);
    }
private:
    // SpecularTransmission Private Data
    Spectrum<rho_d> T;
    float etai, etat;
    FresnelDielectric fresnel;
};


class Lambertian : public BxDF {
public:
    // Lambertian Public Methods
    Lambertian(const Spectrum<rho_d> &reflectance)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(reflectance) { }
    Spectrum<brdf_d> f(const Direction<local_s> &wo, const Direction<local_s> &wi) const;
    Spectrum<rho_d> rho(const Direction<local_s> &, int, const float *) const { return R; }
    Spectrum<rho_d> rho(int, const float *, const float *) const { return R; }
private:
    // Lambertian Private Data
    Spectrum<rho_d> R;
};


class OrenNayar : public BxDF {
public:
    // OrenNayar Public Methods
    Spectrum<brdf_d> f(const Direction<local_s> &wo, const Direction<local_s> &wi) const;
    OrenNayar(const Spectrum<rho_d> &reflectance, float sig)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
          R(reflectance) {
        float sigma = _Radians(sig);
        float sigma2 = sigma*sigma;
        A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
        B = 0.45f * sigma2 / (sigma2 + 0.09f);
    }
private:
    // OrenNayar Private Data
    Spectrum<rho_d> R;
    float A;
    float B;
};


class MicrofacetDistribution {
public:
    // MicrofacetDistribution Interface
    virtual ~MicrofacetDistribution() { }
    virtual mfloat<invsolidangle_d> D(const Direction<local_s> &wh) const = 0;
    virtual void Sample_f(const Direction<local_s> &wo, Direction<local_s> *wi,
                          float u1, float u2, mfloat<invsolidangle_d> *pdf) const = 0;
    virtual mfloat<invsolidangle_d> Pdf(const Direction<local_s> &wo, const Direction<local_s> &wi) const = 0;
};


class Microfacet : public BxDF {
public:
    // Microfacet Public Methods
    Microfacet(const Spectrum<rho_d> &reflectance, Fresnel *f,
        MicrofacetDistribution *d);
    Spectrum<brdf_d> f(const Direction<local_s> &wo, const Direction<local_s> &wi) const;
    float G(const Direction<local_s> &wo, const Direction<local_s> &wi, const Direction<local_s> &wh) const {
        float NdotWh = AbsCosTheta(wh);
        float NdotWo = AbsCosTheta(wo);
        float NdotWi = AbsCosTheta(wi);
        float WOdotWh = abs(Dot(wo, wh));
        return min((1.f), min((2.f * NdotWh * NdotWo / WOdotWh),
                            (2.f * NdotWh * NdotWi / WOdotWh)));
    }
    Spectrum<brdf_d> Sample_f(const Direction<local_s> &wo, Direction<local_s> *wi,
                              float u1, float u2, mfloat<invsolidangle_d> *pdf) const;
    mfloat<invsolidangle_d> Pdf(const Direction<local_s> &wo, const Direction<local_s> &wi) const;
private:
    // Microfacet Private Data
    Spectrum<rho_d> R;
    MicrofacetDistribution *distribution;
    Fresnel *fresnel;
};


class Blinn : public MicrofacetDistribution {
public:
    Blinn(float e) { if (e > 10000.f || isnan(e)) e = 10000.f;
                     exponent = e; }
    // Blinn Public Methods
    mfloat<invsolidangle_d> D(const Direction<local_s> &wh) const {
        float costhetah = AbsCosTheta(wh);
        return (exponent+2) / hemisphericalAngle * powf(costhetah, exponent);
    }
    virtual void Sample_f(const Direction<local_s> &wi, Direction<local_s> *sampled_f, float u1, float u2, mfloat<invsolidangle_d> *pdf) const;
    virtual mfloat<invsolidangle_d> Pdf(const Direction<local_s> &wi, const Direction<local_s> &wo) const;
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
    mfloat<invsolidangle_d> D(const Direction<local_s> &wh) const {
        float costhetah = AbsCosTheta(wh);
        float d = 1.f - costhetah * costhetah;
        if (d == 0.f) return mfloat<invsolidangle_d>(0.0f);
        float e = (ex * wh.x * wh.x + ey * wh.y * wh.y) / d;
        return sqrtf((ex+2.f) * (ey+2.f)) / hemisphericalAngle * powf(costhetah, e);
    }
    void Sample_f(const Direction<local_s> &wi, Direction<local_s> *sampled_f, float u1, float u2, mfloat<invsolidangle_d> *pdf) const;
    mfloat<invsolidangle_d> Pdf(const Direction<local_s> &wi, const Direction<local_s> &wo) const;
    void sampleFirstQuadrant(float u1, float u2, mfloat<angle_d> *phi, float *costheta) const;
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
