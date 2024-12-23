
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


// core/reflection.cpp*
#include "reflection.h"
#include "spectrum.h"
#include "sampler.h"
#include "montecarlo.h"
#include <stdarg.h>


 //BxDF Utility Functions
Spectrum<rho_d> FrDiel(float cosi, float cost, const float &etai,
                const float &etat) {
    Spectrum<rho_d> Rparl = Spectrum<rho_d>(((etat * cosi) - (etai * cost)) /
                     ((etat * cosi) + (etai * cost)));
    Spectrum<rho_d> Rperp = Spectrum<rho_d>(((etai * cosi) - (etat * cost)) /
                     ((etai * cosi) + (etat * cost)));
    return (Rparl*Rparl + Rperp*Rperp) / 2.f;
}


Spectrum<rho_d> FrCond(float cosi, const Spectrum<unit_d> &eta, const Spectrum<unit_d> &k) {
    Spectrum<rho_d> tmp = (eta*eta + k*k) * cosi*cosi;
    Spectrum<rho_d> Rparl2 = (tmp - (2.f * eta * cosi) + Spectrum<rho_d>(1.0f)) /
                      (tmp + (2.f * eta * cosi) + Spectrum<rho_d>(1.0f));
    Spectrum<rho_d> tmp_f = eta*eta + k*k;
    Spectrum<rho_d> Rperp2 =
        (tmp_f - (2.f * eta * cosi) + Spectrum<rho_d>(cosi*cosi)) /
        (tmp_f + (2.f * eta * cosi) + Spectrum<rho_d>(cosi*cosi));
    return (Rparl2 + Rperp2) / 2.f;
}



// BxDF Method Definitions
Spectrum<brdf_d> BRDFToBTDF::f(const Direction<local_s> &wo, const Direction<local_s> &wi) const {
    return brdf->f(wo, otherHemisphere(wi));
}


Spectrum<brdf_d> BRDFToBTDF::Sample_f(const Direction<local_s> &wo, Direction<local_s> *wi,
                              float u1, float u2, mfloat<invsolidangle_d> *pdf) const {
    Spectrum<brdf_d> f = brdf->Sample_f(wo, wi, u1, u2, pdf);
    *wi = otherHemisphere(*wi);
    return f;
}

Fresnel::~Fresnel() { }

Spectrum<rho_d> FresnelConductor::Evaluate(float cosi) const {
    return FrCond(fabsf(cosi), eta, k);
}


Spectrum<rho_d> FresnelDielectric::Evaluate(float cosi) const {
    //Compute Fresnel reflectance for dielectric
    cosi = Clamp(cosi, -1.f, 1.f);

    //Compute indices of refraction for dielectric
    bool entering = cosi > 0.;
    float ei = eta_i, et = eta_t;
    if (!entering)
        swap(ei, et);

    //Compute _sint_ using Snell's law
    float sint = ei/et * sqrtf(max(0.f, 1.f - cosi*cosi));
    if (sint >= 1.) {
        //Handle total internal reflection
        return Spectrum<rho_d>(1.0f);
    }
    else {
        float cost = sqrtf(max(0.f, 1.f - sint*sint));
        return FrDiel(fabsf(cosi), cost, ei, et);
    }
}

Spectrum<brdf_d> SpecularReflection::Sample_f(const Direction<local_s> &wo,
        Direction<local_s> *wi, float u1, float u2, mfloat<invsolidangle_d> *pdf) const {
    //Compute perfect specular reflection direction
    *wi = Direction<local_s>(-wo.x, -wo.y, wo.z);
    *pdf = mfloat<invsolidangle_d>(1.0f);
    return fresnel->Evaluate(wo.z) * R / AbsProjectTheta(*wi) / steradians;
}


Spectrum<brdf_d> SpecularTransmission::Sample_f(const Direction<local_s> &wo,
        Direction<local_s> *wi, float u1, float u2, mfloat<invsolidangle_d> *pdf) const {
    // Figure out which $\eta$ is incident and which is transmitted
    bool entering = CosTheta(wo) > 0.;
    float ei = etai, et = etat;
    if (!entering)
        swap(ei, et);

    // Compute transmitted ray direction
    float sini2 = SinTheta2(wo);
    float eta = ei / et;
    float sint2 = eta * eta * sini2;

    // Handle total internal reflection for transmission
    if (sint2 >= 1.0f) return Spectrum<brdf_d>();
    float cost = sqrtf(max(0.f, 1.f - sint2));
    if (entering) cost = -cost;
    float sintOverSini = eta;
    *wi = Normalize(Vector<local_s>(sintOverSini * -wo.x * meters, sintOverSini * -wo.y * meters, cost*meters));
    *pdf = mfloat<invsolidangle_d>(1.0f);
    Spectrum<rho_d> F = fresnel.Evaluate(CosTheta(wo));
    return (et*et)/(ei*ei) * (Spectrum<rho_d>(1.0f)-F) * T / AbsProjectTheta(*wi) / steradians;
}


Spectrum<brdf_d> Lambertian::f(const Direction<local_s> &wo, const Direction<local_s> &wi) const {
    return R / hemisphericalProjectedAngle;
}


Spectrum<brdf_d> OrenNayar::f(const Direction<local_s> &wo, const Direction<local_s> &wi) const {
    float sinthetai = SinTheta(wi);
    float sinthetao = SinTheta(wo);
    // Compute cosine term of Oren-Nayar model
    float maxcos = 0.0f;
    if (sinthetai > 1e-4 && sinthetao > 1e-4) {
        float sinphii = SinPhi(wi), cosphii = CosPhi(wi);
        float sinphio = SinPhi(wo), cosphio = CosPhi(wo);
        float dcos = cosphii * cosphio + sinphii * sinphio;
        maxcos = max(0.0f, dcos);
    }

    // Compute sine and tangent terms of Oren-Nayar model
    float sinalpha;
    float tanbeta;
    if (AbsCosTheta(wi) > AbsCosTheta(wo)) {
        sinalpha = sinthetao;
        tanbeta = sinthetai / AbsCosTheta(wi);
    }
    else {
        sinalpha = sinthetai;
        tanbeta = sinthetao / AbsCosTheta(wo);
    }
    return R / hemisphericalProjectedAngle * (A + B * maxcos * sinalpha * tanbeta);
}


Microfacet::Microfacet(const Spectrum<rho_d> &reflectance, Fresnel *f,
                       MicrofacetDistribution *d)
    : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
     R(reflectance), distribution(d), fresnel(f) {
}


Spectrum<brdf_d> Microfacet::f(const Direction<local_s> &wo, const Direction<local_s> &wi) const {
    float cosThetaO = AbsCosTheta(wo);
    mfloat<proj_d> cosThetaI = AbsProjectTheta(wi);
    if (cosThetaI == mfloat<proj_d>(0.0f) || cosThetaO == (0.0f)) return Spectrum<brdf_d>();
    Direction<local_s> wh = Bisector(wi, wo);
    float cosThetaH = Dot(wi, wh);
    Spectrum<rho_d> F = fresnel->Evaluate(cosThetaH);
    return R * distribution->D(wh) * G(wo, wi, wh) * F /
               (4.f * cosThetaI * cosThetaO);
}

Spectrum<brdf_d> BxDF::Sample_f(const Direction<local_s> &wo, Direction<local_s> *wi,
                        float u1, float u2, mfloat<invsolidangle_d> *pdf) const {
    // Cosine-sample the hemisphere, flipping the direction if necessary
    *wi = CosineSampleHemisphere<local_s>(u1, u2);
    if (wo.z < 0.) wi->z *= -1.f;
    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}


mfloat<invsolidangle_d> BxDF::Pdf(const Direction<local_s> &wo, const Direction<local_s> &wi) const {
    return SameHemisphere(wo, wi) ? AbsProjectTheta(wi) / (hemisphericalProjectedAngle) 
        : mfloat<invsolidangle_d>();
}


mfloat<invsolidangle_d> BRDFToBTDF::Pdf(const Direction<local_s> &wo,
        const Direction<local_s> &wi) const {
    return brdf->Pdf(wo, -wi);
}


Spectrum<brdf_d> Microfacet::Sample_f(const Direction<local_s> &wo, Direction<local_s> *wi,
                              float u1, float u2, mfloat<invsolidangle_d> *pdf) const {
    distribution->Sample_f(wo, wi, u1, u2, pdf);
    if (!SameHemisphere(wo, *wi)) return Spectrum<brdf_d>();
    return f(wo, *wi);
}


mfloat<invsolidangle_d> Microfacet::Pdf(const Direction<local_s> &wo, const Direction<local_s> &wi) const {
    if (!SameHemisphere(wo, wi)) return mfloat<invsolidangle_d>(0.f);
    return distribution->Pdf(wo, wi);
}


void Blinn::Sample_f(const Direction<local_s> &wo, Direction<local_s> *wi, float u1, float u2,
                     mfloat<invsolidangle_d> *pdf) const {
    // Compute sampled half-angle vector $\wh$ for Blinn distribution
    float costheta = powf(u1, 1.f / (exponent+1));
    float sintheta = sqrtf(max(0.f, 1.f - costheta*costheta));
    mfloat<angle_d> phi = u2 * circleAngle;
    Direction<local_s> wh = SphericalDirection<local_s>(sintheta, costheta, phi);
    if (!SameHemisphere(wo, wh)) wh = -wh;

    // Compute incident direction by reflecting about $\wh$
    *wi = Reflect(wo, wh);

    // Compute PDF for $\wi$ from Blinn distribution
    mfloat<invsolidangle_d> blinn_pdf = ((exponent + 1.f) * powf(costheta, exponent)) /
                      (hemisphericalAngle * 4.f * Dot(wo, wh));
    if (Dot(wo, wh) <= 0.f) blinn_pdf = mfloat<invsolidangle_d>(0.0f);
    *pdf = blinn_pdf;
}


mfloat<invsolidangle_d> Blinn::Pdf(const Direction<local_s> &wo, const Direction<local_s> &wi) const {
    Direction<local_s> wh = Bisector(wo, wi);
    float costheta = AbsCosTheta(wh);
    // Compute PDF for $\wi$ from Blinn distribution
    mfloat<invsolidangle_d> blinn_pdf = ((exponent + 1.f) * powf(costheta, exponent)) /
                      (hemisphericalAngle * 4.f * Dot(wo, wh));
    if (Dot(wo, wh) <= 0.f) blinn_pdf = mfloat<invsolidangle_d>(0.0f);
    return blinn_pdf;
}


void Anisotropic::Sample_f(const Direction<local_s> &wo, Direction<local_s> *wi,
                           float u1, float u2, mfloat<invsolidangle_d> *pdf) const {
    // Sample from first quadrant and remap to hemisphere to sample $\wh$
    mfloat<angle_d> phi;
    float costheta;
    if (u1 < .25f) {
        sampleFirstQuadrant(4.f * u1, u2, &phi, &costheta);
    } else if (u1 < .5f) {
        u1 = 4.f * (.5f - u1);
        sampleFirstQuadrant(u1, u2, &phi, &costheta);
        phi = halfCircleAngle - phi;
    } else if (u1 < .75f) {
        u1 = 4.f * (u1 - .5f);
        sampleFirstQuadrant(u1, u2, &phi, &costheta);
        phi += halfCircleAngle;
    } else {
        u1 = 4.f * (1.f - u1);
        sampleFirstQuadrant(u1, u2, &phi, &costheta);
        phi = circleAngle - phi;
    }
    float sintheta = sqrtf(max(0.f, 1.f - costheta*costheta));
    Direction<local_s> wh = SphericalDirection<local_s>(sintheta, costheta, phi);
    if (!SameHemisphere(wo, wh)) wh = -wh;

    // Compute incident direction by reflecting about $\wh$
    *wi = Reflect(wo, wh);

    // Compute PDF for $\wi$ from anisotropic distribution
    float costhetah = AbsCosTheta(wh);
    float ds = 1.f - costhetah * costhetah;
    mfloat<invsolidangle_d> anisotropic_pdf = mfloat<invsolidangle_d>(0.f);
    if (ds > 0.f && Dot(wo, wh) > 0.f) {
        float e = (ex * wh.x * wh.x + ey * wh.y * wh.y) / ds;
        anisotropic_pdf = sqrtf((ex+1.f) * (ey+1.f)) / hemisphericalAngle *
            powf(costhetah, e) / (4.f * Dot(wo, wh));
    }
    *pdf = anisotropic_pdf;
}


void Anisotropic::sampleFirstQuadrant(float u1, float u2,
        mfloat<angle_d> *phi, float *costheta) const {
    if (ex == ey)
        *phi = halfCircleAngle * u1;
    else
        *phi = _atan(sqrtf((ex+1.f) / (ey+1.f)) *
                     tanf(M_PI * u1 * 0.5f));
    float cosphi = cos(*phi), sinphi = sin(*phi);
    *costheta = powf(u2, 1.f/(ex * cosphi * cosphi +
                              ey * sinphi * sinphi + 1));
}


mfloat<invsolidangle_d> Anisotropic::Pdf(const Direction<local_s> &wo, const Direction<local_s> &wi) const {
    Direction<local_s> wh = Bisector(wo, wi);
    // Compute PDF for $\wi$ from anisotropic distribution
    float costhetah = AbsCosTheta(wh);
    float ds = 1.f - costhetah * costhetah;
    mfloat<invsolidangle_d> anisotropic_pdf = mfloat<invsolidangle_d>(0.f);
    if (ds > 0.f && Dot(wo, wh) > 0.f) {
        float e = (ex * wh.x * wh.x + ey * wh.y * wh.y) / ds;
        anisotropic_pdf = sqrtf((ex+1.f) * (ey+1.f)) / hemisphericalAngle *
                  powf(costhetah, e)  / (4.f * Dot(wo, wh));
    }
    return anisotropic_pdf;
}

Spectrum<rho_d> BxDF::rho(const Direction<local_s> &w, int nSamples,
                   const float *samples) const {
    Spectrum<rho_d> r;
    for (int i = 0; i < nSamples; ++i) {
        // Estimate one term of $\rho_\roman{hd}$
        Direction<local_s> wi;
        mfloat<invsolidangle_d> pdf(0.0f);
        Spectrum<brdf_d> f = Sample_f(w, &wi, samples[2*i], samples[2*i+1], &pdf);
        if (pdf > mfloat<invsolidangle_d>(0.0f)) 
            r += f * AbsProjectTheta(wi) / pdf;
    }
    return r / float(nSamples);
}


Spectrum<rho_d> BxDF::rho(int nSamples, const float *samples1,
                   const float *samples2) const {
    Spectrum<projsolidangle_d> r;
    for (int i = 0; i < nSamples; ++i) {
        // Estimate one term of $\rho_\roman{hh}$
        Direction<local_s> wo, wi;
        wo = UniformSampleHemisphere<local_s>(samples1[2*i], samples1[2*i+1]);
        mfloat<invsolidangle_d> pdf_o = mfloat<unit_d>(1.0f) / hemisphericalAngle;
        mfloat<invsolidangle_d> pdf_i = mfloat<invsolidangle_d>(0.0f);
        Spectrum<brdf_d> f = Sample_f(wo, &wi, samples2[2*i], samples2[2*i+1], &pdf_i);
        if (pdf_i > mfloat<invsolidangle_d>(0.0f))
            r += f * AbsProjectTheta(wi) * AbsProjectTheta(wo) / (pdf_o * pdf_i);
    }
    return r / (hemisphericalProjectedAngle * nSamples);
}



// BSDF Method Definitions
BSDFSampleOffsets::BSDFSampleOffsets(int count, Sample *sample) {
    nSamples = count;
    componentOffset = sample->Add1D(nSamples);
    dirOffset = sample->Add2D(nSamples);
}


BSDFSample::BSDFSample(const Sample *sample,
                       const BSDFSampleOffsets &offsets, uint32_t n) {
    Assert(n < sample->n2D[offsets.dirOffset]);
    Assert(n < sample->n1D[offsets.componentOffset]);
    uDir[0] = sample->twoD[offsets.dirOffset][2*n];
    uDir[1] = sample->twoD[offsets.dirOffset][2*n+1];
    uComponent = sample->oneD[offsets.componentOffset][n];
}


Spectrum<brdf_d> BSDF::Sample_f(const Direction<world_s> &woW, Direction<world_s> *wiW,
                        const BSDFSample &bsdfSample, mfloat<invsolidangle_d> *pdf,
                        BxDFType flags, BxDFType *sampledType) const {
    // Choose which _BxDF_ to sample
    int matchingComps = NumComponents(flags);
    if (matchingComps == 0) {
        *pdf = mfloat<invsolidangle_d>(0.0f);
        return Spectrum<brdf_d>();
    }
    int which = min(Floor2Int(bsdfSample.uComponent * matchingComps),
                    matchingComps-1);
    BxDF *bxdf = NULL;
    int count = which;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags) && count-- == 0) {
            bxdf = bxdfs[i];
            break;
        }
    Assert(bxdf);

    // Sample chosen _BxDF_
    Direction<local_s> wo = WorldToLocal(woW);
    Direction<local_s> wi;
    *pdf = mfloat<invsolidangle_d>(0.0f);
    Spectrum<brdf_d> f = bxdf->Sample_f(wo, &wi, bsdfSample.uDir[0],
                                bsdfSample.uDir[1], pdf);
    if (*pdf == mfloat<invsolidangle_d>(0.0f)) 
        return Spectrum<brdf_d>();
    if (sampledType) *sampledType = bxdf->type;
    *wiW = LocalToWorld(wi);

    // Compute overall PDF with all matching _BxDF_s
    if (!(bxdf->type & BSDF_SPECULAR) && matchingComps > 1)
        for (int i = 0; i < nBxDFs; ++i)
            if (bxdfs[i] != bxdf && bxdfs[i]->MatchesFlags(flags))
                *pdf += bxdfs[i]->Pdf(wo, wi);
    if (matchingComps > 1) *pdf /= static_cast<float>(matchingComps);

    // Compute value of BSDF for sampled direction
    if (!(bxdf->type & BSDF_SPECULAR)) {
        f = Spectrum<brdf_d>();
        if (Dot(*wiW, ng) * Dot(woW, ng) > 0) // ignore BTDFs
            flags = BxDFType(flags & ~BSDF_TRANSMISSION);
        else // ignore BRDFs
            flags = BxDFType(flags & ~BSDF_REFLECTION);
        for (int i = 0; i < nBxDFs; ++i)
            if (bxdfs[i]->MatchesFlags(flags))
                f += bxdfs[i]->f(wo, wi);
    }
    return f;
}


mfloat<invsolidangle_d> BSDF::Pdf(const Direction<world_s> &woW, const Direction<world_s> &wiW,
        BxDFType flags) const {
    if (nBxDFs == 0.) return mfloat<invsolidangle_d>(0.0f);
    Direction<local_s> wo = WorldToLocal(woW), wi = WorldToLocal(wiW);
    mfloat<invsolidangle_d> pdf(0.0f);
    int matchingComps = 0;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags)) {
            ++matchingComps;
            pdf += bxdfs[i]->Pdf(wo, wi);
        }
    return matchingComps > 0 ? pdf / matchingComps : mfloat<invsolidangle_d>(0.0f);
}


BSDF::BSDF(const DifferentialGeometry &dg, const Normal<world_s> &ngeom,
           float e)
    : dgShading(dg), eta(e) {
    ng = ngeom;
    nn = dgShading.nn;
    sn = Normalize(dgShading.dpdu);
    tn = Orthogonal(nn, sn);
    nBxDFs = 0;
}


Spectrum<brdf_d> BSDF::f(const Direction<world_s> &woW, const Direction<world_s> &wiW,
                 BxDFType flags) const {
    Direction<local_s> wi = WorldToLocal(wiW), wo = WorldToLocal(woW);
    if (Dot(wiW, ng) * Dot(woW, ng) > 0) // ignore BTDFs
        flags = BxDFType(flags & ~BSDF_TRANSMISSION);
    else // ignore BRDFs
        flags = BxDFType(flags & ~BSDF_REFLECTION);
    Spectrum<brdf_d> f;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags))
            f += bxdfs[i]->f(wo, wi);
    return f;
}


Spectrum<rho_d> BSDF::rho(RNG &rng, BxDFType flags, int sqrtSamples) const {
    int nSamples = sqrtSamples * sqrtSamples;
    float *s1 = ALLOCA(float, 2 * nSamples);
    StratifiedSample2D(s1, sqrtSamples, sqrtSamples, rng);
    float *s2 = ALLOCA(float, 2 * nSamples);
    StratifiedSample2D(s2, sqrtSamples, sqrtSamples, rng);

    Spectrum<rho_d> ret;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags))
            ret += bxdfs[i]->rho(nSamples, s1, s2);
    return ret;
}


Spectrum<rho_d> BSDF::rho(const Direction<local_s> &wo, RNG &rng, BxDFType flags,
                   int sqrtSamples) const {
    int nSamples = sqrtSamples * sqrtSamples;
    float *s1 = ALLOCA(float, 2 * nSamples);
    StratifiedSample2D(s1, sqrtSamples, sqrtSamples, rng);

    Spectrum<rho_d> ret;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags))
            ret += bxdfs[i]->rho(wo, nSamples, s1);
    return ret;
}


