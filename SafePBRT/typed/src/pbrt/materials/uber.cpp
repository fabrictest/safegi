
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


// materials/uber.cpp*
#include "materials/uber.h"
#include "spectrum.h"
#include "reflection.h"
#include "texture.h"
#include "paramset.h"

// UberMaterial Method Definitions
BSDF *UberMaterial::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);

    Spectrum<rho_d> op = opacity->Evaluate(dgs).Clamp();
    if (op != Spectrum<rho_d>(1.0f)) {
        BxDF *tr = BSDF_ALLOC(arena, SpecularTransmission)(-op + Spectrum<rho_d>(1.0f), 1., 1.);
        bsdf->Add(tr);
    }

    Spectrum<rho_d> kd = op * Kd->Evaluate(dgs).Clamp();
    if (!kd.IsBlack()) {
        BxDF *diff = BSDF_ALLOC(arena, Lambertian)(kd);
        bsdf->Add(diff);
    }

    Spectrum<rho_d> ks = op * Ks->Evaluate(dgs).Clamp();
    if (!ks.IsBlack()) {
        Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(1.5f, 1.f);
        float rough = roughness->Evaluate(dgs);
        BxDF *spec = BSDF_ALLOC(arena, Microfacet)(ks, fresnel, BSDF_ALLOC(arena, Blinn)(1.f / rough));
        bsdf->Add(spec);
    }

    Spectrum<rho_d> kr = op * Kr->Evaluate(dgs).Clamp();
    if (!kr.IsBlack()) {
        Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(1.5f, 1.f);
        bsdf->Add(BSDF_ALLOC(arena, SpecularReflection)(kr, fresnel));
    }

    return bsdf;
}


UberMaterial *CreateUberMaterial(const Transform<material_s, world_s> &xform,
        const TextureParams &mp) {
    Reference<Texture<Spectrum<rho_d> > > Kd = mp.GetSpectrumTexture("Kd", Spectrum<rho_d>(1.0f));
    Reference<Texture<Spectrum<rho_d> > > Ks = mp.GetSpectrumTexture("Ks", Spectrum<rho_d>(1.0f));
    Reference<Texture<Spectrum<rho_d> > > Kr = mp.GetSpectrumTexture("Kr", Spectrum<rho_d>(0.0f));
    Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
    Reference<Texture<Spectrum<rho_d> > > opacity = mp.GetSpectrumTexture("opacity", Spectrum<rho_d>(1.0f));
    Reference<Texture<mfloat<length_d> > > bumpMap = mp.GetFloatTexture("bumpmap", mfloat<length_d>(0.0f));
    return new UberMaterial(Kd, Ks, Kr, roughness, opacity, bumpMap);
}


