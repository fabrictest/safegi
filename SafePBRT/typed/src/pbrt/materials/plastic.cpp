
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


// materials/plastic.cpp*
#include "materials/plastic.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"

// PlasticMaterial Method Definitions
BSDF *PlasticMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
                               const DifferentialGeometry &dgShading,
                               MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);
    Spectrum<rho_d> kd = Kd->Evaluate(dgs).Clamp();
    BxDF *diff = BSDF_ALLOC(arena, Lambertian)(kd);
    Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(1.5f, 1.f);
    Spectrum<rho_d> ks = Ks->Evaluate(dgs).Clamp();
    float rough = roughness->Evaluate(dgs);
    BxDF *spec = BSDF_ALLOC(arena, Microfacet)
                   (ks, fresnel, BSDF_ALLOC(arena, Blinn)(1.f / rough));
    bsdf->Add(diff);
    bsdf->Add(spec);
    return bsdf;
}


PlasticMaterial *CreatePlasticMaterial(const Transform<material_s, world_s> &xform,
        const TextureParams &mp) {
    Reference<Texture<Spectrum<rho_d> > > Kd = mp.GetSpectrumTexture("Kd", Spectrum<rho_d>(1.f));
    Reference<Texture<Spectrum<rho_d> > > Ks = mp.GetSpectrumTexture("Ks", Spectrum<rho_d>(1.f));
    Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
    Reference<Texture<mfloat<length_d> > > bumpMap = mp.GetFloatTexture<length_d>("bumpmap", mfloat<length_d>(0.0f));
    return new PlasticMaterial(Kd, Ks, roughness, bumpMap);
}


