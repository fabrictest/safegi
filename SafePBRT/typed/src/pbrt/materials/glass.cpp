
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


// materials/glass.cpp*
#include "materials/glass.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"

// GlassMaterial Method Definitions
BSDF *GlassMaterial::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    float ior = index->Evaluate(dgs);
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn, ior);
    Spectrum<rho_d> R = Kr->Evaluate(dgs).Clamp();
    Spectrum<rho_d> T = Kt->Evaluate(dgs).Clamp();
    if (!R.IsBlack())
        bsdf->Add(BSDF_ALLOC(arena, SpecularReflection)(R,
            BSDF_ALLOC(arena, FresnelDielectric)(1., ior)));
    if (!T.IsBlack())
        bsdf->Add(BSDF_ALLOC(arena, SpecularTransmission)(T, 1., ior));
    return bsdf;
}


GlassMaterial *CreateGlassMaterial(const Transform<material_s, world_s> &xform,
        const TextureParams &mp) {
    Reference<Texture<Spectrum<rho_d> > > Kr = mp.GetSpectrumTexture("Kr", Spectrum<rho_d>(1.f));
    Reference<Texture<Spectrum<rho_d> > > Kt = mp.GetSpectrumTexture("Kt", Spectrum<rho_d>(1.f));
    Reference<Texture<float> > index = mp.GetFloatTexture("index", 1.5f);
    Reference<Texture<mfloat<length_d> > > bumpMap = mp.GetFloatTexture("bumpmap", mfloat<length_d>(0.0f));
    return new GlassMaterial(Kr, Kt, index, bumpMap);
}


