
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

#ifndef PBRT_MATERIALS_UBER_H
#define PBRT_MATERIALS_UBER_H

// materials/uber.h*
#include "pbrt.h"
#include "material.h"

// UberMaterial Declarations
class UberMaterial : public Material {
public:
    UberMaterial(Reference<Texture<Spectrum<rho_d> > > kd,
        Reference<Texture<Spectrum<rho_d> > > ks,
        Reference<Texture<Spectrum<rho_d> > > kr,
        Reference<Texture<float> > rough,
        Reference<Texture<Spectrum<rho_d> > > op,
        Reference<Texture<mfloat<length_d> > > bump) {
        Kd = kd;
        Ks = ks;
        Kr = kr;
        roughness = rough;
        opacity = op;
        bumpMap = bump;
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const;
private:
    // UberMaterial Private Data
    Reference<Texture<Spectrum<rho_d> > > Kd, Ks, Kr, opacity;
    Reference<Texture<float> > roughness;
    Reference<Texture<mfloat<length_d> > > bumpMap;
};


UberMaterial *CreateUberMaterial(const Transform<material_s, world_s> &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_UBER_H
