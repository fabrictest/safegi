
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


// textures/constant.cpp*
#include "textures/constant.h"
#include "paramset.h"

ConstantTexture<mfloat<unit_d> > *CreateConstantFloatTexture(const Transform<texture_s, world_s> &tex2world,
                                                             const TextureParams &tp)
{
    return new ConstantTexture<mfloat<unit_d> >(tp.FindFloat("value", mfloat<unit_d>(1.0f)));
}


ConstantTexture<Spectrum<rho_d> > *CreateConstantSpectrumTexture(const Transform<texture_s, world_s> &tex2world,
                                                                 const TextureParams &tp) 
{
    return new ConstantTexture<Spectrum<rho_d> >(tp.FindSpectrum("value", Spectrum<rho_d>(mfloat<rho_d>(1.0))));
}
// ConstantTexture Method Definitions


