
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


// core/material.cpp*
#include "material.h"
#include "primitive.h"
#include "texture.h"
#include "spectrum.h"
#include "reflection.h"

// Material Method Definitions
Material::~Material() {
}


void Material::Bump(const Reference<Texture<mfloat<length_d> > > &d,
                    const DifferentialGeometry &dgGeom,
                    const DifferentialGeometry &dgs,
                    DifferentialGeometry *dgBump) {
    // Compute offset positions and evaluate displacement texture
    DifferentialGeometry dgEval = dgs;

    // Shift _dgEval_ _du_ in the $u$ direction
    float du = .5f * (fabsf(dgs.dudx) + fabsf(dgs.dudy));
    if (du == 0.f) du = .01f;
    dgEval.p = dgs.p + du * dgs.dpdu;
    dgEval.u = dgs.u + du;

    dgEval.nn = 
        Normalize<world_s>(Cross(dgs.dpdu, dgs.dpdv)/sqr_meters + du * dgs.dndu);
    mfloat<length_d> uDisplace = d->Evaluate(dgEval);
    if (isnan(__asfloat(uDisplace)))
        uDisplace = d->Evaluate(dgEval);

    // Shift _dgEval_ _dv_ in the $v$ direction
    float dv = .5f * (fabsf(dgs.dvdx) + fabsf(dgs.dvdy));
    if (dv == 0.f) dv = .01f;
    dgEval.p = dgs.p + mfloat<length_d>(dv) * dgs.dndv;
    dgEval.u = dgs.u;
    dgEval.v = dgs.v + dv;
    dgEval.nn = 
        Normalize<world_s>(Cross(dgs.dpdu, dgs.dpdv)/sqr_meters + dv * dgs.dndu);
    mfloat<length_d> vDisplace = d->Evaluate(dgEval);
    mfloat<length_d> displace = d->Evaluate(dgs);

    // Compute bump-mapped differential geometry
    *dgBump = dgs;
    dgBump->dpdu = dgs.dpdu + (uDisplace - displace) / du * (dgs.nn) +
                   displace * dgs.dndu;
    dgBump->dpdv = dgs.dpdv + (vDisplace - displace) / dv * (dgs.nn) +
                   displace * dgs.dndv;
    dgBump->nn = Normalize<world_s>(Cross(dgBump->dpdu, dgBump->dpdv));
    if (dgs.shape->ReverseOrientation ^ dgs.shape->TransformSwapsHandedness)
        dgBump->nn = -dgBump->nn;

    // Orient shading normal to match geometric normal
    dgBump->nn = Faceforward(dgBump->nn, dgGeom.nn);
}


