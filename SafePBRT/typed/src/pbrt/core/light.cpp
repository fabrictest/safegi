
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


// core/light.cpp*
#include "light.h"
#include "scene.h"
#include "montecarlo.h"
#include "paramset.h"
#include "sampler.h"

// Light Method Definitions
Light::~Light() {
}


Spectrum<radiance_d> Light::Le(const RayDifferential<world_s> &) const {
    return Spectrum<radiance_d>();
}


LightSampleOffsets::LightSampleOffsets(int count, Sample *sample) {
    nSamples = count;
    componentOffset = sample->Add1D(nSamples);
    posOffset = sample->Add2D(nSamples);
}


LightSample::LightSample(const Sample *sample,
        const LightSampleOffsets &offsets, uint32_t n) {
    Assert(n < sample->n2D[offsets.posOffset]);
    Assert(n < sample->n1D[offsets.componentOffset]);
    uPos[0] = sample->twoD[offsets.posOffset][2*n];
    uPos[1] = sample->twoD[offsets.posOffset][2*n+1];
    uComponent = sample->oneD[offsets.componentOffset][n];
}

// ShapeSet Method Definitions
ShapeSet::ShapeSet(const Reference<Shape> &s) {
    std::vector<Reference<Shape> > todo;
    todo.push_back(s);
    while (todo.size()) {
        Reference<Shape> sh = todo.back();
        todo.pop_back();
        if (sh->CanIntersect())
            shapes.push_back(sh);
        else
            sh->Refine(todo);
    }
    if (shapes.size() > 64)
        Warning("Area light geometry turned into %d shapes; "
            "may be very inefficient.", (int)shapes.size());

    // Compute total area of shapes in _ShapeSet_ and area CDF
    sumArea = mfloat<area_d>(0.f);
    for (uint32_t i = 0; i < shapes.size(); ++i) {
        mfloat<area_d> a = shapes[i]->Area();
        areas.push_back(a);
        sumArea += a;
    }
    areaDistribution = new Distribution1D<area_d>(&areas[0], areas.size());
}


ShapeSet::~ShapeSet() {
    delete areaDistribution;
}


Point<world_s> ShapeSet::Sample(const Point<world_s> &p, const LightSample &ls,
                       Normal<world_s> *Ns) const {
    int sn = areaDistribution->SampleDiscrete(ls.uComponent, NULL);
    return shapes[sn]->Sample(p, ls.uPos[0], ls.uPos[1], Ns);
}


Point<world_s> ShapeSet::Sample(const LightSample &ls, Normal<world_s> *Ns) const {
    int sn = areaDistribution->SampleDiscrete(ls.uComponent, NULL);
    return shapes[sn]->Sample(ls.uPos[0], ls.uPos[1], Ns);
}


mfloat<invsolidangle_d> ShapeSet::Pdf(const Point<world_s> &p, const Direction<world_s> &wi) const {
    mfloat<invsolidangle_d> pdf(0.0f);
    for (uint32_t i = 0; i < shapes.size(); ++i)
        pdf += shapes[i]->Pdf(p, wi) * (areas[i] / sumArea);
    return pdf;
}


mfloat<invarea_d> ShapeSet::Pdf(const Point<world_s> &p) const {
    mfloat<unit_d> pdf(0.0f);
    for (uint32_t i = 0; i < shapes.size(); ++i)
        pdf += areas[i] * shapes[i]->Pdf(p);
    return pdf / sumArea;
}


