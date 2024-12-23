
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


// core/shape.cpp*
#include "shape.h"

// Shape Method Definitions
Shape::~Shape() {
}


Shape::Shape(const Transform<object_s, world_s> *o2w, const Transform<world_s, object_s> *w2o, bool ro)
    : ObjectToWorld(o2w), WorldToObject(w2o), ReverseOrientation(ro),
      TransformSwapsHandedness(o2w->SwapsHandedness()),
      shapeId(nextshapeId++) {
    // Update shape creation statistics
    PBRT_CREATED_SHAPE(this);
}


uint32_t Shape::nextshapeId = 1;
BBox<world_s> Shape::WorldBound() const {
    return (*ObjectToWorld)(ObjectBound());
}


bool Shape::CanIntersect() const {
    return true;
}


void Shape::Refine(std::vector<Reference<Shape> > &refined) const {
    Severe("Unimplemented Shape::Refine() method called");
}


bool Shape::Intersect(const Ray<world_s> &ray, mfloat<length_d> *tHit, mfloat<length_d> *rayEpsilon,
                      DifferentialGeometry *dg) const {
    Severe("Unimplemented Shape::Intersect() method called");
    return false;
}


bool Shape::IntersectP(const Ray<world_s> &ray) const {
    Severe("Unimplemented Shape::IntersectP() method called");
    return false;
}


mfloat<area_d> Shape::Area() const {
    Severe("Unimplemented Shape::Area() method called");
    return 0.0f*sqr_meters;
}


mfloat<invsolidangle_d> Shape::Pdf(const Point<world_s> &p, const Direction<world_s> &wi) const {
    // Intersect sample ray with area light geometry
    DifferentialGeometry dgLight;
    Ray<world_s> ray(p, wi, mfloat<length_d>(1e-3f));
    mfloat<length_d> thit, rayEpsilon;
    if (!Intersect(ray, &thit, &rayEpsilon, &dgLight)) return mfloat<invsolidangle_d>(0.0f);

    // Convert light sample weight to solid angle measure
    mfloat<invsolidangle_d> pdf = mfloat<invsolidangle_d>(DistanceSquared(p, ray(thit)) /
                (abs(Dot(dgLight.nn, -wi)) * Area()));
    if (isinf(pdf.v)) pdf = mfloat<invsolidangle_d>(0.0f);
    return pdf;
}


