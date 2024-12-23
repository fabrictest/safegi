
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

#ifndef PBRT_CORE_SHAPE_H
#define PBRT_CORE_SHAPE_H

// core/shape.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
#include "diffgeom.h"
#include "memory.h"

// Shape Declarations
class Shape : public ReferenceCounted {
public:
    // Shape Interface
    Shape(const Transform<object_s, world_s> *o2w, const Transform<world_s, object_s> *w2o, bool ro);
    virtual ~Shape();
    virtual BBox<object_s> ObjectBound() const = 0;
    virtual BBox<world_s> WorldBound() const;
    virtual bool CanIntersect() const;
    virtual void Refine(std::vector<Reference<Shape> > &refined) const;
    virtual bool Intersect(const Ray<world_s> &ray, mfloat<length_d> *tHit,
                           mfloat<length_d> *rayEpsilon, DifferentialGeometry *dg) const;
    virtual bool IntersectP(const Ray<world_s> &ray) const;
    virtual void GetShadingGeometry(const Transform<object_s, world_s> &obj2world,
            const DifferentialGeometry &dg,
            DifferentialGeometry *dgShading) const {
        *dgShading = dg;
    }
    virtual mfloat<area_d> Area() const;
    virtual Point<world_s> Sample(float u1, float u2, Normal<world_s> *Ns) const {
        Severe("Unimplemented Shape::Sample() method called");
        return Point<world_s>();
    }
    virtual mfloat<invarea_d> Pdf(const Point<world_s> &Pshape) const {
        return mfloat<rho_d>(1.0f) / Area();
    }
    virtual Point<world_s> Sample(const Point<world_s> &P, float u1, float u2,
                         Normal<world_s> *Ns) const {
        return Sample(u1, u2, Ns);
    }
    virtual mfloat<invsolidangle_d> Pdf(const Point<world_s> &p, const Direction<world_s> &wi) const;

    // Shape Public Data
    const Transform<object_s, world_s> *ObjectToWorld;
    const Transform<world_s, object_s> *WorldToObject;
    const bool ReverseOrientation, TransformSwapsHandedness;
    const uint32_t shapeId;
    static uint32_t nextshapeId;
};



#endif // PBRT_CORE_SHAPE_H
