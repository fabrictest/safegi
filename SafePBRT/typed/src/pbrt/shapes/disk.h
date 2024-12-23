
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

#ifndef PBRT_SHAPES_DISK_H
#define PBRT_SHAPES_DISK_H

// shapes/disk.h*
#include "shape.h"

// Disk Declarations
class Disk : public Shape {
public:
    // Disk Public Methods
    Disk(const Transform<object_s, world_s> *o2w, const Transform<world_s, object_s> *w2o, 
        bool ro, const mfloat<length_d> &height, const mfloat<length_d> &radius, 
        const mfloat<length_d> &innerRadius, float phiMax);
    BBox<object_s> ObjectBound() const;
    bool Intersect(const Ray<world_s> &ray, mfloat<length_d> *tHit, 
        mfloat<length_d> *rayEpsilon, DifferentialGeometry *dg) const;
    bool IntersectP(const Ray<world_s> &ray) const;
    mfloat<area_d> Area() const;
    Point<world_s> Sample(float u1, float u2, Normal<world_s> *Ns) const;
private:
    // Disk Private Data
    mfloat<length_d> height, radius, innerRadius;
    mfloat<angle_d> phiMax;
};


Disk *CreateDiskShape(const Transform<object_s, world_s> *o2w, 
    const Transform<world_s, object_s> *w2o, bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_DISK_H
