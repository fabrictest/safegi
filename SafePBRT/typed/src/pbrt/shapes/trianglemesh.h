
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

#ifndef PBRT_SHAPES_TRIANGLEMESH_H
#define PBRT_SHAPES_TRIANGLEMESH_H

// shapes/trianglemesh.h*
#include "shape.h"

// TriangleMesh Declarations
class TriangleMesh : public Shape {
public:
    // TriangleMesh Public Methods
    TriangleMesh(const Transform<object_s, world_s> *o2w, const Transform<world_s, object_s> *w2o, bool ro,
        int ntris, int nverts, const int *vptr,
        const Point<object_s> *P, const Normal<object_s> *N,
        const Vector<object_s> *S, const float *uv);
    ~TriangleMesh();
    BBox<object_s> ObjectBound() const;
    BBox<world_s> WorldBound() const;
    bool CanIntersect() const { return false; }
    void Refine(std::vector<Reference<Shape> > &refined) const;
    friend class Triangle;
    template <typename T> friend class VertexTexture;
protected:
    // TriangleMesh Proteted Data
    int ntris, nverts;
    int *vertexIndex;
    Point<world_s> *p;
    Normal<object_s> *n;
    Direction<object_s> *s;
    float *uvs;
};


class Triangle : public Shape {
public:
    // Triangle Public Methods
    Triangle(const Transform<object_s, world_s> *o2w, const Transform<world_s, object_s> *w2o, bool ro,
        TriangleMesh *m, int n)
        : Shape(o2w, w2o, ro) {
            mesh = m;
            v = &mesh->vertexIndex[3*n];
            PBRT_CREATED_TRIANGLE(this);
    }
    BBox<object_s> ObjectBound() const;
    BBox<world_s> WorldBound() const;
    bool Intersect(const Ray<world_s> &ray, mfloat<length_d> *tHit, mfloat<length_d> *rayEpsilon,
        DifferentialGeometry *dg) const;
    bool IntersectP(const Ray<world_s> &ray) const;
    void GetUVs(float uv[3][2]) const {
        if (mesh->uvs) {
            uv[0][0] = mesh->uvs[2*v[0]];
            uv[0][1] = mesh->uvs[2*v[0]+1];
            uv[1][0] = mesh->uvs[2*v[1]];
            uv[1][1] = mesh->uvs[2*v[1]+1];
            uv[2][0] = mesh->uvs[2*v[2]];
            uv[2][1] = mesh->uvs[2*v[2]+1];
        }
        else {
            uv[0][0] = 0.; uv[0][1] = 0.;
            uv[1][0] = 1.; uv[1][1] = 0.;
            uv[2][0] = 1.; uv[2][1] = 1.;
        }
    }
    mfloat<area_d> Area() const;
    virtual void GetShadingGeometry(const Transform<object_s, world_s> &obj2world,
        const DifferentialGeometry &dg,
        DifferentialGeometry *dgShading) const;
    Point<world_s> Sample(float u1, float u2, Normal<world_s> *Ns) const;
private:
    // Triangle Private Data
    Reference<TriangleMesh> mesh;
    int *v;
};


TriangleMesh *CreateTriangleMeshShape(const Transform<object_s, world_s> *o2w, const Transform<world_s, object_s> *w2o,
                                      bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_TRIANGLEMESH_H
