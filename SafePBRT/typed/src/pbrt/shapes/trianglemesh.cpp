
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


// shapes/trianglemesh.cpp*
#include "shapes/trianglemesh.h"
#include "paramset.h"
#include "montecarlo.h"

// TriangleMesh Method Definitions
TriangleMesh::TriangleMesh(const Transform<object_s, world_s> *o2w, const Transform<world_s, object_s> *w2o,
                           bool ro, int nt, int nv, const int *vi, const Point<object_s> *P,
                           const Normal<object_s> *N, const Vector<object_s> *S, const float *uv)
                           : Shape(o2w, w2o, ro) 
{
    ntris = nt;
    nverts = nv;
    vertexIndex = new int[3 * ntris];
    memcpy(vertexIndex, vi, 3 * ntris * sizeof(int));
    // Copy _uv_, _N_, and _S_ vertex data, if present
    if (uv) {
        uvs = new float[2*nverts];
        memcpy(uvs, uv, 2*nverts*sizeof(float));
    }
    else uvs = NULL;
    p = new Point<world_s>[nverts];
    if (N) {
        n = new Normal<object_s>[nverts];
        memcpy(n, N, nverts*sizeof(Normal<object_s>));
    }
    else n = NULL;
    if (S) {
        s = new Direction<object_s>[nverts];
        for (int i  = 0; i < nverts; ++i)
            s[i] = Normalize(S[i]);
        //memcpy(s, S, nverts*sizeof(Direction<object_s>));
    }
    else s = NULL;

    // Transform mesh vertices to world space
    for (int i  = 0; i < nverts; ++i)
        p[i] = (*ObjectToWorld)(P[i]);
}


TriangleMesh::~TriangleMesh() 
{
    delete[] vertexIndex;
    delete[] p;
    delete[] s;
    delete[] n;
    delete[] uvs;
}


BBox<object_s> TriangleMesh::ObjectBound() const 
{
    BBox<object_s> bobj;
    for (int i = 0; i < nverts; i++)
        bobj = Union(bobj, (*WorldToObject)(p[i]));
    return bobj;
}


BBox<world_s> TriangleMesh::WorldBound() const 
{
    BBox<world_s> worldBounds;
    for (int i = 0; i < nverts; i++)
        worldBounds = Union(worldBounds, p[i]);
    return worldBounds;
}


void TriangleMesh::Refine(std::vector<Reference<Shape> > &refined) const 
{
    for (int i = 0; i < ntris; ++i)
        refined.push_back(new Triangle(ObjectToWorld,
        WorldToObject, ReverseOrientation,
        (TriangleMesh *)this, i));
}


BBox<object_s> Triangle::ObjectBound() const 
{
    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point<world_s> &p1 = mesh->p[v[0]];
    const Point<world_s> &p2 = mesh->p[v[1]];
    const Point<world_s> &p3 = mesh->p[v[2]];
    return Union(BBox<object_s>((*WorldToObject)(p1), (*WorldToObject)(p2)),
        (*WorldToObject)(p3));
}


BBox<world_s> Triangle::WorldBound() const 
{
    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point<world_s> &p1 = mesh->p[v[0]];
    const Point<world_s> &p2 = mesh->p[v[1]];
    const Point<world_s> &p3 = mesh->p[v[2]];
    return Union(BBox<world_s>(p1, p2), p3);
}


bool Triangle::Intersect(const Ray<world_s> &ray, mfloat<length_d> *tHit, mfloat<length_d> *rayEpsilon,
                         DifferentialGeometry *dg) const 
{
    PBRT_RAY_TRIANGLE_INTERSECTION_TEST(const_cast<Ray *>(&ray), const_cast<Triangle *>(this));
    // Compute $\VEC{s}_1$

    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point<world_s> &p1 = mesh->p[v[0]];
    const Point<world_s> &p2 = mesh->p[v[1]];
    const Point<world_s> &p3 = mesh->p[v[2]];
    Vector<world_s> e1 = p2 - p1;
    Vector<world_s> e2 = p3 - p1;
    Vector<world_s> s1 = Cross(ray.d, e2);
    auto divisor = Dot(s1, e1);
    if (divisor == mfloat<area_d>(0.0))
        return false;
    mfloat<invarea_d> invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector<world_s> d = ray.o - p1;
    float b1 = Dot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    auto s2 = Cross(d, e1);
    auto b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    auto t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;

    // Fill in _DifferentialGeometry_ from triangle hit

    // Compute triangle partial derivatives
    Vector<world_s> dpdu, dpdv;
    float uvs[3][2];
    GetUVs(uvs);

    // Compute deltas for triangle partial derivatives
    float du1 = uvs[0][0] - uvs[2][0];
    float du2 = uvs[1][0] - uvs[2][0];
    float dv1 = uvs[0][1] - uvs[2][1];
    float dv2 = uvs[1][1] - uvs[2][1];
    auto dp1 = p1 - p3, dp2 = p2 - p3;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for triangle partial derivative matrix
        Direction<world_s> dddu, dddv;
        CoordinateSystem(Normalize<world_s>(Cross(e2, e1)), &dddu, &dddv);
        dpdu = dddu*meters;
        dpdv = dddv*meters;
    }
    else {
        auto invdet = 1.f / determinant;
        dpdu = (dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }

    // Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - b1 - b2;
    float tu = b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0];
    float tv = b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1];
    *dg = DifferentialGeometry(ray(t), dpdu, dpdv,
        Normal<world_s>(), Normal<world_s>(),
        tu, tv, this);
    *tHit = t;
    *rayEpsilon = 1e-4f * *tHit;
    PBRT_RAY_TRIANGLE_INTERSECTION_HIT(const_cast<Ray *>(&ray), t);
    return true;
}


bool Triangle::IntersectP(const Ray<world_s> &ray) const 
{
    PBRT_RAY_TRIANGLE_INTERSECTIONP_TEST(const_cast<Ray *>(&ray), const_cast<Triangle *>(this));
    // Compute $\VEC{s}_1$

    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point<world_s> &p1 = mesh->p[v[0]];
    const Point<world_s> &p2 = mesh->p[v[1]];
    const Point<world_s> &p3 = mesh->p[v[2]];
    Vector<world_s> e1 = p2 - p1;
    Vector<world_s> e2 = p3 - p1;
    Vector<world_s> s1 = Cross(ray.d, e2);
    auto divisor = Dot(s1, e1);
    if (divisor == 0.0 * sqr_meters)
        return false;
    auto invDivisor = mfloat<unit_d>(1.0f) / divisor;

    // Compute first barycentric coordinate
    Vector<world_s> d = ray.o - p1;
    float b1 = Dot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    auto s2 = Cross(Normalize(d), e1);
    float b2 = Dot(ray.d, s2) * d.Length() * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    auto t = Dot(e2, s2) * d.Length() * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;
    PBRT_RAY_TRIANGLE_INTERSECTIONP_HIT(const_cast<Ray *>(&ray), t);
    return true;
}


mfloat<area_d> Triangle::Area() const {
    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point<world_s> &p1 = mesh->p[v[0]];
    const Point<world_s> &p2 = mesh->p[v[1]];
    const Point<world_s> &p3 = mesh->p[v[2]];
    const Vector<world_s> e0 = p2 - p1;
    const Vector<world_s> e1 = p3 - p1;
    return 0.5f * Cross(Normalize(e0), e1).Length() * e0.Length();
}


void Triangle::GetShadingGeometry(const Transform<object_s, world_s> &obj2world,
                                  const DifferentialGeometry &dg,
                                  DifferentialGeometry *dgShading) const 
{
    if (!mesh->n && !mesh->s) {
        *dgShading = dg;
        return;
    }
    // Initialize _Triangle_ shading geometry with _n_ and _s_

    // Compute barycentric coordinates for point
    float b[3];

    // Initialize _A_ and _C_ matrices for barycentrics
    float uv[3][2];
    GetUVs(uv);
    float A[2][2] =
    { { uv[1][0] - uv[0][0], uv[2][0] - uv[0][0] },
    { uv[1][1] - uv[0][1], uv[2][1] - uv[0][1] } };
    float C[2] = { dg.u - uv[0][0], dg.v - uv[0][1] };
    if (!SolveLinearSystem2x2(A, C, &b[1], &b[2])) {
        // Handle degenerate parametric mapping
        b[0] = b[1] = b[2] = 1.f/3.f;
    }
    else
        b[0] = 1.f - b[1] - b[2];

    // Use _n_ and _s_ to compute shading tangents for triangle, _ss_ and _ts_
    Normal<world_s> ns;
    Direction<world_s> ss;
    Direction<world_s> ts;
    if (mesh->n) ns = obj2world(Normalize<object_s>(b[0] * (mesh->n[v[0]]) +
        b[1] * (mesh->n[v[1]]) +
        b[2] * (mesh->n[v[2]])));
    else   ns = dg.nn;
    if (mesh->s) ss = obj2world(Normalize<object_s>(b[0] * (mesh->s[v[0]]) +
        b[1] * (mesh->s[v[1]]) +
        b[2] * (mesh->s[v[2]])));
    else   ss = Normalize(dg.dpdu);
    ts = Orthogonal(ss, ns);
    ss = Orthogonal(ts, ns);
    Normal<object_s> dndu, dndv;

    // Compute $\dndu$ and $\dndv$ for triangle shading geometry
    if (mesh->n) {
        float uvs[3][2];
        GetUVs(uvs);
        // Compute deltas for triangle partial derivatives of normal
        float du1 = uvs[0][0] - uvs[2][0];
        float du2 = uvs[1][0] - uvs[2][0];
        float dv1 = uvs[0][1] - uvs[2][1];
        float dv2 = uvs[1][1] - uvs[2][1];
        auto dn1 = (mesh->n[v[0]]) - (mesh->n[v[2]]);
        auto dn2 = (mesh->n[v[1]]) - (mesh->n[v[2]]);
        auto determinant = du1 * dv2 - dv1 * du2;
        if (determinant == 0.f)
            dndu = dndv = Normal<object_s>();
        else {
            auto invdet = 1.f / determinant;
            dndu = AsNormalDiff<object_s>((dv2 * dn1 - dv1 * dn2) * invdet);
            dndv = AsNormalDiff<object_s>((-du2 * dn1 + du1 * dn2) * invdet);
        }
    }
    else
        dndu = dndv = Normal<object_s>();
    *dgShading = DifferentialGeometry(dg.p, ss*meters, ts*meters,
        (*ObjectToWorld)(dndu), (*ObjectToWorld)(dndv),
        dg.u, dg.v, dg.shape);
    dgShading->dudx = dg.dudx;  dgShading->dvdx = dg.dvdx;
    dgShading->dudy = dg.dudy;  dgShading->dvdy = dg.dvdy;
    dgShading->dpdx = dg.dpdx;  dgShading->dpdy = dg.dpdy;
}


TriangleMesh *CreateTriangleMeshShape(const Transform<object_s, world_s> *o2w, const Transform<world_s, object_s> *w2o,
                                      bool reverseOrientation, const ParamSet &params) 
{
    int nvi, npi, nuvi, nsi, nni;
    const int *vi = params.FindInt("indices", &nvi);
    const Point<object_s> *P = params.FindPoint<object_s>("P", &npi);
    const float *uvs = params.FindFloat("uv", &nuvi);
    if (!uvs) uvs = params.FindFloat("st", &nuvi);
    bool discardDegnerateUVs = params.FindOneBool("discarddegenerateUVs", false);
    // XXX should complain if uvs aren't an array of 2...
    if (uvs) {
        if (nuvi < 2 * npi) {
            Error("Not enough of \"uv\"s for triangle mesh.  Expencted %d, "
                "found %d.  Discarding.\n", 2*npi, nuvi);
            uvs = NULL;
        }
        else if (nuvi > 2 * npi)
            Warning("More \"uv\"s provided than will be used for triangle "
            "mesh.  (%d expcted, %d found)\n", 2*npi, nuvi);
    }
    if (!vi || !P) return NULL;
    const Vector<object_s> *S = params.FindVector<object_s>("S", &nsi);
    if (S && nsi != npi) {
        Error("Number of \"S\"s for triangle mesh must match \"P\"s");
        S = NULL;
    }
    const Normal<object_s> *N = params.FindNormal<object_s>("N", &nni);
    if (N && nni != npi) {
        Error("Number of \"N\"s for triangle mesh must match \"P\"s");
        N = NULL;
    }
    if (discardDegnerateUVs && uvs && N) {
        // if there are normals, check for bad uv's that
        // give degenerate mappings; discard them if so
        const int *vp = vi;
        for (int i = 0; i < nvi; i += 3, vp += 3) {
            Vector<object_s> e = P[vp[0]]-P[vp[1]];
            mfloat<area_d> area = .5f * Cross(Normalize(e), P[vp[2]]-P[vp[1]]).Length() * e.Length();
            if (area < 1e-7 * sqr_meters) continue; // ignore degenerate tris.
            if ((uvs[2*vp[0]] == uvs[2*vp[1]] &&
                uvs[2*vp[0]+1] == uvs[2*vp[1]+1]) ||
                (uvs[2*vp[1]] == uvs[2*vp[2]] &&
                uvs[2*vp[1]+1] == uvs[2*vp[2]+1]) ||
                (uvs[2*vp[2]] == uvs[2*vp[0]] &&
                uvs[2*vp[2]+1] == uvs[2*vp[0]+1])) {
                    Warning("Degenerate uv coordinates in triangle mesh.  Discarding all uvs.");
                    uvs = NULL;
                    break;
            }
        }
    }
    for (int i = 0; i < nvi; ++i)
        if (vi[i] >= npi) {
            Error("trianglemesh has out of-bounds vertex index %d (%d \"P\" values were given",
                vi[i], npi);
            return NULL;
        }
        return new TriangleMesh(o2w, w2o, reverseOrientation, nvi/3, npi, vi, P,
            N, S, uvs);
}


Point<world_s> Triangle::Sample(float u1, float u2, Normal<world_s> *Ns) const {
    float b1, b2;
    UniformSampleTriangle(u1, u2, &b1, &b2);
    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point<world_s> &p1 = mesh->p[v[0]];
    const Point<world_s> &p2 = mesh->p[v[1]];
    const Point<world_s> &p3 = mesh->p[v[2]];
    Point<world_s> p = Point<world_s>() + (b1 * AsVector(p1) + b2 * AsVector(p2) + (1.f - b1 - b2) * AsVector(p3));
    *Ns = Orthogonal(Normalize(p2-p1), Normalize(p3-p1));
    if (ReverseOrientation) *Ns = -(*Ns);
    return p;
}


