
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


// cameras/perspective.cpp*
#include "cameras/perspective.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"

// PerspectiveCamera Method Definitions
PerspectiveCamera:: PerspectiveCamera(const AnimatedTransform<camera_s, world_s> &cam2world,
                                      const mfloat<length_d> screenWindow[4], const mfloat<time_d> &sopen, const mfloat<time_d> &sclose, 
                                      const mfloat<length_d> &lensr, const mfloat<length_d> &focald, float fov, Film *f)
                                      : ProjectiveCamera(cam2world, Perspective(fov, 1e-2f*meters, 1000.f*meters),
                                      screenWindow, sopen, sclose, lensr, focald, f) 
{
    // Compute differential changes in origin for perspective camera rays
    dxCamera = RasterToCamera(Point<raster_s>(1*meters,0*meters,0*meters)) - RasterToCamera(Point<raster_s>());
    dyCamera = RasterToCamera(Point<raster_s>(0*meters,1*meters,0*meters)) - RasterToCamera(Point<raster_s>());
}


float PerspectiveCamera::GenerateRay(const CameraSample &sample,
                                     Ray<world_s> *ray) const 
{
    // Generate raster and camera samples
    Point<raster_s> Pras(sample.imageX, sample.imageY, mfloat<length_d>(0.0f));
    Point<camera_s> Pcamera;
    RasterToCamera(Pras, &Pcamera);
    Ray<camera_s> r(Point<camera_s>(), Normalize(AsVector(Pcamera)), 0.f*meters, INFINITY*meters);
    // Modify ray for depth of field
    if (lensRadius > mfloat<length_d>(0.0f)) {
        // Sample point on lens
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);

        // Compute point on plane of focus
        mfloat<length_d> ft = focalDistance / r.d.z;
        Point<camera_s> Pfocus = r(ft);

        // Update ray for effect of lens
        r.o = Point<camera_s>(lensU * lensRadius, lensV * lensRadius, 0.0f*meters);
        r.d = Normalize(Pfocus - r.o);
    }
    r.time = Lerp(__asfloat(sample.time), shutterOpen, shutterClose);
    CameraToWorld(r, ray);
    return 1.0f;
}


float PerspectiveCamera::GenerateRayDifferential(const CameraSample &sample,
                                                 RayDifferential<world_s> *ray) const 
{
    // Generate raster and camera samples
    Point<raster_s> Pras(sample.imageX, sample.imageY, mfloat<length_d>(0.0f));
    Point<camera_s> Pcamera;
    RasterToCamera(Pras, &Pcamera);
    Direction<camera_s> dir = Normalize(Vector<camera_s>(Pcamera.x, Pcamera.y, Pcamera.z));
    RayDifferential<camera_s> r(Point<camera_s>(), dir, 0.f*meters, INFINITY*meters);
    // Modify ray for depth of field
    if (lensRadius > mfloat<length_d>(0.0)) {
        // Sample point on lens
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);

        // Compute point on plane of focus
        mfloat<length_d> ft = focalDistance / r.d.z;
        Point<camera_s> Pfocus = r(ft);

        // Update ray for effect of lens
        r.o = Point<camera_s>(lensU * lensRadius, lensV * lensRadius, 0.0f * meters);
        r.d = Normalize(Pfocus - r.o);
    }

    // Compute offset rays for \use{PerspectiveCamera} ray differentials
    r.rxOrigin = r.ryOrigin = r.o;
    r.rxDirection = Normalize(Vector<camera_s>(Pcamera) + dxCamera);
    r.ryDirection = Normalize(Vector<camera_s>(Pcamera) + dyCamera);
    r.time = Lerp(__asfloat(sample.time), shutterOpen, shutterClose);
    r.hasDifferentials = true;
    CameraToWorld(r, ray);
    return 1.f;
}


PerspectiveCamera *CreatePerspectiveCamera(const ParamSet &params,
                                           const AnimatedTransform<camera_s, world_s> &cam2world, Film *film) 
{
    // Extract common camera parameters from _ParamSet_
    mfloat<time_d> shutteropen = params.FindOneFloat<time_d>("shutteropen", mfloat<time_d>(0.f));
    mfloat<time_d> shutterclose = params.FindOneFloat<time_d>("shutterclose", mfloat<time_d>(1.f));
    mfloat<length_d> lensradius = params.FindOneFloat<length_d>("lensradius", mfloat<length_d>(0.f));
    mfloat<length_d> focaldistance = params.FindOneFloat<length_d>("focaldistance", mfloat<length_d>(1e30f));
    float frame = params.FindOneFloat("frameaspectratio",
        float(film->xResolution)/float(film->yResolution));
    mfloat<length_d> screen[4];
    if (frame > 1.f) {
        screen[0] = mfloat<length_d>(-frame);
        screen[1] = mfloat<length_d>(frame);
        screen[2] = mfloat<length_d>(-1.f);
        screen[3] = mfloat<length_d>(1.f);
    }
    else {
        screen[0] = mfloat<length_d>(-1.f);
        screen[1] = mfloat<length_d>(1.f);
        screen[2] = mfloat<length_d>(-1.f / frame);
        screen[3] = mfloat<length_d>(1.f / frame);
    }
    int swi;
    const mfloat<length_d> *sw = params.FindFloat<length_d>("screenwindow", &swi);
    if (sw && swi == 4)
        memcpy(screen, sw, 4*sizeof(mfloat<length_d>));
    float fov = params.FindOneFloat("fov", 90.);
    return new PerspectiveCamera(cam2world, screen, shutteropen,
        shutterclose, lensradius, focaldistance, fov, film);
}


