
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


// core/camera.cpp*
#include "pbrt.h"
#include "camera.h"
#include "film.h"
#include "montecarlo.h"
#include "sampler.h"

// Camera Method Definitions
Camera::~Camera() {
    delete film;
}


Camera::Camera(const AnimatedTransform<camera_s, world_s> &cam2world,
               const mfloat<time_d> &sopen, const mfloat<time_d> &sclose, Film *f)
    : CameraToWorld(cam2world), shutterOpen(sopen), shutterClose(sclose) {
    film = f;
    if (CameraToWorld.HasScale())
        Warning("Scaling detected in world-to-camera transformation!\n"
                "The system has numerous assumptions, implicit and explicit,\n"
                "that this transform will have no scale factors in it.\n"
                "Proceed at your own risk; your image may have errors or\n"
                "the system may crash as a result of this.");
}


float Camera::GenerateRayDifferential(const CameraSample &sample,
                                      RayDifferential<world_s> *rd) const {
    float wt = GenerateRay(sample, rd);
    // Find ray after shifting one pixel in the $x$ direction
    CameraSample sshift = sample;
    (sshift.imageX) += 1*meters;
    Ray<world_s> rx;
    float wtx = GenerateRay(sshift, &rx);
    rd->rxOrigin = rx.o;
    rd->rxDirection = rx.d;

    // Find ray after shifting one pixel in the $y$ direction
    (sshift.imageX) -= 1*meters;
    (sshift.imageY) += 1*meters;
    Ray<world_s> ry;
    float wty = GenerateRay(sshift, &ry);
    rd->ryOrigin = ry.o;
    rd->ryDirection = ry.d;
    if (wtx == 0.f || wty == 0.f) return 0.f;
    rd->hasDifferentials = true;
    return wt;
}


ProjectiveCamera::ProjectiveCamera(const AnimatedTransform<camera_s, world_s> &cam2world,
        const Transform<camera_s, screen_s> &proj, const mfloat<length_d> screenWindow[4], 
        const mfloat<time_d> &sopen, const mfloat<time_d> &sclose, 
        const mfloat<length_d> &lensr, const mfloat<length_d> &focald, Film *f)
    : Camera(cam2world, sopen, sclose, f) {
    // Initialize depth of field parameters
    lensRadius = lensr;
    focalDistance = focald;

    // Compute projective camera transformations
    CameraToScreen = proj;

    // Compute projective camera screen transformations
    ScreenToRaster = 
        Scale<screen_s, raster_s>(mfloat<length_d>(film->xResolution) / (screenWindow[1] - screenWindow[0]),
              mfloat<length_d>(film->yResolution) / (screenWindow[2] - screenWindow[3]), 1.f) *
        Translate<screen_s, screen_s>(Vector<screen_s>(-screenWindow[0], -screenWindow[3], mfloat<length_d>(0.0f)));
    RasterToScreen = Inverse(ScreenToRaster);
    RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;
}


