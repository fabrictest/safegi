
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

#ifndef PBRT_CORE_CAMERA_H
#define PBRT_CORE_CAMERA_H

// core/camera.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"

// Camera Declarations
class Camera {
public:
    // Camera Interface
    Camera(const AnimatedTransform<camera_s, world_s> &cam2world, const mfloat<time_d> &sopen, const mfloat<time_d> &sclose,
           Film *film);
    virtual ~Camera();
    virtual float GenerateRay(const CameraSample &sample,
                              Ray<world_s> *ray) const = 0;
    virtual float GenerateRayDifferential(const CameraSample &sample, RayDifferential<world_s> *rd) const;

    // Camera Public Data
    AnimatedTransform<camera_s, world_s> CameraToWorld;
    const mfloat<time_d> shutterOpen, shutterClose;
    Film *film;
};


class ProjectiveCamera : public Camera {
public:
    // ProjectiveCamera Public Methods
    ProjectiveCamera(const AnimatedTransform<camera_s, world_s> &cam2world,
        const Transform<camera_s, screen_s> &proj, const mfloat<length_d> screenWindow[4],
        const mfloat<time_d> &sopen, const mfloat<time_d> &sclose, 
        const mfloat<length_d> &lensr, const mfloat<length_d> &focald, Film *film);
protected:
    // ProjectiveCamera Protected Data
    Transform<camera_s, screen_s> CameraToScreen;
    Transform<raster_s, camera_s> RasterToCamera;
    Transform<screen_s, raster_s> ScreenToRaster;
    Transform<raster_s, screen_s> RasterToScreen;
    mfloat<length_d> lensRadius, focalDistance;
};



#endif // PBRT_CORE_CAMERA_H
