
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


// core/scene.cpp*
#include "scene.h"
#include "camera.h"
#include "film.h"
#include "sampler.h"
#include "volume.h"
#include "parallel.h"
#include "progressreporter.h"
#include "renderer.h"
#include "light.h"

// Scene Method Definitions
Scene::~Scene() {
    delete aggregate;
    //delete volumeRegion;
    for (uint32_t i = 0; i < lights.size(); ++i)
        delete lights[i];
}


Scene::Scene(Primitive *accel, const std::vector<Light *> &lts,
             VolumeRegion *vr) {
    lights = lts;
    aggregate = accel;
    volumeRegion = vr;
    if (lts.size() == 0)
        Warning("No light sources defined in scene; "
            "possibly rendering a black image.");
    // Scene Constructor Implementation
    bound = aggregate->WorldBound();
    //if (volumeRegion) bound = Union(bound, volumeRegion->WorldBound());
}


const BBox<world_s> &Scene::WorldBound() const {
    return bound;
}


