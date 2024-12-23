
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


// volumes/homogeneous.cpp*
#include "volumes/homogeneous.h"
#include "paramset.h"

// HomogeneousVolumeDensity Method Definitions
HomogeneousVolumeDensity *CreateHomogeneousVolumeDensityRegion(const Transform<volume_s, world_s> &volume2world,
        const ParamSet &params) {
    // Initialize common volume region parameters
    Spectrum<invlength_d> sigma_a = params.FindOneSpectrum<invlength_d>("sigma_a", Spectrum<invlength_d>());
    Spectrum<invlength_d> sigma_s = params.FindOneSpectrum<invlength_d>("sigma_s", Spectrum<invlength_d>());
    float g = params.FindOneFloat("g", 0.);
    Spectrum<radiancediff_d> Le = params.FindOneSpectrum<radiancediff_d>("Le", Spectrum<radiancediff_d>());
    Point<volume_s> p0 = params.FindOnePoint("p0", Point<volume_s>(0*meters,0*meters,0*meters));
    Point<volume_s> p1 = params.FindOnePoint("p1", Point<volume_s>(1*meters,1*meters,1*meters));
    return new HomogeneousVolumeDensity(sigma_a, sigma_s, g, Le, BBox<volume_s>(p0, p1),
        volume2world);
}


