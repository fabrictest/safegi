
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

#ifndef PBRT_VOLUMES_HOMOGENEOUS_H
#define PBRT_VOLUMES_HOMOGENEOUS_H

// volumes/homogeneous.h*
#include "volume.h"

// HomogeneousVolumeDensity Declarations
class HomogeneousVolumeDensity : public VolumeRegion {
public:
    // HomogeneousVolumeDensity Public Methods
    HomogeneousVolumeDensity(const Spectrum<invlength_d> &sa, const Spectrum<invlength_d> &ss, float gg,
            const Spectrum<radiancediff_d> &emit, const BBox<volume_s> &e, const Transform<volume_s, world_s> &v2w) {
        WorldToVolume = Inverse(v2w);
        sig_a = sa;
        sig_s = ss;
        g = gg;
        le = emit;
        extent = e;
    }
    BBox<world_s> WorldBound() const {
        return Inverse(WorldToVolume)(extent);
    }
    bool IntersectP(const Ray<world_s> &r, mfloat<length_d> *t0, mfloat<length_d> *t1) const {
        Ray<volume_s> ray = WorldToVolume(r);
        return extent.IntersectP(ray, t0, t1);
    }
    Spectrum<invlength_d> sigma_a(const Point<world_s> &p, const Direction<world_s> &, const mfloat<time_d> &) const {
        return extent.Inside(WorldToVolume(p)) ? sig_a : Spectrum<invlength_d>();
    }
    Spectrum<invlength_d> sigma_s(const Point<world_s> &p, const Direction<world_s> &, const mfloat<time_d> &) const {
        return extent.Inside(WorldToVolume(p)) ? sig_s : Spectrum<invlength_d>();
    }
    Spectrum<invlength_d> sigma_t(const Point<world_s> &p, const Direction<world_s> &, const mfloat<time_d> &) const {
        return extent.Inside(WorldToVolume(p)) ? (sig_a + sig_s) : Spectrum<invlength_d>();
    }
    Spectrum<radiancediff_d> Lve(const Point<world_s> &p, const Direction<world_s> &, const mfloat<time_d> &) const {
        return extent.Inside(WorldToVolume(p)) ? le : Spectrum<radiancediff_d>();
    }
    mfloat<invsolidangle_d> p(const Point<world_s> &p, const Direction<world_s> &wi, const Direction<world_s> &wo, const mfloat<time_d> &) const {
        if (!extent.Inside(WorldToVolume(p))) return mfloat<invsolidangle_d>(0.0f);
        return PhaseHG(wi, wo, g);
    }
    Spectrum<unit_d> tau(const Ray<world_s> &ray, const mfloat<length_d> &, float) const {
        mfloat<length_d> t0, t1;
        if (!IntersectP(ray, &t0, &t1)) return Spectrum<rho_d>();
        return Distance(ray(t0), ray(t1)) * (sig_a + sig_s);
    }
private:
    // HomogeneousVolumeDensity Private Data
    Spectrum<invlength_d> sig_a, sig_s;
    Spectrum<radiancediff_d> le;
    float g;
    BBox<volume_s> extent;
    Transform<world_s, volume_s> WorldToVolume;
};


HomogeneousVolumeDensity *CreateHomogeneousVolumeDensityRegion(const Transform<volume_s, world_s> &volume2world,
        const ParamSet &params);

#endif // PBRT_VOLUMES_HOMOGENEOUS_H
