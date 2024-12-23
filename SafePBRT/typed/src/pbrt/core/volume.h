
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

#ifndef PBRT_CORE_VOLUME_H
#define PBRT_CORE_VOLUME_H

// core/volume.h*
#include "pbrt.h"
#include "spectrum.h"
#include "geometry.h"
#include "transform.h"
#include "integrator.h"

 //Volume Scattering Definitions
template<typename S>
mfloat<invsolidangle_d> PhaseIsotropic(const Direction<S> &, const Direction<S> &) {
    return 1.f / sphericalAngle;
}

template<typename S>
mfloat<invsolidangle_d> PhaseRayleigh(const Direction<S> &w, const Direction<S> &wp) {
    float costheta = Dot(w, wp);
    return  3.f/(4.f * sphericalAngle) * (1 + costheta * costheta);
}

template<typename S>
mfloat<invsolidangle_d> PhaseMieHazy(const Direction<S> &w, const Direction<S> &wp) {
    float costheta = Dot(w, wp);
    return (0.5f + 4.5f * powf(0.5 * (1.f + costheta), 8.f)) / (sphericalAngle);
}

template<typename S>
mfloat<invsolidangle_d> PhaseMieMurky(const Direction<S> &w, const Direction<S> &wp) {
    float costheta = Dot(w, wp);
    return (0.5f + 16.5f * powf(0.5 * (1.f + costheta), 32.f)) / (sphericalAngle);
}

template<typename S>
mfloat<invsolidangle_d> PhaseHG(const Direction<S> &w, const Direction<S> &wp, float g) {
    float costheta = Dot(w, wp);
    return 1.f / (sphericalAngle) *
        (1.f - g*g) / powf(1.f + g*g - 2.f * g * costheta, 1.5f);
}

template<typename S>
mfloat<invsolidangle_d> PhaseSchlick(const Direction<S> &w, const Direction<S> &wp, float g) {
    float k = 1.55f * g - .55f * g * g * g;
    float kcostheta = k * Dot(w, wp);
    return 1.f / (sphericalAngle) *
        (1.f - k*k) / ((1.f - kcostheta) * (1.f - kcostheta));
}

class VolumeRegion {
public:
    // VolumeRegion Interface
    virtual ~VolumeRegion();
    virtual BBox<world_s> WorldBound() const = 0;
    virtual bool IntersectP(const Ray<world_s> &ray, mfloat<length_d> *t0, mfloat<length_d> *t1) const = 0;
    virtual Spectrum<invlength_d> sigma_a(const Point<world_s> &, const Direction<world_s> &,
                            const mfloat<time_d> &time) const = 0;
    virtual Spectrum<invlength_d> sigma_s(const Point<world_s> &, const Direction<world_s> &,
                            const mfloat<time_d> &time) const = 0;
    virtual Spectrum<radiancediff_d> Lve(const Point<world_s> &, const Direction<world_s> &,
                         const mfloat<time_d> &time) const = 0;
    virtual mfloat<invsolidangle_d> p(const Point<world_s> &, const Direction<world_s> &,
                    const Direction<world_s> &, const mfloat<time_d> &time) const = 0;
    virtual Spectrum<invlength_d> sigma_t(const Point<world_s> &p, const Direction<world_s> &wo, const mfloat<time_d> & time) const;
    virtual Spectrum<rho_d> tau(const Ray<world_s> &ray, const mfloat<length_d> &step = mfloat<length_d>(1.f),
                         float offset = 0.5f) const = 0;
};


class DensityRegion : public VolumeRegion {
public:
    // DensityRegion Public Methods
    DensityRegion(const Spectrum<invlength_d> &sa, const Spectrum<invlength_d> &ss, float gg,
                  const Spectrum<radiancediff_d> &emit, const Transform<volume_s, world_s> &VolumeToWorld)
        : sig_a(sa), sig_s(ss), le(emit), g(gg),
          WorldToVolume(Inverse(VolumeToWorld)) { }
    virtual float Density(const Point<volume_s> &Pobj) const = 0;
    Spectrum<invlength_d> sigma_a(const Point<world_s> &p, const Direction<world_s> &, const mfloat<time_d> &) const {
        return Density(WorldToVolume(p)) * sig_a;
    }
    Spectrum<invlength_d> sigma_s(const Point<world_s> &p, const Direction<world_s> &, const mfloat<time_d> &) const {
        return Density(WorldToVolume(p)) * sig_s;
    }
    Spectrum<invlength_d> sigma_t(const Point<world_s> &p, const Direction<world_s> &, const mfloat<time_d> &) const {
        return Density(WorldToVolume(p)) * (sig_a + sig_s);
    }
    Spectrum<radiancediff_d> Lve(const Point<world_s> &p, const Direction<world_s> &, const mfloat<time_d> &) const {
        return Density(WorldToVolume(p)) * le;
    }
    mfloat<invsolidangle_d> p(const Point<world_s> &p, const Direction<world_s> &w, const Direction<world_s> &wp, const mfloat<time_d> &) const {
        return PhaseHG(w, wp, g);
    }
    Spectrum<rho_d> tau(const Ray<world_s> &r, const mfloat<length_d> &stepSize, float offset) const;
protected:
    // DensityRegion Protected Data
    Spectrum<invlength_d> sig_a, sig_s;
    Spectrum<radiancediff_d> le;
    float g;
    Transform<world_s, volume_s> WorldToVolume;
};


class AggregateVolume : public VolumeRegion {
public:
    // AggregateVolume Public Methods
    AggregateVolume(const std::vector<VolumeRegion *> &r);
    ~AggregateVolume();
    BBox<world_s> WorldBound() const;
    bool IntersectP(const Ray<world_s> &ray, mfloat<length_d> *t0, mfloat<length_d> *t1) const;
    Spectrum<invlength_d> sigma_a(const Point<world_s> &, const Direction<world_s> &, const mfloat<time_d> &) const;
    Spectrum<invlength_d> sigma_s(const Point<world_s> &, const Direction<world_s> &, const mfloat<time_d> &) const;
    Spectrum<radiancediff_d> Lve(const Point<world_s> &, const Direction<world_s> &, const mfloat<time_d> &) const;
    mfloat<invsolidangle_d> p(const Point<world_s> &, const Direction<world_s> &, const Direction<world_s> &, const mfloat<time_d> &) const;
    Spectrum<invlength_d> sigma_t(const Point<world_s> &, const Direction<world_s> &, const mfloat<time_d> &) const;
    Spectrum<rho_d> tau(const Ray<world_s> &ray, const mfloat<length_d> &, float) const;
private:
    // AggregateVolume Private Data
    std::vector<VolumeRegion *> regions;
    BBox<world_s> bound;
};


//bool GetVolumeScatteringProperties(const string &name, Spectrum *sigma_a,
//                                   Spectrum *sigma_prime_s);
class VolumeIntegrator : public Integrator {
public:
    // VolumeIntegrator Interface
    virtual Spectrum<radiance_d> Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential<world_s> &ray, const Sample *sample, RNG &rng,
        Spectrum<rho_d> *transmittance, MemoryArena &arena) const = 0;
    virtual Spectrum<rho_d> Transmittance(const Scene *scene,
        const Renderer *renderer, const RayDifferential<world_s> &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const = 0;
};


//void SubsurfaceFromDiffuse(const Spectrum<rho_d> &Kd, const mfloat<length_d> &meanPathLength, float eta,
//        Spectrum<invlength_d> *sigma_a, Spectrum<invlength_d> *sigma_prime_s);

#endif // PBRT_CORE_VOLUME_H
