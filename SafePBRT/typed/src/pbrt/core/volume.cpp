
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


// core/volume.cpp*
#include "volume.h"

// Volume Scattering Local Definitions
struct MeasuredSS {
    const char *name;
    float sigma_prime_s[3], sigma_a[3]; // mm^-1
};


static MeasuredSS mss[] = {
// From "A Practical Model for Subsurface Light Transport"
// Jensen, Marschner, Levoy, Hanrahan
// Proc SIGGRAPH 2001
{ "Apple", { 2.29, 2.39, 1.97 }, { 0.0030, 0.0034, 0.046 }, },
{ "Chicken1", { 0.15, 0.21, 0.38 }, { 0.015, 0.077, 0.19 }, },
{ "Chicken2", { 0.19, 0.25, 0.32 }, { 0.018, 0.088, 0.20 }, },
{ "Cream", { 7.38, 5.47, 3.15 }, { 0.0002, 0.0028, 0.0163 }, },
{ "Ketchup", { 0.18, 0.07, 0.03 }, { 0.061, 0.97, 1.45 }, },
{ "Marble", { 2.19, 2.62, 3.00 }, { 0.0021, 0.0041, 0.0071 }, },
{ "Potato", { 0.68, 0.70, 0.55 }, { 0.0024, 0.0090, 0.12 }, },
{ "Skimmilk", { 0.70, 1.22, 1.90 }, { 0.0014, 0.0025, 0.0142 }, },
{ "Skin1", { 0.74, 0.88, 1.01 }, { 0.032, 0.17, 0.48 }, },
{ "Skin2", { 1.09, 1.59, 1.79 }, { 0.013, 0.070, 0.145 }, },
{ "Spectralon", { 11.6, 20.4, 14.9 }, { 0.00, 0.00, 0.00 }, },
{ "Wholemilk", { 2.55, 3.21, 3.77 }, { 0.0011, 0.0024, 0.014 }, },

// From "Acquiring Scattering Properties of Participating Media by Dilution",
// Narasimhan, Gupta, Donner, Ramamoorthi, Nayar, Jensen
// Proc SIGGRAPH 2006
{ "Lowfat Milk", { 0.912600, 1.074800, 1.250000 }, { 0.000200, 0.000400, 0.000800 } },
{ "Reduced Milk", { 1.075000, 1.221300, 1.394100 }, { 0.000200, 0.000400, 0.001000 } },
{ "Regular Milk", { 1.187400, 1.329600, 1.460200 }, { 0.000100, 0.000300, 0.001300 } },
{ "Espresso", { 0.437600, 0.511500, 0.604800 }, { 0.166900, 0.228700, 0.307800 } },
{ "Mint Mocha Coffee", { 0.190000, 0.260000, 0.350000 }, { 0.098400, 0.151900, 0.204000 } },
{ "Lowfat Soy Milk", { 0.141900, 0.162500, 0.274000 }, { 0.000100, 0.000500, 0.002500 } },
{ "Regular Soy Milk", { 0.243400, 0.271900, 0.459700 }, { 0.000100, 0.000500, 0.003400 } },
{ "Lowfat Chocolate Milk", { 0.428200, 0.501400, 0.579100 }, { 0.000500, 0.001600, 0.006800 } },
{ "Regular Chocolate Milk", { 0.735900, 0.917200, 1.068800 }, { 0.000700, 0.003000, 0.010000 } },
{ "Coke", { 0.714300, 1.168800, 1.716900 }, { 0.696600, 1.148000, 1.716900 } },
{ "Pepsi", { 0.643300, 0.999000, 1.442000 }, { 0.637500, 0.984900, 1.442000 } },
{ "Sprite", { 0.129900, 0.128300, 0.139500 }, { 0.123000, 0.119400, 0.130600 } },
{ "Gatorade", { 0.400900, 0.418500, 0.432400 }, { 0.161700, 0.125800, 0.057900 } },
{ "Chardonnay", { 0.157700, 0.174800, 0.351200 }, { 0.154700, 0.170100, 0.344300 } },
{ "White Zinfandel", { 0.176300, 0.237000, 0.291300 }, { 0.173200, 0.232200, 0.284700 } },
{ "Merlot", { 0.763900, 1.642900, 1.919600 }, { 0.758600, 1.642900, 1.919600 } },
{ "Budweiser Beer", { 0.148600, 0.321000, 0.736000 }, { 0.144900, 0.314100, 0.728600 } },
{ "Coors Light Beer", { 0.029500, 0.066300, 0.152100 }, { 0.026800, 0.060800, 0.152100 } },
{ "Clorox", { 0.160000, 0.250000, 0.330000 }, { 0.017500, 0.077700, 0.137200 } },
{ "Apple Juice", { 0.121500, 0.210100, 0.440700 }, { 0.101400, 0.185800, 0.408400 } },
{ "Cranberry Juice", { 0.270000, 0.630000, 0.830000 }, { 0.257200, 0.614500, 0.810400 } },
{ "Grape Juice", { 0.550000, 1.250000, 1.530000 }, { 0.542800, 1.250000, 1.530000 } },
{ "Ruby Grapefruit Juice", { 0.251300, 0.351700, 0.430500 }, { 0.089600, 0.191100, 0.263600 } },
{ "White Grapefruit Juice", { 0.360900, 0.380000, 0.563200 }, { 0.009600, 0.013100, 0.039500 } },
{ "Shampoo", { 0.028800, 0.071000, 0.095200 }, { 0.018400, 0.059600, 0.080500 } },
{ "Strawberry Shampoo", { 0.021700, 0.078800, 0.102200 }, { 0.018900, 0.075600, 0.098900 } },
{ "Head & Shoulders Shampoo", { 0.367400, 0.452700, 0.521100 }, { 0.088300, 0.163700, 0.212500 } },
{ "Lemon Tea", { 0.340000, 0.580000, 0.880000 }, { 0.260200, 0.490200, 0.772700 } },
{ "Orange Juice Powder", { 0.337700, 0.557300, 1.012200 }, { 0.144900, 0.344100, 0.786300 } },
{ "Pink Lemonade", { 0.240000, 0.370000, 0.450000 }, { 0.116500, 0.236600, 0.319500 } },
{ "Cappuccino Powder", { 0.257400, 0.353600, 0.484000 }, { 0.192000, 0.265400, 0.327200 } },
{ "Salt Powder", { 0.760000, 0.868500, 0.936300 }, { 0.511500, 0.586300, 0.614700 } },
{ "Sugar Powder", { 0.079500, 0.175900, 0.278000 }, { 0.065000, 0.159700, 0.257800 } },
{ "Suisse Mocha", { 0.509800, 0.647600, 0.794400 }, { 0.187500, 0.289300, 0.379600 } },
{ "Pacific Ocean Surface Water", { 3.364500, 3.315800, 3.242800 }, { 3.184500, 3.132400, 3.014700 } },

};


static float RdIntegral(float alphap, float A) {
    float sqrtTerm = sqrtf(3.f * (1.f - alphap));
    return alphap / 2.f * (1.f + expf(-4.f/3.f * A * sqrtTerm)) *
        expf(-sqrtTerm);
}


static float RdToAlphap(float reflectance, float A) {
    float alphaLow = 0., alphaHigh = 1.f;
    float kd0 = RdIntegral(alphaLow, A);
    float kd1 = RdIntegral(alphaHigh, A);
    for (int i = 0; i < 16; ++i) {
        Assert(kd0 <= reflectance && kd1 >= reflectance);
        float alphaMid = (alphaLow + alphaHigh) * 0.5f;
        float kd = RdIntegral(alphaMid, A);
        if (kd < reflectance) { alphaLow = alphaMid;  kd0 = kd; }
        else                  { alphaHigh = alphaMid; kd1 = kd; }
    }
    return (alphaLow + alphaHigh) * 0.5f;
}






VolumeRegion::~VolumeRegion() {
}


Spectrum<invlength_d> VolumeRegion::sigma_t(const Point<world_s> &p, const Direction<world_s> &w,
                               const mfloat<time_d> & time) const {
    return sigma_a(p, w, time) + sigma_s(p, w, time);
}


AggregateVolume::AggregateVolume(const std::vector<VolumeRegion *> &r) {
    regions = r;
    for (uint32_t i = 0; i < regions.size(); ++i)
        bound = Union(bound, regions[i]->WorldBound());
}


Spectrum<invlength_d> AggregateVolume::sigma_a(const Point<world_s> &p, const Direction<world_s> &w,
                                  const mfloat<time_d> & time) const {
    Spectrum<invlength_d> s;
    for (uint32_t i = 0; i < regions.size(); ++i)
        s += regions[i]->sigma_a(p, w, time);
    return s;
}


Spectrum<invlength_d> AggregateVolume::sigma_s(const Point<world_s> &p, const Direction<world_s> &w, const mfloat<time_d> & time) const {
    Spectrum<invlength_d> s;
    for (uint32_t i = 0; i < regions.size(); ++i)
        s += regions[i]->sigma_s(p, w, time);
    return s;
}


Spectrum<radiancediff_d> AggregateVolume::Lve(const Point<world_s> &p, const Direction<world_s> &w, const mfloat<time_d> & time) const {
    Spectrum<radiancediff_d> L;
    for (uint32_t i = 0; i < regions.size(); ++i)
        L += regions[i]->Lve(p, w, time);
    return L;
}


mfloat<invsolidangle_d> AggregateVolume::p(const Point<world_s> &p, const Direction<world_s> &w, const Direction<world_s> &wp,
        const mfloat<time_d> & time) const {
    mfloat<invsolidangle_d> ph(0.0f);
    mfloat<rho_d> sumWt = 0;
    for (uint32_t i = 0; i < regions.size(); ++i) {
        mfloat<invlength_d> sigt = regions[i]->sigma_t(p, w, time).y();
        if (sigt != mfloat<invlength_d>(0.0f)) {
            mfloat<rho_d> wt = regions[i]->sigma_s(p, w, time).y() / sigt;
            sumWt += wt;
            ph += wt * regions[i]->p(p, w, wp, time);
        }
    }
    return ph / sumWt;
}


Spectrum<invlength_d> AggregateVolume::sigma_t(const Point<world_s> &p, const Direction<world_s> &w, const mfloat<time_d> & time) const {
    Spectrum<invlength_d> s;
    for (uint32_t i = 0; i < regions.size(); ++i)
        s += regions[i]->sigma_t(p, w, time);
    return s;
}


Spectrum<rho_d> AggregateVolume::tau(const Ray<world_s> &ray, const mfloat<length_d> &step, float offset) const {
    Spectrum<rho_d> t;
    for (uint32_t i = 0; i < regions.size(); ++i)
        t += regions[i]->tau(ray, step, offset);
    return t;
}


bool AggregateVolume::IntersectP(const Ray<world_s> &ray,
                                 mfloat<length_d> *t0, mfloat<length_d> *t1) const {
    *t0 = mfloat<length_d>(INFINITY);
    *t1 = mfloat<length_d>(-INFINITY);
    for (uint32_t i = 0; i < regions.size(); ++i) {
        mfloat<length_d> tr0, tr1;
        if (regions[i]->IntersectP(ray, &tr0, &tr1)) {
            *t0 = min(*t0, tr0);
            *t1 = max(*t1, tr1);
        }
    }
    return (*t0 < *t1);
}


AggregateVolume::~AggregateVolume() {
    for (uint32_t i = 0; i < regions.size(); ++i)
        delete regions[i];
}


BBox<world_s> AggregateVolume::WorldBound() const {
    return bound;
}

Spectrum<rho_d> DensityRegion::tau(const Ray<world_s> &r, const mfloat<length_d> &stepSize,
                                   float u) const 
{
    mfloat<length_d> t0, t1;
    Ray<world_s> rn(r.o, r.d, r.mint, r.maxt, r.time);
    if (!IntersectP(rn, &t0, &t1)) return Spectrum<rho_d>();
    Spectrum<invlength_d> tau;
    t0 += u * stepSize;
    while (t0 < t1) {
        tau += sigma_t(rn(t0), -rn.d, r.time);
        t0 += stepSize;
    }
    return tau * stepSize;
}


