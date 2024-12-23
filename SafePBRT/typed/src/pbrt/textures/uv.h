
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

#ifndef PBRT_TEXTURES_UV_H
#define PBRT_TEXTURES_UV_H

// textures/uv.h*
#include "pbrt.h"
#include "mtexture.h"
#include "paramset.h"

// UVTexture Declarations
template<typename D>
class UVTexture : public Texture<Spectrum<D> > {
public:
    // UVTexture Public Methods
    UVTexture(TextureMapping2D *m) {
        mapping = m;
    }
    ~UVTexture() {
        delete mapping;
    }
    Spectrum<D> Evaluate(const DifferentialGeometry &dg) const {
        float s, t, dsdx, dtdx, dsdy, dtdy;
        mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
        float rgb[3] = { s - Floor2Int(s), t - Floor2Int(t), 0.f };
        return Spectrum<D>::FromRGB(rgb);
    }
private:
    TextureMapping2D *mapping;
};


inline Texture<mfloat<unit_d> > *CreateUVFloatTexture(const Transform<texture_s, world_s> &tex2world,
        const TextureParams &tp)
{
    return NULL;
}

UVTexture<Spectrum<rho_d> > *CreateUVSpectrumTexture(const Transform<texture_s, world_s> &tex2world,
                                      const TextureParams &tp) 
{
    // Initialize 2D texture mapping _map_ from _tp_
    TextureMapping2D *map = NULL;
    string type = tp.FindString("mapping", "uv");
    if (type == "uv") {
        float su = tp.FindFloat("uscale", 1.);
        float sv = tp.FindFloat("vscale", 1.);
        float du = tp.FindFloat("udelta", 0.);
        float dv = tp.FindFloat("vdelta", 0.);
        map = new UVMapping2D(su, sv, du, dv);
    }
    else if (type == "spherical") map = new SphericalMapping2D(Inverse(tex2world));
    else if (type == "cylindrical") map = new CylindricalMapping2D(Inverse(tex2world));
    else if (type == "planar")
        map = new PlanarMapping2D(tp.FindVector("v1", Vector<world_s>(1*meters,0*meters,0*meters)),
        tp.FindVector("v2", Vector<world_s>(0*meters,1*meters,0*meters)),
        tp.FindFloat("udelta", 0.f), tp.FindFloat("vdelta", 0.f));
    else {
        Error("2D texture mapping \"%s\" unknown", type.c_str());
        map = new UVMapping2D();
    }
    return new UVTexture<Spectrum<rho_d> >(map);
}


#endif // PBRT_TEXTURES_UV_H
