
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


// textures/imagemap.cpp*
#include "textures/imagemap.h"

// ImageTexture Method Definitions

ImageTexture<mfloat<unit_d>, mfloat<unit_d>> *CreateImageFloatTexture(const Transform<texture_s, world_s> &tex2world,
        const TextureParams &tp) {
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
        map = new PlanarMapping2D(
        Normalize(tp.FindVector("v1", Vector<world_s>(1*meters,0*meters,0*meters))),
            Normalize(tp.FindVector("v2", Vector<world_s>(0*meters,1*meters,0*meters))),
            tp.FindFloat("udelta", 0.f), tp.FindFloat("vdelta", 0.f));
    else {
        Error("2D texture mapping \"%s\" unknown", type.c_str());
        map = new UVMapping2D;
    }

    // Initialize _ImageTexture_ parameters
    float maxAniso = tp.FindFloat("maxanisotropy", 8.f);
    bool trilerp = tp.FindBool("trilinear", false);
    string wrap = tp.FindString("wrap", "repeat");
    ImageWrap wrapMode = TEXTURE_REPEAT;
    if (wrap == "black") wrapMode = TEXTURE_BLACK;
    else if (wrap == "clamp") wrapMode = TEXTURE_CLAMP;
    float scale = tp.FindFloat("scale", 1.f);
    float gamma = tp.FindFloat("gamma", 1.f);
    return new ImageTexture<mfloat<unit_d>, mfloat<unit_d> >(map, tp.FindString("filename"),
        trilerp, maxAniso, wrapMode, scale, gamma);
}

ImageTexture<RGBSpectrum<rho_d>, Spectrum<rho_d>> *CreateImageSpectrumTexture(const Transform<texture_s, world_s> &tex2world,
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
        map = new PlanarMapping2D(
        Normalize(tp.FindVector("v1", Vector<world_s>(mfloat<length_d>(1.0f),mfloat<length_d>(0.0f),mfloat<length_d>(0.0f)))),
        Normalize(tp.FindVector("v2", Vector<world_s>(mfloat<length_d>(0.0f),mfloat<length_d>(1.0f),mfloat<length_d>(0.0f)))),
        tp.FindFloat("udelta", 0.f), tp.FindFloat("vdelta", 0.f));
    else {
        Error("2D texture mapping \"%s\" unknown", type.c_str());
        map = new UVMapping2D;
    }

    // Initialize _ImageTexture_ parameters
    float maxAniso = tp.FindFloat("maxanisotropy", 8.f);
    bool trilerp = tp.FindBool("trilinear", false);
    string wrap = tp.FindString("wrap", "repeat");
    ImageWrap wrapMode = TEXTURE_REPEAT;
    if (wrap == "black") wrapMode = TEXTURE_BLACK;
    else if (wrap == "clamp") wrapMode = TEXTURE_CLAMP;
    float scale = tp.FindFloat("scale", 1.f);
    float gamma = tp.FindFloat("gamma", 1.f);
    return new ImageTexture<RGBSpectrum<rho_d>, Spectrum<rho_d>>(map, tp.FindString("filename"),
        trilerp, maxAniso, wrapMode, scale, gamma);
}


