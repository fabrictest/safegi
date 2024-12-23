
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

#ifndef PBRT_TEXTURES_IMAGEMAP_H
#define PBRT_TEXTURES_IMAGEMAP_H

// textures/imagemap.h*
#include "pbrt.h"
#include "texture.h"
#include "mipmap.h"
#include "paramset.h"
#include <map>
#include "imageio.h"

// TexInfo Declarations
struct TexInfo {
    TexInfo(const string &f, bool dt, float ma, ImageWrap wm, float sc, float ga)
        : filename(f), doTrilinear(dt), maxAniso(ma), wrapMode(wm), scale(sc), gamma(ga) { }
    string filename;
    bool doTrilinear;
    float maxAniso;
    ImageWrap wrapMode;
    float scale, gamma;
    bool operator<(const TexInfo &t2) const {
        if (filename != t2.filename) return filename < t2.filename;
        if (doTrilinear != t2.doTrilinear) return doTrilinear < t2.doTrilinear;
        if (maxAniso != t2.maxAniso) return maxAniso < t2.maxAniso;
        if (scale != t2.scale) return scale < t2.scale;
        if (gamma != t2.gamma) return gamma < t2.gamma;
        return wrapMode < t2.wrapMode;
    }
};



// ImageTexture Declarations
template <typename Tmemory, typename Treturn>
    class ImageTexture : public Texture<Treturn> {
public:
    // ImageTexture Public Methods
    ImageTexture(TextureMapping2D *m, const string &filename, bool doTri,
                 float maxAniso, ImageWrap wm, float scale, float gamma);
    Treturn Evaluate(const DifferentialGeometry &) const;
    ~ImageTexture();
    static void ClearCache() {
        typename std::map<TexInfo, MIPMap<Tmemory> *>::iterator iter;
        iter = textures.begin();
        while (iter != textures.end()) {
            delete iter->second;
            ++iter;
        }
        textures.erase(textures.begin(), textures.end());
    }
private:
    // ImageTexture Private Methods
    static MIPMap<Tmemory> *GetTexture(const string &filename,
        bool doTrilinear, float maxAniso, ImageWrap wm, float scale, float gamma);

    template<typename D>
    static void convertIn(const RGBSpectrum<D> &from, RGBSpectrum<D> *to,
                          float scale, float gamma) 
    {
        RGBSpectrum<D> s = scale * from;
        float rgb[3];
        s.ToRGB(rgb);
        rgb[0] = pow(rgb[0], gamma);
        rgb[1] = pow(rgb[1], gamma);
        rgb[2] = pow(rgb[2], gamma);
        *to = RGBSpectrum<D>::FromRGB(rgb);
    }
    template<typename D>
    static void convertIn(const RGBSpectrum<D> &from, mfloat<D> *to,
                          float scale, float gamma) {
        *to = mfloat<D>(pow(scale * from.y(), gamma));
    }
    template<typename D>
    static void convertOut(const RGBSpectrum<D> &from, Spectrum<D> *to) {
        float rgb[3];
        from.ToRGB(rgb);
        *to = Spectrum<D>::FromRGB(rgb);
    }
    template<typename D>
    static void convertOut(mfloat<D> from, mfloat<D> *to) {
        *to = from;
    }

    // ImageTexture Private Data
    MIPMap<Tmemory> *mipmap;
    TextureMapping2D *mapping;
    static std::map<TexInfo, MIPMap<Tmemory> *> textures;
};

template <typename Tmemory, typename Treturn>
ImageTexture<Tmemory, Treturn>::ImageTexture(TextureMapping2D *m,
        const string &filename, bool doTrilinear, float maxAniso,
        ImageWrap wrapMode, float scale, float gamma) {
    mapping = m;
    mipmap = GetTexture(filename, doTrilinear, maxAniso,
                        wrapMode, scale, gamma);
}


template <typename Tmemory, typename Treturn>
    ImageTexture<Tmemory, Treturn>::~ImageTexture() {
    delete mapping;
}


template <typename Tmemory, typename Treturn> MIPMap<Tmemory> *
ImageTexture<Tmemory, Treturn>::GetTexture(const string &filename,
        bool doTrilinear, float maxAniso, ImageWrap wrap,
        float scale, float gamma) {
    // Look for texture in texture cache
    TexInfo texInfo(filename, doTrilinear, maxAniso, wrap, scale, gamma);
    if (textures.find(texInfo) != textures.end())
        return textures[texInfo];
    int width, height;
    RGBSpectrum<unit_d> *texels = ReadImage<unit_d>(filename, &width, &height);
    MIPMap<Tmemory> *ret = NULL;
    if (texels) {
        // Convert texels to type _Tmemory_ and create _MIPMap_
        Tmemory *convertedTexels = new Tmemory[width*height];
        for (int i = 0; i < width*height; ++i)
            convertIn(texels[i], &convertedTexels[i], scale, gamma);
        ret = new MIPMap<Tmemory>(width, height, convertedTexels, doTrilinear,
                                  maxAniso, wrap);
        delete[] texels;
        delete[] convertedTexels;
    }
    else {
        // Create one-valued _MIPMap_
        Tmemory *oneVal = new Tmemory[1];
        oneVal[0] = texval<Tmemory>::convert_float(powf(scale, gamma));
        ret = new MIPMap<Tmemory>(1, 1, oneVal);
        delete[] oneVal;
    }
    textures[texInfo] = ret;
    PBRT_LOADED_IMAGE_MAP(const_cast<char *>(filename.c_str()), width, height, sizeof(Tmemory), ret);
    return ret;
}


template <typename Tmemory, typename Treturn>
    std::map<TexInfo,
             MIPMap<Tmemory> *> ImageTexture<Tmemory, Treturn>::textures;
template <typename Tmemory, typename Treturn> Treturn
ImageTexture<Tmemory,
             Treturn>::Evaluate(const DifferentialGeometry &dg) const {
    float s, t, dsdx, dtdx, dsdy, dtdy;
    mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
    Tmemory mem = mipmap->Lookup(s, t, dsdx, dtdx, dsdy, dtdy);
    Treturn ret;
    convertOut(mem, &ret);
    return ret;
}

ImageTexture<mfloat<unit_d>, mfloat<unit_d>> *CreateImageFloatTexture(const Transform<texture_s, world_s> &tex2world,
        const TextureParams &tp);
ImageTexture<RGBSpectrum<rho_d>, Spectrum<rho_d>> *CreateImageSpectrumTexture(const Transform<texture_s, world_s> &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_IMAGEMAP_H
