
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

#ifndef PBRT_CORE_PARAMSET_H
#define PBRT_CORE_PARAMSET_H

// core/paramset.h*
#include "m_float.h"
#include "m_units.h"
#include "m_spaces.h"
#include "pbrt.h"
#include "geometry.h"
#include "texture.h"
#include "spectrum.h"
 #if (_MSC_VER >= 1400)
 #include <stdio.h>
 #define snprintf _snprintf
 #endif
#include <map>
using std::map;

template <typename T> struct ParamSetItem : public ReferenceCounted {
    // ParamSetItem Public Methods
    ParamSetItem(const string &name, const T *val, int nItems = 1);
    ~ParamSetItem() {
        delete[] data;
    }

    // ParamSetItem Data
    string name;
    int nItems;
    T *data;
    mutable bool lookedUp;
};


template<typename T, typename U>
T FindOneCast(const std::vector<U> &collection, const string &name, const T& d)
{
    for (uint32_t i = 0; i < collection.size(); ++i)
    {
        if (collection[i]->name == name && collection[i]->nItems == 1) 
        {
            collection[i]->lookedUp = true;
            return *(reinterpret_cast<T*>(collection[i]->data));
        }
    }
    return d;
}

template<typename T, typename U>
const T* FindCast(const std::vector<U> &collection, const string &name, int *nItems)
{
    for (uint32_t i = 0; i < collection.size(); ++i)
    {
        if (collection[i]->name == name) {
            *nItems = collection[i]->nItems;
            collection[i]->lookedUp = true;
            return reinterpret_cast<T*>(collection[i]->data);
        }
    }
    return NULL;
}

// ParamSet Declarations
class ParamSet {
public:
    // ParamSet Public Methods
    ParamSet() { }
    void AddFloat(const string &, const float *, int nItems = 1);
    void AddInt(const string &, const int *, int nItems);
    void AddBool(const string &, const bool *, int nItems);
    void AddPoint(const string &, const tuple3 *, int nItems);
    void AddVector(const string &, const tuple3 *, int nItems);
    void AddNormal(const string &, const tuple3 *, int nItems);
    void AddString(const string &, const string *, int nItems);
    void AddTexture(const string &, const string &);
    void AddRGBSpectrum(const string &, const float *, int nItems);
    void AddXYZSpectrum(const string &, const float *, int nItems);
    void AddBlackbodySpectrum(const string &, const float *, int nItems);
    void AddSampledSpectrumFiles(const string &, const char **, int nItems);
    void AddSampledSpectrum(const string &, const float *, int nItems);
    bool EraseInt(const string &);
    bool EraseBool(const string &);
    bool EraseFloat(const string &);
    bool ErasePoint(const string &);
    bool EraseVector(const string &);
    bool EraseNormal(const string &);
    bool EraseSpectrum(const string &);
    bool EraseString(const string &);
    bool EraseTexture(const string &);
    float FindOneFloat(const string &, float d) const;
    template<typename D>
    mfloat<D> FindOneFloat(const string &name, const mfloat<D> &d) const
    {
        return FindOneCast(floats, name, d);
    }
    int FindOneInt(const string &, int d) const;
    bool FindOneBool(const string &, bool d) const;
    template<typename S>
    Point<S> FindOnePoint(const string &name, const Point<S> &d) const
    {
        return FindOneCast(points, name, d);
    }
    template<typename S>
    Vector<S> FindOneVector(const string &name, const Vector<S> &d) const
    {
        return FindOneCast(vectors, name, d);
    }
    template<typename S>
    Normal<S> FindOneNormal(const string &name, const Normal<S> &d) const
    {
        return FindOneCast(normals, name, d);
    }
    template<typename D>
    Spectrum<D> FindOneSpectrum(const string &name,
                             const Spectrum<D> &d) const
    {
        return FindOneCast(spectra, name, d);
    }
    string FindOneString(const string &, const string &d) const;
    string FindTexture(const string &) const;
    const float *FindFloat(const string &, int *nItems) const;
    template<typename D>
    const mfloat<D> *FindFloat(const string &name, int *nItems) const
    {
        for (uint32_t i = 0; i < floats.size(); ++i)
        {
            if (floats[i]->name == name) {
                *nItems = floats[i]->nItems;
                floats[i]->lookedUp = true;
                return reinterpret_cast<mfloat<D>*>(floats[i]->data);
            }
        }
        return NULL;
    }
    const int *FindInt(const string &, int *nItems) const;
    const bool *FindBool(const string &, int *nItems) const;
    template<typename S>
    const Point<S> *FindPoint(const string &name, int *nItems) const
    {
        return FindCast<Point<S> >(points, name, nItems);
    }
    template<typename S>
    const Vector<S> *FindVector(const string &name, int *nItems) const
    {
        return FindCast<Vector<S> >(vectors, name, nItems);
    }
    template<typename S>
    const Normal<S> *FindNormal(const string &name, int *nItems) const
    {
        return FindCast<Normal<S> >(normals, name, nItems);
    }
    template<typename D>
    const Spectrum<D> *FindSpectrum(const string &name, int *nItems) const
    {
        return FindCast<Spectrum<D> >(spectra, name, nItems);
    }
    const string *FindString(const string &, int *nItems) const;
    void ReportUnused() const;
    void Clear();
    string ToString() const;
private:
    // ParamSet Private Data
    std::vector<Reference<ParamSetItem<bool> > > bools;
    std::vector<Reference<ParamSetItem<int> > > ints;
    std::vector<Reference<ParamSetItem<float> > > floats;
    std::vector<Reference<ParamSetItem<tuple3> > > points;
    std::vector<Reference<ParamSetItem<tuple3> > > vectors;
    std::vector<Reference<ParamSetItem<tuple3> > > normals;
    std::vector<Reference<ParamSetItem<tuple3> > > spectra;
    std::vector<Reference<ParamSetItem<string> > > strings;
    std::vector<Reference<ParamSetItem<string> > > textures;
};





// ParamSetItem Methods
template <typename T>
ParamSetItem<T>::ParamSetItem(const string &n, const T *v, int ni) {
    name = n;
    nItems = ni;
    data = new T[nItems];
    for (int i = 0; i < nItems; ++i) data[i] = v[i];
    lookedUp = false;
}

// TextureParams Declarations
class TextureParams {
public:
    // TextureParams Public Methods
    TextureParams(const ParamSet &geomp, const ParamSet &matp,
                  map<string, Reference<Texture<mfloat<unit_d> > > > &ft,
                  map<string, Reference<Texture<Spectrum<rho_d> > > > &st)
        : floatTextures(ft), spectrumTextures(st),
          geomParams(geomp), materialParams(matp) {
    }
    template<typename D>
    Reference<Texture<Spectrum<D> > > GetSpectrumTexture(const string &n,
            const Spectrum<D> &def) const;
    template<typename D>
    Reference<Texture<mfloat<D> > > GetFloatTexture(const string &name,
        const mfloat<D> & def) const;
    Reference<Texture<float> > GetFloatTexture(const string &name,
            float def) const;

    float FindFloat(const string &n, float d) const {
        return geomParams.FindOneFloat(n, materialParams.FindOneFloat(n, d));
    }
    string FindString(const string &n, const string &d = "") const {
           return geomParams.FindOneString(n, materialParams.FindOneString(n, d));
    }
    int FindInt(const string &n, int d) const {
           return geomParams.FindOneInt(n, materialParams.FindOneInt(n, d));
    }
    bool FindBool(const string &n, bool d) const {
           return geomParams.FindOneBool(n, materialParams.FindOneBool(n, d));
    }
    template<typename S>
    Point<S> FindPoint(const string &n, const Point<S> &d) const {
           return geomParams.FindOnePoint(n, materialParams.FindOnePoint(n, d));
    }
    template<typename S>
    Vector<S> FindVector(const string &n, const Vector<S> &d) const {
           return geomParams.FindOneVector(n, materialParams.FindOneVector(n, d));
    }
    template<typename S>
    Normal<S> FindNormal(const string &n, const Normal<S> &d) const {
           return geomParams.FindOneNormal(n, materialParams.FindOneNormal(n, d));
    }
    template<typename D>
    Spectrum<D> FindSpectrum(const string &n, const Spectrum<D> &d) const {
           return geomParams.FindOneSpectrum(n, materialParams.FindOneSpectrum(n, d));
    }
    void ReportUnused() const {
        geomParams.ReportUnused();
        materialParams.ReportUnused();
    }
    const ParamSet &GetGeomParams() const { return geomParams; }
    const ParamSet &GetMaterialParams() const { return materialParams; }
private:
    // TextureParams Private Data
    map<string, Reference<Texture<mfloat<unit_d> > > > &floatTextures;
    map<string, Reference<Texture<Spectrum<rho_d> > > > &spectrumTextures;
    const ParamSet &geomParams, &materialParams;
};

#include "textures/constant.h"

template<typename D>
Reference<Texture<Spectrum<D> > > TextureParams::GetSpectrumTexture(const string &n,
                                                     const Spectrum<D> &def) const
{
    string name = geomParams.FindTexture(n);
    if (name == "") name = materialParams.FindTexture(n);
    if (name != "") {
        if (spectrumTextures.find(name) != spectrumTextures.end())
            return *reinterpret_cast<Reference<Texture<Spectrum<D> > > *>(&spectrumTextures[name]);
        else
            Error("Couldn't find spectrum texture \"%s\"", n.c_str());
    }
    Spectrum<D> val = geomParams.FindOneSpectrum<D>(n,
        materialParams.FindOneSpectrum(n, def));
    return new ConstantTexture<Spectrum<D> >(val);
}

template<typename D>
Reference<Texture<mfloat<D> > > TextureParams::GetFloatTexture(const string &n,
                                                          const mfloat<D> &def) const 
{
    string name = geomParams.FindTexture(n);
    if (name == "") name = materialParams.FindTexture(n);
    if (name != "") {
        if (floatTextures.find(name) != floatTextures.end())
            return *reinterpret_cast<Reference<Texture<mfloat<D> > > *>(&floatTextures[name]);
        else
            Error("Couldn't find float texture named \"%s\"", n.c_str());
    }
    mfloat<D> val = geomParams.FindOneFloat<D>(n,
        materialParams.FindOneFloat<D>(n, def));
    return new ConstantTexture<mfloat<D> >(val);
}
#endif // PBRT_CORE_PARAMSET_H
