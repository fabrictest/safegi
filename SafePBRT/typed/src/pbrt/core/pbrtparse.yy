
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

%{
#include "api.h"
#include "pbrt.h"
#include "paramset.h"
#include <stdarg.h>

#ifdef WIN32
#pragma warning(disable:4065)
#pragma warning(disable:4996)
#pragma warning(disable:4018)
#endif // WIN32

extern int yylex();
int line_num = 0;
string current_file;

#define YYMAXDEPTH 100000000

void yyerror(const char *str) {
    Severe("Parsing error: %s", str);
}



struct ParamArray {
    ParamArray() {
        isString = false;
        element_size = allocated = nelems = 0;
        array = NULL;
    }
    bool isString;
    int element_size;
    int allocated;
    int nelems;
    void *array;
};



struct ParamListItem {
    ParamListItem(const char *t, ParamArray *array) {
        arg = array->array;
        name = t;
        size = array->nelems;
        isString = array->isString;
        array->allocated = 0;
        array->nelems = 0;
        array->array = NULL;
    }
    const char *name;
    void *arg;
    int size;
    bool isString;
};



static std::vector<ParamListItem> cur_paramlist;

static ParamArray *cur_array = NULL;

static void AddArrayElement(void *elem) {
    if (cur_array->nelems >= cur_array->allocated) {
        cur_array->allocated = 2*cur_array->allocated + 1;
        cur_array->array = realloc(cur_array->array,
            cur_array->allocated*cur_array->element_size);
    }
    char *next = ((char *)cur_array->array) + cur_array->nelems * cur_array->element_size;
    Assert(cur_array->element_size == 4 || cur_array->element_size == 8);
    if (cur_array->element_size == 4)
        *((uint32_t *)next) = *((uint32_t *)elem);
    else
        *((uint64_t *)next) = *((uint64_t *)elem);
    cur_array->nelems++;
}



static void ArrayFree(ParamArray *ra) {
    if (ra->isString && ra->array)
        for (int i = 0; i < ra->nelems; ++i) free(((char **)ra->array)[i]);
    free(ra->array);
    delete ra;
}



static void FreeArgs() {
    for (uint32_t i = 0; i < cur_paramlist.size(); ++i)
        free((char *)cur_paramlist[i].arg);
    cur_paramlist.erase(cur_paramlist.begin(), cur_paramlist.end());
}



static bool VerifyArrayLength(ParamArray *arr, int required,
    const char *command) {
    if (arr->nelems != required) {
        Error("\"%s\" requires a %d element array! (%d found)",
                    command, required, arr->nelems);
        return false;
    }
    return true;
}


enum { PARAM_TYPE_INT, PARAM_TYPE_BOOL, PARAM_TYPE_FLOAT, PARAM_TYPE_POINT,
    PARAM_TYPE_VECTOR, PARAM_TYPE_NORMAL, PARAM_TYPE_RGB, PARAM_TYPE_XYZ,
    PARAM_TYPE_BLACKBODY, PARAM_TYPE_SPECTRUM,
    PARAM_TYPE_STRING, PARAM_TYPE_TEXTURE };
static const char *paramTypeToName(int type);
static void InitParamSet(ParamSet &ps, SpectrumType);
static bool lookupType(const char *name, int *type, string &sname);
#define YYPRINT(file, type, value)  { \
    if ((type) == ID || (type) == STRING) \
        fprintf ((file), " %s", (value).string); \
    else if ((type) == NUM) \
        fprintf ((file), " %f", (value).num); \
}


%}

%union {
char string[1024];
float num;
ParamArray *ribarray;
}


%token <string> STRING ID
%token <num> NUM
%token LBRACK RBRACK

%token ACCELERATOR ACTIVETRANSFORM ALL AREALIGHTSOURCE ATTRIBUTEBEGIN
%token ATTRIBUTEEND CAMERA CONCATTRANSFORM COORDINATESYSTEM COORDSYSTRANSFORM
%token ENDTIME FILM IDENTITY LIGHTSOURCE LOOKAT MAKENAMEDMATERIAL MATERIAL
%token NAMEDMATERIAL OBJECTBEGIN OBJECTEND OBJECTINSTANCE PIXELFILTER RENDERER
%token REVERSEORIENTATION ROTATE SAMPLER SCALE SHAPE STARTTIME
%token SURFACEINTEGRATOR TEXTURE TRANSFORMBEGIN TRANSFORMEND TRANSFORMTIMES
%token TRANSFORM TRANSLATE VOLUME VOLUMEINTEGRATOR WORLDBEGIN WORLDEND

%token HIGH_PRECEDENCE

%type<ribarray> array num_array string_array
%%
start: pbrt_stmt_list
{
};



array_init: %prec HIGH_PRECEDENCE
{
    if (cur_array) Severe("MUH");
    cur_array = new ParamArray;
};



string_array_init: %prec HIGH_PRECEDENCE
{
    cur_array->element_size = sizeof(const char *);
    cur_array->isString = true;
};



num_array_init: %prec HIGH_PRECEDENCE
{
    cur_array->element_size = sizeof(float);
    cur_array->isString = false;
};



array: string_array
{
    $$ = $1;
}


| num_array
{
    $$ = $1;
};



string_array: array_init LBRACK string_list RBRACK
{
    $$ = cur_array;
    cur_array = NULL;
}


| single_element_string_array
{
    $$ = cur_array;
    cur_array = NULL;
};



single_element_string_array: array_init string_list_entry
{
};



string_list: string_list string_list_entry
{
}


| string_list_entry
{
};



string_list_entry: string_array_init STRING
{
    char *to_add = strdup($2);
    AddArrayElement(&to_add);
};



num_array: array_init LBRACK num_list RBRACK
{
    $$ = cur_array;
    cur_array = NULL;
}


| single_element_num_array
{
    $$ = cur_array;
    cur_array = NULL;
};



single_element_num_array: array_init num_list_entry
{
};



num_list: num_list num_list_entry
{
}


| num_list_entry
{
};



num_list_entry: num_array_init NUM
{
    float to_add = $2;
    AddArrayElement(&to_add);
};



paramlist: paramlist_init paramlist_contents
{
};



paramlist_init: %prec HIGH_PRECEDENCE
{
    for (uint32_t i = 0; i < cur_paramlist.size(); ++i) {
        if (cur_paramlist[i].isString) {
            for (uint32_t j = 0; j < cur_paramlist[i].size; ++j)
                free(((char **)cur_paramlist[i].arg)[j]);
        }
    }
    cur_paramlist.erase(cur_paramlist.begin(), cur_paramlist.end());
};



paramlist_contents: paramlist_entry paramlist_contents
{
}


|
{
};



paramlist_entry: STRING array
{
    cur_paramlist.push_back(ParamListItem($1, $2));
    ArrayFree($2);
};



pbrt_stmt_list: pbrt_stmt_list pbrt_stmt
{
}


| pbrt_stmt
{
};



pbrt_stmt: ACCELERATOR STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtAccelerator($2, params);
    FreeArgs();
}


| ACTIVETRANSFORM ALL
{
    pbrtActiveTransformAll();
}


| ACTIVETRANSFORM ENDTIME
{
    pbrtActiveTransformEndTime();
}


| ACTIVETRANSFORM STARTTIME
{
    pbrtActiveTransformStartTime();
}


| AREALIGHTSOURCE STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_ILLUMINANT);
    pbrtAreaLightSource($2, params);
    FreeArgs();
}


| ATTRIBUTEBEGIN
{
    pbrtAttributeBegin();
}


| ATTRIBUTEEND
{
    pbrtAttributeEnd();
}


| CAMERA STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtCamera($2, params);
    FreeArgs();
}


| CONCATTRANSFORM num_array
{
    if (VerifyArrayLength($2, 16, "ConcatTransform"))
        pbrtConcatTransform((float *) $2->array);
    ArrayFree($2);
}


| COORDINATESYSTEM STRING
{
    pbrtCoordinateSystem($2);
}


| COORDSYSTRANSFORM STRING
{
    pbrtCoordSysTransform($2);
}


| FILM STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtFilm($2, params);
    FreeArgs();
}


| IDENTITY
{
    pbrtIdentity();
}


| LIGHTSOURCE STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_ILLUMINANT);
    pbrtLightSource($2, params);
    FreeArgs();
}


| LOOKAT NUM NUM NUM NUM NUM NUM NUM NUM NUM
{
    pbrtLookAt($2, $3, $4, $5, $6, $7, $8, $9, $10);
}


| MAKENAMEDMATERIAL STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtMakeNamedMaterial($2, params);
    FreeArgs();
}


| MATERIAL STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtMaterial($2, params);
    FreeArgs();
}


| NAMEDMATERIAL STRING
{
    pbrtNamedMaterial($2);
}


| OBJECTBEGIN STRING
{
    pbrtObjectBegin($2);
}


| OBJECTEND
{
    pbrtObjectEnd();
}


| OBJECTINSTANCE STRING
{
    pbrtObjectInstance($2);
}


| PIXELFILTER STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtPixelFilter($2, params);
    FreeArgs();
}


| RENDERER STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtRenderer($2, params);
    FreeArgs();
}


| REVERSEORIENTATION
{
    pbrtReverseOrientation();
}


| ROTATE NUM NUM NUM NUM
{
    pbrtRotate($2, $3, $4, $5);
}


| SAMPLER STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtSampler($2, params);
    FreeArgs();
}


| SCALE NUM NUM NUM
{
    pbrtScale($2, $3, $4);
}


| SHAPE STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtShape($2, params);
    FreeArgs();
}


| SURFACEINTEGRATOR STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtSurfaceIntegrator($2, params);
    FreeArgs();
}


| TEXTURE STRING STRING STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtTexture($2, $3, $4, params);
    FreeArgs();
}


| TRANSFORMBEGIN
{
    pbrtTransformBegin();
}


| TRANSFORMEND
{
    pbrtTransformEnd();
}


| TRANSFORMTIMES NUM NUM
{
    pbrtTransformTimes($2, $3);
}


| TRANSFORM num_array
{
    if (VerifyArrayLength( $2, 16, "Transform" ))
        pbrtTransform( (float *) $2->array );
    ArrayFree($2);
}


| TRANSLATE NUM NUM NUM
{
    pbrtTranslate($2, $3, $4);
}


| VOLUMEINTEGRATOR STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtVolumeIntegrator($2, params);
    FreeArgs();
}


| VOLUME STRING paramlist
{
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtVolume($2, params);
    FreeArgs();
}


| WORLDBEGIN
{
    pbrtWorldBegin();
}


| WORLDEND
{
    pbrtWorldEnd();
};


%%
static const char *paramTypeToName(int type) {
    switch (type) {
    case PARAM_TYPE_INT: return "int";
    case PARAM_TYPE_BOOL: return "bool";
    case PARAM_TYPE_FLOAT: return "float";
    case PARAM_TYPE_POINT: return "point";
    case PARAM_TYPE_VECTOR: return "vector";
    case PARAM_TYPE_NORMAL: return "normal";
    case PARAM_TYPE_RGB: return "rgb/color";
    case PARAM_TYPE_XYZ: return "xyz";
    case PARAM_TYPE_BLACKBODY: return "blackbody";
    case PARAM_TYPE_SPECTRUM: return "spectrum";
    case PARAM_TYPE_STRING: return "string";
    case PARAM_TYPE_TEXTURE: return "texture";
    default: Severe("Error in paramTypeToName"); return NULL;
    }
}


static void InitParamSet(ParamSet &ps, SpectrumType type) {
    ps.Clear();
    for (uint32_t i = 0; i < cur_paramlist.size(); ++i) {
        int type;
        string name;
        if (lookupType(cur_paramlist[i].name, &type, name)) {
            if (type == PARAM_TYPE_TEXTURE || type == PARAM_TYPE_STRING ||
                type == PARAM_TYPE_BOOL) {
                if (!cur_paramlist[i].isString) {
                    Error("Expected string parameter value for parameter \"%s\" with type \"%s\". Ignoring.",
                          name.c_str(), paramTypeToName(type));
                    continue;
                }
            }
            else if (type != PARAM_TYPE_SPECTRUM) { /* spectrum can be either... */
                if (cur_paramlist[i].isString) {
                    Error("Expected numeric parameter value for parameter \"%s\" with type \"%s\".  Ignoring.",
                          name.c_str(), paramTypeToName(type));
                    continue;
                }
            }
            void *data = cur_paramlist[i].arg;
            int nItems = cur_paramlist[i].size;
            if (type == PARAM_TYPE_INT) {
                // parser doesn't handle ints, so convert from floats here....
                int nAlloc = nItems;
                int *idata = new int[nAlloc];
                float *fdata = (float *)cur_paramlist[i].arg;
                for (int j = 0; j < nAlloc; ++j)
                    idata[j] = int(fdata[j]);
                ps.AddInt(name, idata, nItems);
                delete[] idata;
            }
            else if (type == PARAM_TYPE_BOOL) {
                // strings -> bools
                int nAlloc = cur_paramlist[i].size;
                bool *bdata = new bool[nAlloc];
                for (int j = 0; j < nAlloc; ++j) {
                    string s(((const char **)data)[j]);
                    if (s == "true") bdata[j] = true;
                    else if (s == "false") bdata[j] = false;
                    else {
                        Warning("Value \"%s\" unknown for boolean parameter \"%s\"."
                            "Using \"false\".", s.c_str(), cur_paramlist[i].name);
                        bdata[j] = false;
                    }
                }
                ps.AddBool(name, bdata, nItems);
                delete[] bdata;
            }
            else if (type == PARAM_TYPE_FLOAT) {
                ps.AddFloat(name, (float *)data, nItems);
            } else if (type == PARAM_TYPE_POINT) {
                if ((nItems % 3) != 0)
                    Warning("Excess values given with point parameter \"%s\". "
                            "Ignoring last %d of them", cur_paramlist[i].name, nItems % 3);
                ps.AddPoint(name, (tuple3 *)data, nItems / 3);
            } else if (type == PARAM_TYPE_VECTOR) {
                if ((nItems % 3) != 0)
                    Warning("Excess values given with vector parameter \"%s\". "
                            "Ignoring last %d of them", cur_paramlist[i].name, nItems % 3);
                ps.AddVector(name, (tuple3 *)data, nItems / 3);
            } else if (type == PARAM_TYPE_NORMAL) {
                if ((nItems % 3) != 0)
                    Warning("Excess values given with normal parameter \"%s\". "
                            "Ignoring last %d of them", cur_paramlist[i].name, nItems % 3);
                ps.AddNormal(name, (tuple3 *)data, nItems / 3);
            } else if (type == PARAM_TYPE_RGB) {
                if ((nItems % 3) != 0)
                    Warning("Excess RGB values given with parameter \"%s\". "
                            "Ignoring last %d of them", cur_paramlist[i].name, nItems % 3);
                ps.AddRGBSpectrum(name, (float *)data, nItems);
            } else if (type == PARAM_TYPE_XYZ) {
                if ((nItems % 3) != 0)
                    Warning("Excess XYZ values given with parameter \"%s\". "
                            "Ignoring last %d of them", cur_paramlist[i].name, nItems % 3);
                ps.AddXYZSpectrum(name, (float *)data, nItems);
            } else if (type == PARAM_TYPE_BLACKBODY) {
                if ((nItems % 2) != 0)
                    Warning("Excess value given with blackbody parameter \"%s\". "
                            "Ignoring extra one.", cur_paramlist[i].name);
                ps.AddBlackbodySpectrum(name, (float *)data, nItems);
            } else if (type == PARAM_TYPE_SPECTRUM) {
                if (cur_paramlist[i].isString) {
                    ps.AddSampledSpectrumFiles(name, (const char **)data, nItems);
                }
                else {
                    if ((nItems % 2) != 0)
                        Warning("Non-even number of values given with sampled spectrum "
                                "parameter \"%s\". Ignoring extra.", cur_paramlist[i].name);
                    ps.AddSampledSpectrum(name, (float *)data, nItems);
                }
            } else if (type == PARAM_TYPE_STRING) {
                string *strings = new string[nItems];
                for (int j = 0; j < nItems; ++j)
                    strings[j] = string(((const char **)data)[j]);
                ps.AddString(name, strings, nItems);
                delete[] strings;
            }
            else if (type == PARAM_TYPE_TEXTURE) {
                if (nItems == 1) {
                    string val(*((const char **)data));
                    ps.AddTexture(name, val);
                }
                else
                    Error("Only one string allowed for \"texture\" parameter \"%s\"",
                        name.c_str());
            }
        }
        else
            Warning("Type of parameter \"%s\" is unknown",
                cur_paramlist[i].name);
    }
}


static bool lookupType(const char *name, int *type, string &sname) {
    Assert(name != NULL);
    *type = 0;
    const char *strp = name;
    while (*strp && isspace(*strp))
        ++strp;
    if (!*strp) {
        Error("Parameter \"%s\" doesn't have a type declaration?!", name);
        return false;
    }
#define TRY_DECODING_TYPE(name, mask) \
        if (strncmp(name, strp, strlen(name)) == 0) { \
            *type = mask; strp += strlen(name); \
        }
         TRY_DECODING_TYPE("float",     PARAM_TYPE_FLOAT)
    else TRY_DECODING_TYPE("integer",   PARAM_TYPE_INT)
    else TRY_DECODING_TYPE("bool",      PARAM_TYPE_BOOL)
    else TRY_DECODING_TYPE("point",     PARAM_TYPE_POINT)
    else TRY_DECODING_TYPE("vector",    PARAM_TYPE_VECTOR)
    else TRY_DECODING_TYPE("normal",    PARAM_TYPE_NORMAL)
    else TRY_DECODING_TYPE("string",    PARAM_TYPE_STRING)
    else TRY_DECODING_TYPE("texture",   PARAM_TYPE_TEXTURE)
    else TRY_DECODING_TYPE("color",     PARAM_TYPE_RGB)
    else TRY_DECODING_TYPE("rgb",       PARAM_TYPE_RGB)
    else TRY_DECODING_TYPE("xyz",       PARAM_TYPE_XYZ)
    else TRY_DECODING_TYPE("blackbody", PARAM_TYPE_BLACKBODY)
    else TRY_DECODING_TYPE("spectrum",  PARAM_TYPE_SPECTRUM)
    else {
        Error("Unable to decode type for name \"%s\"", name);
        return false;
    }
    while (*strp && isspace(*strp))
        ++strp;
    sname = string(strp);
    return true;
}


