
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

#ifndef PBRT_CORE_IMAGEIO_H
#define PBRT_CORE_IMAGEIO_H

// core/imageio.h*
#include "pbrt.h"

// ImageIO Declarations
template<typename D>
RGBSpectrum<D> *ReadImage(const string &name, int *xSize, int *ySize);
inline void WriteImage(const string &name, float *pixels, float *alpha,
    int XRes, int YRes, int totalXRes, int totalYRes, int xOffset,
    int yOffset);

#include "imageio.hpp"
#endif // PBRT_CORE_IMAGEIO_H
