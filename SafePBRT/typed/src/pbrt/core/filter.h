
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

#ifndef PBRT_CORE_FILTER_H
#define PBRT_CORE_FILTER_H

// core/filter.h*
#include "pbrt.h"

// Filter Declarations
class Filter {
public:
    // Filter Interface
    virtual ~Filter();
    Filter(mfloat<length_d> xw, mfloat<length_d> yw)
        : xWidth(xw), yWidth(yw), invXWidth(1.f/xw), invYWidth(1.f/yw) {
    }
    virtual float Evaluate(const mfloat<length_d> &x, const mfloat<length_d> &y) const = 0;

    // Filter Public Data
    const mfloat<length_d> xWidth, yWidth;
    const mfloat<invlength_d> invXWidth, invYWidth;
};



#endif // PBRT_CORE_FILTER_H
