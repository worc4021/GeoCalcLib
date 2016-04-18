/*	File: mainFunctions.h
	
	Program: Geometric calculation library
	
	Description: The geometric calculation library provides a MATLAB callable
	functions to compute the vertex enumeration, facet enumeration and projection
	for convex polyhedra. In addition it provides functions to compute minimal
	representations of convex polyhedra. The geometric library uses the LRS library,
	see http://cgm.cs.mcgill.ca/~avis/C/lrs.html for details.

    Copyright (C) 2015  Rainer Manuel Schaich (rainer.schaich@eng.ox.ac.uk)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
*/

#include "mex.h"
#include "matrix.h"

mxArray *vertexEnumeration(const mxArray *A, const mxArray *b);
mxArray *facetEnumeration(const mxArray *V, const mxArray *type);
mxArray *vertexReduction(const mxArray *V, const mxArray *type);
mxArray *ineqReduction(const mxArray *A, const mxArray *b);
mxArray *projection(const mxArray *A, const mxArray *b, const mxArray *dim);
mxArray *lrsCall(const mxArray *dat);