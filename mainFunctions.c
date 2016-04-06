/*	File: mainFunctions.c
	
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

#include "mainFunctions.h"
#include "translation_functions.h"


mxArray *vertexEnumeration(const mxArray *A, const mxArray *b)
{
	if( mxGetM(A) != mxGetM(b) )
		mexErrMsgIdAndTxt("GeoCalcLib:Dimensions:vertex:Enumeration","Matrices A and b must have same number of rows.\n" );

	if( mxIsEmpty(A) || mxIsEmpty(b) )
		mexErrMsgIdAndTxt("GeoCalcLib:Dimensions:vertex:Enumeration","Matrices must not be empty.\n" );

	mxArray *helper, *retVal;
	struct GMPmat *myHelp, *calcMat;

	helper = VertConcat( b , A );
	myHelp = GMPmat_fromMXArray( helper );
	
	mxDestroyArray(helper);

	GMPmat_invertSignForFacetEnumeration(myHelp);
	calcMat = H2V(myHelp);

	if (!calcMat->empty) {
		helper = MXArray_fromGMPmat(calcMat);
		GMPmat_destroy(calcMat);

		retVal = VertBreakdown(helper);
		mxDestroyArray(helper);
	} else {
		GMPmat_destroy(calcMat);
		retVal = createEmptyCell();
	}

	return retVal;
}

mxArray *facetEnumeration(const mxArray *V, const mxArray *type)
{
	if( mxGetM(V) != mxGetM(type)) 
		mexErrMsgIdAndTxt("GeoCalcLib:Dimensions:facet:Enumeration","Matrices V and type must have same number of rows.\n" );

	if( mxIsEmpty(V) || mxIsEmpty(type) )
		mexErrMsgIdAndTxt("GeoCalcLib:Dimensions:vertex:Enumeration","Matrices must not be empty.\n" );

	mxArray *helper, *retVal;
	struct GMPmat *myHelp;

	helper = VertConcat( type , V );
	myHelp = GMPmat_fromMXArray( helper );

	mxDestroyArray(helper);

	myHelp = reducevertices(myHelp);
	myHelp = GMPmat_lift(myHelp);
	myHelp = H2V(myHelp);
	myHelp = GMPmat_unlift(myHelp);
	if (!myHelp->empty) {

		myHelp = reducemat(myHelp);
		GMPmat_invertSignForFacetEnumeration(myHelp);

		helper = MXArray_fromGMPmat(myHelp);
		GMPmat_destroy(myHelp);

		retVal = VertBreakdown(helper);
		mxDestroyArray(helper);
	} else {
		GMPmat_destroy(myHelp);
		retVal = createEmptyCell();
	}
	return retVal;
}

mxArray *vertexReduction(const mxArray *V, const mxArray *type)
{
	if( mxGetM(V) != mxGetM(type)) 
		mexErrMsgIdAndTxt("GeoCalcLib:Dimensions:vertex:Reduction","Matrices V and type must have same number of rows.\n" );

	if( mxIsEmpty(V) || mxIsEmpty(type) )
		mexErrMsgIdAndTxt("GeoCalcLib:Dimensions:vertex:Enumeration","Matrices must not be empty.\n" );

	mxArray *helper, *retVal;
	struct GMPmat *myHelp;

	helper = VertConcat( type , V );
	myHelp = GMPmat_fromMXArray( helper );

	mxDestroyArray(helper);

	myHelp = reducevertices(myHelp);
	if (!myHelp->empty) {
		helper = MXArray_fromGMPmat(myHelp);
		GMPmat_destroy(myHelp);

		retVal = VertBreakdown(helper);
		mxDestroyArray(helper);

	} else {
		GMPmat_destroy(myHelp);
		retVal = createEmptyCell();
	}
	return retVal;
}

mxArray *ineqReduction(const mxArray *A, const mxArray *b)
{
	if( mxGetM(A) != mxGetM(b)) 
		mexErrMsgIdAndTxt("GeoCalcLib:Dimensions:ineq:Reduction","Matrices A and b must have same number of rows.\n" );

	if( mxIsEmpty(A) || mxIsEmpty(b) )
		mexErrMsgIdAndTxt("GeoCalcLib:Dimensions:vertex:Enumeration","Matrices must not be empty.\n" );

	mxArray *helper, *retVal;
	struct GMPmat *myHelp;

	helper = VertConcat( b , A );
	myHelp = GMPmat_fromMXArray( helper );

	mxDestroyArray(helper);

	GMPmat_invertSignForFacetEnumeration(myHelp);
	myHelp = reducemat(myHelp);
	if (!myHelp->empty) {
	GMPmat_invertSignForFacetEnumeration(myHelp);

	helper = MXArray_fromGMPmat(myHelp);
	GMPmat_destroy(myHelp);

	retVal = VertBreakdown(helper);
	mxDestroyArray(helper);
	} else {
		GMPmat_destroy(myHelp);
		retVal = createEmptyCell();
	}
	return retVal;

}

mxArray *projection(const mxArray *A, const mxArray *b, const mxArray *dim)
{
	if( mxGetM(A) != mxGetM(b)) 
		mexErrMsgIdAndTxt("GeoCalcLib:Dimensions:projection","Matrices A and b must have same number of rows.\n" );

	if( mxIsEmpty(A) || mxIsEmpty(b) )
		mexErrMsgIdAndTxt("GeoCalcLib:Dimensions:vertex:Enumeration","Matrices must not be empty.\n" );

	mxArray *helper, *retVal;
	struct GMPmat *myHelp;

	char BREAK = 0;

	size_t d = MXArray_to_integer(dim);

	if ( mxGetN(A) <= d)
		mexErrMsgIdAndTxt("GeoCalcLib:Codimension:projection","dim must smaller than number of columns of A.\n");

	helper = VertConcat( b , A );
	myHelp = GMPmat_fromMXArray( helper );

	mxDestroyArray(helper);

	GMPmat_invertSignForFacetEnumeration(myHelp);
	myHelp = H2V(myHelp);
	if (!myHelp->empty) 
	{
		myHelp = GMPmat_dropCols(myHelp, d);
		myHelp = reducevertices(myHelp);
		
		if (!myHelp->empty) 
		{
			myHelp = V2H(myHelp);
			
			if (!myHelp->empty) 
			{
				myHelp = reducemat(myHelp);
				GMPmat_invertSignForFacetEnumeration(myHelp);
				helper = MXArray_fromGMPmat(myHelp);
				GMPmat_destroy(myHelp);

				retVal = VertBreakdown(helper);
				mxDestroyArray(helper);
			}
			else
			{
			BREAK = 1;
			}
		} 
		else
		{
			BREAK = 1;
		}
	} 
	else
	{ 
		BREAK = 1;
	}
	
	if (BREAK) 
	{
		GMPmat_destroy(myHelp);
		retVal = createEmptyCell();
	}
	
	return retVal;

}
