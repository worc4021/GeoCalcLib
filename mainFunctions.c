#include "mainFunctions.h"
#include "translation_functions.h"


mxArray *vertexEnumeration(const mxArray *A, const mxArray *b)
{
	mxArray *helper, *retVal;
	struct GMPmat *myHelp, *calcMat;

	helper = VertConcat( b , A );
	myHelp = GMPmat_fromMXArray( helper );
	
	mxDestroyArray(helper);
	

	GMPmat_invertSignForFacetEnumeration(myHelp);
	calcMat = H2V(myHelp);
	GMPmat_invertSignForFacetEnumeration(calcMat);

	helper = MXArray_fromGMPmat(calcMat);
	GMPmat_destroy(calcMat);

	retVal = VertBreakdown(helper);
	mxDestroyArray(helper);
	return retVal;
}

mxArray *facetEnumeration(const mxArray *V, const mxArray *type)
{
	mxArray *helper, *retVal;
	struct GMPmat *myHelp;

	helper = VertConcat( type , V );
	myHelp = GMPmat_fromMXArray( helper );

	mxDestroyArray(helper);

	myHelp = reducevertices(myHelp);
	myHelp = V2H(myHelp);
	myHelp = reducemat(myHelp);
	GMPmat_invertSignForFacetEnumeration(myHelp);

	helper = MXArray_fromGMPmat(myHelp);
	GMPmat_destroy(myHelp);

	retVal = VertBreakdown(helper);
	mxDestroyArray(helper);
	return retVal;
}

mxArray *vertexReduction(const mxArray *V, const mxArray *type)
{
	mxArray *helper, *retVal;
	struct GMPmat *myHelp;

	helper = VertConcat( type , V );
	myHelp = GMPmat_fromMXArray( helper );

	mxDestroyArray(helper);

	GMPmat_invertSignForFacetEnumeration(myHelp);
	myHelp = reducevertices(myHelp);
	GMPmat_invertSignForFacetEnumeration(myHelp);

	helper = MXArray_fromGMPmat(myHelp);
	GMPmat_destroy(myHelp);

	retVal = VertBreakdown(helper);
	mxDestroyArray(helper);
	return retVal;
}

mxArray *ineqReduction(const mxArray *A, const mxArray *b)
{
	mxArray *helper, *retVal;
	struct GMPmat *myHelp;

	helper = VertConcat( b , A );
	myHelp = GMPmat_fromMXArray( helper );

	mxDestroyArray(helper);

	GMPmat_invertSignForFacetEnumeration(myHelp);
	myHelp = reducemat(myHelp);
	GMPmat_invertSignForFacetEnumeration(myHelp);

	helper = MXArray_fromGMPmat(myHelp);
	GMPmat_destroy(myHelp);

	retVal = VertBreakdown(helper);
	mxDestroyArray(helper);
	return retVal;

}

mxArray *projection(const mxArray *A, const mxArray *b, const mxArray *dim)
{
	mxArray *helper, *retVal;
	struct GMPmat *myHelp;

	size_t d = MXArray_to_integer(dim);

	helper = VertConcat( b , A );
	myHelp = GMPmat_fromMXArray( helper );

	mxDestroyArray(helper);

	GMPmat_invertSignForFacetEnumeration(myHelp);
	myHelp = H2V(myHelp);
	myHelp = GMPmat_dropCols(myHelp, d);
	myHelp = reducevertices(myHelp);
	myHelp = V2H(myHelp);
	myHelp = reducemat(myHelp);
	GMPmat_invertSignForFacetEnumeration(myHelp);

	helper = MXArray_fromGMPmat(myHelp);
	GMPmat_destroy(myHelp);

	retVal = VertBreakdown(helper);
	mxDestroyArray(helper);
	return retVal;

}
