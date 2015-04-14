#include "mex.h"
#include "matrix.h"

mxArray *vertexEnumeration(const mxArray *A, const mxArray *b);
mxArray *facetEnumeration(const mxArray *V, const mxArray *type);
mxArray *vertexReduction(const mxArray *V, const mxArray *type);
mxArray *ineqReduction(const mxArray *A, const mxArray *b);
mxArray *projection(const mxArray *A, const mxArray *b, const mxArray *dim);