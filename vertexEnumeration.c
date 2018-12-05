#include "mex.h"
#include "mainFunctions.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs<2)
        mexErrMsgTxt("Vertex enumeration requires A,b to define {A*x<=b}.");

    const mxArray *A;
    const mxArray *b;
    mxArray *retVal;
    mxArray *V;
    mxArray *type;

    A = prhs[0];
    b = prhs[1];

    retVal = vertexEnumeration(A,b);

    V = mxGetCell(retVal, 0);
    type = mxGetCell(retVal,1);

    if (nlhs==2)
        plhs[1] = type;

    plhs[0] = V;

}