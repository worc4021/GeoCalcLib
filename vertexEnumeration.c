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

    mwSize ndim=2, dims[]={1, 2}, nsubs=2, subs[2];
    mwIndex index;
    subs[0] = 0;
    subs[1] = 0;
    index = mxCalcSingleSubscript(retVal, nsubs, subs);

    V = mxGetCell(retVal, index);
    if (NULL == V)
        mexPrintf("V is null\n");

    subs[1] = 1;
    index = mxCalcSingleSubscript(retVal, nsubs, subs);

    type = mxGetCell(retVal,index);

    if (NULL == type)
        mexPrintf("type is null\n");

    if (nlhs==2)
        plhs[1] = mxDuplicateArray(type);

    plhs[0] = mxDuplicateArray(V);

    mxDestroyArray(retVal);


}