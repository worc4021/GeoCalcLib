#include "mex.h"
#include "mainFunctions.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    const mxArray *V;
    mxArray *type;
    double *typePtr;
    if (nrhs!=2){
        type = mxCreateDoubleMatrix(mxGetM(V),1, mxREAL);
        typePtr = mxGetPr(type);
        for (mwSize i = 0; i != mxGetM(V); i++)
            typePtr[i] = 1.;
    } else {
        type = prhs[1];
    }

    V = prhs[0];

    mxArray *retVal;
    mxArray *A;
    mxArray *b;

    retVal = facetEnumeration(V,type);

    mwSize ndim=2, dims[]={1, 2}, nsubs=2, subs[2];
    mwIndex index;
    subs[0] = 0;
    subs[1] = 0;
    index = mxCalcSingleSubscript(retVal, nsubs, subs);

    A = mxGetCell(retVal, index);
    if (NULL == A)
        mexPrintf("A is null\n");

    subs[1] = 1;
    index = mxCalcSingleSubscript(retVal, nsubs, subs);
    b = mxGetCell(retVal, index);
    if (NULL == b)
        mexPrintf("b is null\n");

    if (nlhs==2)
        plhs[1] = b;

    plhs[0] = A;

}