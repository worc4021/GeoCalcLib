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

    A = mxGetCell(retVal, 0);
    b = mxGetCell(retVal, 1);

    if (nlhs==2)
        plhs[1] = b;

    plhs[0] = A;

}