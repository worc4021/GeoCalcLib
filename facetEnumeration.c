#include "mex.h"
#include "mainFunctions.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs<1)
        mexErrMsgTxt("Not enough inputs.");
        
    const mxArray *V;
    const mxArray *type;
    double *typePtr;
    int doClean = 0;

    V = prhs[0];

    if (nrhs<2){
        doClean = 1;
        type = mxCreateDoubleMatrix(mxGetM(V),1, mxREAL);
        typePtr = mxGetPr(type);
        for (mwSize i = 0; i != mxGetM(V); i++)
            typePtr[i] = 1.;
    } else {
        type = prhs[1];
    }

    
    mxArray *retVal;
    mxArray *A;
    mxArray *b;

    retVal = facetEnumeration(V,type);


    mwSize nsubs=2, subs[2];
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
    plhs[0] = mxDuplicateArray(A);

    if (nlhs==2)
        plhs[1] = mxDuplicateArray(b);
    mxDestroyArray(retVal);
    if (doClean)
        mxDestroyArray((mxArray*)type);

}