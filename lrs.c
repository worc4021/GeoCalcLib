#include "mex.h"
#include "mainFunctions.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs<1)
        mexErrMsgTxt("Not enough arguments");

    const mxArray *inputData;
    mxArray *retVal;
    
    inputData = prhs[0];
    retVal = lrsCall(inputData);

    mwSize nsubs=3, subs[3];
    mwIndex index;
    subs[0] = 0;
    subs[1] = 0;
    index = mxCalcSingleSubscript(retVal, nsubs, subs);
    plhs[0] = mxDuplicateArray(mxGetCell(retVal, index));

    subs[1] = 1;
    index = mxCalcSingleSubscript(retVal, nsubs, subs);
    if (nlhs>1)
        plhs[1] = mxDuplicateArray(mxGetCell(retVal, index));

    subs[1] = 2;
    index = mxCalcSingleSubscript(retVal, nsubs, subs); 
    if (nlhs>2)
        plhs[2] = mxDuplicateArray(mxGetCell(retVal, index));

    mxDestroyArray(retVal);
}