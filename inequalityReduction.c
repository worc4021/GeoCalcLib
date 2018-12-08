#include "mex.h"
#include "mainFunctions.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs<1)
        mexErrMsgTxt("Not enough inputs.");
        
    const mxArray *A;
    const mxArray *b;
    mxArray *retVal;
    

    if (nrhs<2)
        mexErrMsgTxt("Inequality reduction requires A,b to define {A*x<=b}.");

    A = prhs[0];
    b = prhs[1];

    retVal = ineqReduction(A,b);

    mwSize nsubs=2, subs[2];
    mwIndex index;
    subs[0] = 0;
    subs[1] = 0;
    index = mxCalcSingleSubscript(retVal, nsubs, subs);

    plhs[0] = mxDuplicateArray(mxGetCell(retVal, index));

    subs[1] = 1;
    index = mxCalcSingleSubscript(retVal, nsubs, subs);
    if (nlhs>1)
        plhs[1] = mxDuplicateArray(mxGetCell(retVal, index));

    mxDestroyArray(retVal);


}
