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
    mxArray *retVal;
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


    retVal = vertexReduction(V, type);

    mwSize nsubs=2, subs[2];
    mwIndex index;
    subs[0] = 0;
    subs[1] = 0;
    index = mxCalcSingleSubscript(retVal, nsubs, subs);
    plhs[0] = mxDuplicateArray(mxGetCell(retVal, index));

    subs[1] = 1;
    index = mxCalcSingleSubscript(retVal, nsubs, subs);

    if (nlhs==2)
        plhs[1] = mxDuplicateArray(mxGetCell(retVal, index));
    mxDestroyArray(retVal);

    if (doClean)
        mxDestroyArray((mxArray*)type);


}