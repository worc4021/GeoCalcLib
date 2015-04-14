#include <gmp.h>
#include "mex.h"
#include "matrix.h"

struct GMPmat{
    mpq_t *data;
    size_t m, n;
};


struct GMPmat *GMPmat_fromMXArray (const mxArray *pm);
void GMPmat_print(const struct GMPmat *A);



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	struct GMPmat *helper;
	helper = GMPmat_fromMXArray(prhs[0]);
	GMPmat_print(helper);
}


struct GMPmat *GMPmat_fromMXArray (const mxArray *pm)
{
       struct GMPmat *A;
       size_t         i, j;

       A = calloc (1, sizeof (*A));
       // assert (A != NULL);

       double *ptr = mxGetPr(pm);

       A->m    = mxGetM(pm);
       A->n    = mxGetN(pm);
       A->data = calloc (A->m*A->n, sizeof(*A->data));
       // assert ( A->data != NULL );

         for (i=0; i!=A->m; ++i){
           for (j=0; j!=A->n; ++j){
            mpq_init (A->data[i*(A->n) + j]);
            mpq_set_d(A->data[i*(A->n) + j], ptr[i + j*A->m]);
          }
        }

       return (A);
}

void GMPmat_print(const struct GMPmat *A)
{
    const char *val;

    // assert( A!= NULL );
    fprintf(stdout, "\n");
    for (size_t i = 0; i < A->m; ++i)
    {
        for (size_t j = 0; j < A->n; ++j)
        {
            mpq_out_str(stdout, 10, A->data[i*(A->n) + j]);
            fprintf(stdout, " ");
        }
        fprintf(stdout, "\n" );
    }
    fprintf(stdout, "\n" );
}