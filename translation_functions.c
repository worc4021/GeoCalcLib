/*  File: translation_functions.c
  
    Program: Geometric calculation library
    
    Description: The geometric calculation library provides a MATLAB callable
    functions to compute the vertex enumeration, facet enumeration and projection
    for convex polyhedra. In addition it provides functions to compute minimal
    representations of convex polyhedra. The geometric library uses the LRS library,
    see http://cgm.cs.mcgill.ca/~avis/C/lrs.html for details.

    Copyright (C) 2015  Rainer Manuel Schaich (rainer.schaich@eng.ox.ac.uk)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
*/

#include <assert.h>
#include <string.h>
#include "lrslib.h"
#include "tmwtypes.h"
#include "translation_functions.h"

#define OUTPUT_PRECISION mxDOUBLE_CLASS
#define printf mexPrintf

size_t pN = 20; /* Number of rows/vertices/rays found between printing out info */

static long lrs_checkpoint_seconds = 0;

static long lrs_global_count = 0;

#ifdef SIGNALS
static void checkpoint ();
static void die_gracefully ();
static void setup_signals ();
static void timecheck ();
#endif

static void lrs_dump_state () {};

struct GMPmat *GMPmat_fromMXArray (const mxArray *pm)
{
       struct GMPmat *A;
       size_t         i, j;

       A = calloc (1, sizeof (*A));
       assert (A != NULL);

       double *ptr = mxGetPr(pm);

       A->m    = mxGetM(pm);
       A->n    = mxGetN(pm);
       A->data = calloc (A->m*A->n, sizeof(*A->data));
       assert ( A->data != NULL );

         for (i=0; i!=A->m; ++i){
           for (j=0; j!=A->n; ++j){
            mpq_init (A->data[i*(A->n) + j]);
            mpq_set_d(A->data[i*(A->n) + j], ptr[i + j*A->m]);
          }
        }

       return (A);
}


mxArray *MXArray_fromGMPmat(const struct GMPmat *A)
{
    mxArray *retVal;
    size_t i, j;
    double *ptr;
    retVal = mxCreateNumericMatrix(A->m, A->n, OUTPUT_PRECISION, mxREAL);
    ptr = mxGetPr(retVal);
    for (i = 0; i != A->m; ++i)
    {
        for (j = 0; j != A->n; ++j)
        {
            ptr[i + j*A->m] = mpq_get_d(A->data[i*(A->n) + j]);
        }
    }
    return retVal;
}


size_t MXArray_to_integer(const mxArray *pm)
{
  if ( !mxIsUint32( pm ) || !mxIsScalar( pm ) )
    mexErrMsgIdAndTxt("MATLAB:myErrors:scalarInt","Input has to be scalar unsigned integer 32.");
  size_t retVal;
  unsigned int *temp;
  temp =  (unsigned int *)mxGetData( pm );
  retVal = (size_t) *temp;
  return retVal;
}

mxArray *VertConcat(const mxArray *A, const mxArray *b)
{
  size_t mA, mB, nA, nB;
  mA = mxGetM(A);
  mB = mxGetM(b);
  mxAssert(mA == mB, "A and b must have the same number of rows!");
  nA = mxGetN(A);
  nB = mxGetN(b);

  mxArray *retVal;
  retVal = mxCreateNumericMatrix(0, 0, OUTPUT_PRECISION, mxREAL);

  double *data, *dA, *dB;
  dA = mxGetPr(A);
  dB = mxGetPr(b);
  data = mxMalloc( mA*(nA+nB)*sizeof(*data) );
  assert ( data != NULL );
  memcpy(data, dA, mA*nA*sizeof(double) );
  memcpy(data + mA*nA, dB, mB*nB*sizeof(double) );

  mxSetPr(retVal, data);
  mxSetM(retVal,mA);
  mxSetN(retVal,nA+nB);
  return retVal;
}

mxArray *VertBreakdown(const mxArray *res)
{
  size_t m, n;
  m = mxGetM(res);
  n = mxGetN(res) - 1;
  
  mxArray *A, *b;
  A = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
  b = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

  double *A_data, *b_data, *res_data;
  res_data = mxGetPr(res);

  A_data = mxMalloc( m*n*sizeof(*A_data) );
  b_data = mxMalloc( m*sizeof(*b_data) );
  
  memcpy(b_data, res_data, m*sizeof(*res_data) );
  memcpy(A_data, res_data + m, m*n*sizeof(*res_data) );

  mxSetPr(A, A_data);
  mxSetM(A, m);
  mxSetN(A, n);

  mxSetPr(b, b_data);
  mxSetM(b, m);
  mxSetN(b, 1);


  mxArray *retVal;
  mwSize ndim=2, dims[]={1, 2}, nsubs=2, subs[2];
  mwIndex index;
  retVal = mxCreateCellArray( ndim, dims );

  subs[0] = 0;
  subs[1] = 0;
  index = mxCalcSingleSubscript(retVal, nsubs, subs);
  mxSetCell (retVal, index , A);

  subs[1] = 1;
  index = mxCalcSingleSubscript(retVal, nsubs, subs);
  mxSetCell (retVal, index , b);

  return retVal;

}

mxArray *createEmptyCell( )
{
  mxArray *A, *b;
  A = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
  b = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

  mxArray *retVal;
  mwSize ndim=2, dims[]={1, 2}, nsubs=2, subs[2];
  mwIndex index;
  retVal = mxCreateCellArray( ndim, dims );

  subs[0] = 0;
  subs[1] = 0;
  index = mxCalcSingleSubscript(retVal, nsubs, subs);
  mxSetCell (retVal, index , A);

  subs[1] = 1;
  index = mxCalcSingleSubscript(retVal, nsubs, subs);
  mxSetCell (retVal, index , b);

  return retVal;

}

int my_lrs_init()
{
  assert ( lrs_mp_init (ZERO, stdin, stdout) );
  lrs_global_count = 0;
  lrs_checkpoint_seconds = 0;
#ifdef SIGNALS
  setup_signals ();
#endif
  return TRUE;
}

struct GMPmat *H2V(struct GMPmat *inp)
{
    lrs_dic *Pv;
    lrs_dat *Qv;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;

      size_t i;
      long col;

      assert( my_lrs_init () == TRUE );

      Qv = lrs_alloc_dat ("LRS globals");
      assert( Qv!= NULL );

      Qv->m = GMPmat_Rows(inp);
      Qv->n = GMPmat_Cols(inp);
      
      output = lrs_alloc_mp_vector (Qv->n);

      lrs_mp_vector num, den;
      num = lrs_alloc_mp_vector(GMPmat_Cols(inp));
      den = lrs_alloc_mp_vector(GMPmat_Cols(inp));

      Pv = lrs_alloc_dic (Qv);
      assert( Pv != NULL );

      
      struct GMPmat *Helper;
      Helper = GMPmat_create(0, GMPmat_Cols(inp), 1);
      
      mpq_t *curRow;
      curRow = calloc(GMPmat_Cols(inp), sizeof(mpq_t));
      assert( curRow != NULL );
      mpq_row_init(curRow, GMPmat_Cols(inp));



      for (i = 1; i <= GMPmat_Rows(inp); ++i)
      {
        GMPmat_getRow(num, den, inp, i-1);
        lrs_set_row_mp(Pv,Qv,i,num,den,GE);
      }

      if( lrs_getfirstbasis (&Pv, Qv, &Lin, TRUE) ){


      for (col = 0L; col < Qv->nredundcol; col++)
        lrs_printoutput (Qv, Lin[col]); 

      do
        {
          for (col = 0L; col <= Pv->d; col++)
            if (lrs_getsolution (Pv, Qv, output, col)) {
              mpz_to_mpq(curRow, output, GMPmat_Cols(Helper));
              Helper = GMPmat_appendRow(Helper, curRow);
              GMPmat_everyNrows(Helper, pN, "vertices/rays");
            }
        }
        while (lrs_getnextbasis (&Pv, Qv, FALSE));

        mpq_row_clean(curRow, GMPmat_Cols(Helper));
        lrs_clear_mp_vector (output, Qv->n);
        lrs_clear_mp_vector (num, Qv->n);
        lrs_clear_mp_vector (den, Qv->n);
        lrs_free_dic (Pv,Qv);
        lrs_free_dat (Qv);

        GMPmat_destroy(inp);
        
        return Helper;

      } else {

        printf("Empty set.\n");
        mpq_row_clean(curRow, GMPmat_Cols(Helper));
        lrs_clear_mp_vector (output, Qv->n);
        lrs_clear_mp_vector (num, Qv->n);
        lrs_clear_mp_vector (den, Qv->n);
        lrs_free_dic (Pv,Qv);
        lrs_free_dat (Qv);
        
        GMPmat_destroy(inp);

        return Helper;

      }
}

struct GMPmat *V2H(struct GMPmat *inp)
{
    lrs_dic *P;
    lrs_dat *Q;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;

    long i;
    long col;

    assert( my_lrs_init () == TRUE );

    Q = lrs_alloc_dat ("LRS globals");
    assert ( Q != NULL );
    Q->m = GMPmat_Rows(inp);
    Q->n = GMPmat_Cols(inp);
    Q->hull = TRUE;
    Q->polytope = TRUE;
    Q->homogeneous = TRUE;

    output = lrs_alloc_mp_vector (Q->n);

    lrs_mp_vector num, den;
    num = lrs_alloc_mp_vector(GMPmat_Cols(inp));
    den = lrs_alloc_mp_vector(GMPmat_Cols(inp));

    P = lrs_alloc_dic (Q);
    assert ( P != NULL );

    struct GMPmat *retMat;
    retMat = GMPmat_create(0, GMPmat_Cols(inp), 1);
      
    mpq_t *curRow;
    curRow = calloc(GMPmat_Cols(inp), sizeof(mpq_t));
    assert( curRow != NULL );
    mpq_row_init(curRow, GMPmat_Cols(inp));


    for (i = 1; i <= GMPmat_Rows(inp); ++i)
    {
      GMPmat_getRow(num, den, inp, i-1);
      if ( zero(num[0]) ) {
          Q->polytope = FALSE;
          Q->homogeneous = FALSE;
        }
      lrs_set_row_mp(P ,Q ,i ,num ,den , GE);
    }

    if ( lrs_getfirstbasis (&P, Q, &Lin, TRUE) ){
    
    for (col = 0L; col < Q->nredundcol; col++)
      lrs_printoutput (Q, Lin[col]);

    do
    {
      for (col = 0; col <= P->d; col++)
        if (lrs_getsolution (P, Q, output, col)){
          mpz_to_mpq(curRow, output, GMPmat_Cols(retMat));
          retMat = GMPmat_appendRow(retMat, curRow);
          GMPmat_everyNrows(retMat, pN, "inequalities");
        }
    }
    while (lrs_getnextbasis (&P, Q, FALSE));

    mpq_row_clean( curRow, GMPmat_Cols(retMat) );
    lrs_clear_mp_vector ( output, Q->n);
    lrs_clear_mp_vector ( num, Q->n);
    lrs_clear_mp_vector ( den, Q->n);
    lrs_free_dic ( P , Q);
    lrs_free_dat ( Q );

    GMPmat_destroy(inp);

    return retMat;

  } else {

    printf("Empty set.\n");
    retMat->empty = 1;
    mpq_row_clean( curRow, GMPmat_Cols(retMat) );
    lrs_clear_mp_vector ( output, Q->n);
    lrs_clear_mp_vector ( num, Q->n);
    lrs_clear_mp_vector ( den, Q->n);
    lrs_free_dic ( P , Q);
    lrs_free_dat ( Q );

    GMPmat_destroy(inp);

    return retMat;

  }

}

struct GMPmat *reducemat(struct GMPmat *inp)
{
    lrs_dic *P;
    lrs_dat *Q;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;

    long i;

    size_t m = GMPmat_Rows(inp);

    size_t *redRows;
    redRows = malloc( m*sizeof(*redRows) );
    assert( redRows != NULL );

    assert( my_lrs_init () == TRUE );

    Q = lrs_alloc_dat ("LRS globals");
    assert ( Q != NULL );
    Q->m = m;
    Q->n = GMPmat_Cols(inp);

    output = lrs_alloc_mp_vector (Q->n);

    lrs_mp_vector num, den;
    num = lrs_alloc_mp_vector(GMPmat_Cols(inp));
    den = lrs_alloc_mp_vector(GMPmat_Cols(inp));

    P = lrs_alloc_dic (Q);
    assert ( P != NULL );

    struct GMPmat *retMat;
    retMat = GMPmat_create(0, GMPmat_Cols(inp), 1);

    for (i = 1; i <= m; ++i)
    {
      GMPmat_getRow(num, den, inp, i-1);
      lrs_set_row_mp(P ,Q ,i ,num ,den , GE);
    }

    if ( lrs_getfirstbasis (&P, Q, &Lin, TRUE) ){

    size_t lastdv = Q->lastdv;
    size_t d = P->d;
    size_t ineq;
    m = P->m_A;

    for (i = lastdv + 1; i <= m + d; ++i)
    {
      ineq = Q->inequality[i - lastdv] - 1;
      if (!checkindex(P, Q, i))
      {
        retMat = GMPmat_appendRow(retMat, mpq_row_extract(inp, ineq));
        GMPmat_everyNrows(retMat, pN, "irredundant inequalities");
      }
    }

    lrs_clear_mp_vector ( output, Q->n);
    lrs_free_dic ( P , Q);
    lrs_free_dat ( Q );

    GMPmat_destroy(inp);

    return retMat;

  } else {

    printf("Empty set.\n");
    retMat->empty = 1;
    lrs_clear_mp_vector ( output, Q->n);
    lrs_free_dic ( P , Q);
    lrs_free_dat ( Q );

    GMPmat_destroy(inp);

    return retMat;

  }
}

struct GMPmat *reducevertices(struct GMPmat *inp)
{
    lrs_dic *P;
    lrs_dat *Q;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;

    long i;

    size_t m = GMPmat_Rows(inp);

    size_t *redRows;
    redRows = malloc( m*sizeof(*redRows) );
    assert( redRows != NULL );

    assert( my_lrs_init () == TRUE );

    Q = lrs_alloc_dat ("LRS globals");
    assert ( Q != NULL );
    Q->m = m;
    Q->n = GMPmat_Cols(inp);
    Q->hull = TRUE;
    Q->polytope = TRUE;

    output = lrs_alloc_mp_vector (Q->n);

    lrs_mp_vector num, den;
    num = lrs_alloc_mp_vector(GMPmat_Cols(inp));
    den = lrs_alloc_mp_vector(GMPmat_Cols(inp));

    P = lrs_alloc_dic (Q);
    assert ( P != NULL );

    struct GMPmat *retMat;
    retMat = GMPmat_create(0, GMPmat_Cols(inp), 1);

    for (i = 1; i <= m; ++i)
    {
      GMPmat_getRow(num, den, inp, i-1);
      lrs_set_row_mp(P ,Q ,i ,num ,den , GE);
    }

    if ( lrs_getfirstbasis (&P, Q, &Lin, TRUE) ){

    size_t lastdv = Q->lastdv;
    size_t d = P->d;
    size_t ineq;
    m = P->m_A;

    for (i = lastdv + 1; i <= m + d; ++i)
    {
      ineq = Q->inequality[i - lastdv] - 1;
      if (!checkindex(P, Q, i))
      {
        retMat = GMPmat_appendRow(retMat, mpq_row_extract(inp, ineq));
        GMPmat_everyNrows(retMat, pN, "irredundant vertices/rays");
      }
    }

    lrs_clear_mp_vector ( output, Q->n);
    lrs_free_dic ( P , Q);
    lrs_free_dat ( Q );

    GMPmat_destroy(inp);
    
    return retMat;

  } else {

    printf("Empty set.\n");
    retMat->empty = 1;
    lrs_clear_mp_vector ( output, Q->n);
    lrs_free_dic ( P , Q);
    lrs_free_dat ( Q );

    GMPmat_destroy(inp);
    
    return retMat;
    
  }
}


void GMPmat_validate_indices (const struct GMPmat *A, size_t r, size_t c)
     {
       assert (A != NULL);
       assert (r < A->m);
       assert (c < A->n);
     }

void GMPmat_getValue (mpq_t rop, const struct GMPmat *A, size_t r, size_t c)
{
       GMPmat_validate_indices (A, r, c);
       mpq_set (rop, A->data[r*(A->n) + c]);
}

void GMPmat_setValue (const struct GMPmat *A, size_t r, size_t c, double val)
{
       GMPmat_validate_indices (A, r, c);
       mpq_set_d(A->data[r*(A->n) + c], val);
}

size_t GMPmat_Cols (const struct GMPmat *A)
     {
       assert (A != NULL);

       return (A->n);
     }


size_t GMPmat_Rows (const struct GMPmat *A)
     {
       assert (A != NULL);

       return (A->m);
     }


struct GMPmat *GMPmat_create (size_t m, size_t n, int init)
     {
       struct GMPmat *A;
       size_t         i, j;

       A = calloc (1, sizeof (*A));
       assert (A != NULL);

       A->empty = 0;
       A->m    = m;
       A->n    = n;
       A->data = calloc (m*n, sizeof(*A->data));
       assert (A->data != NULL );

       if (init != 0)
       {
         for (i=0; i!=m; ++i){
           for (j=0; j!=n; ++j){
            mpq_init (A->data[i*(A->n) + j]);
            
          }
        }
       }

       return (A);
}

#ifdef DEBUG
void GMPmat_print(const struct GMPmat *A)
{
    const char *val;

    val = getenv ("SILENCE");
    if (val != NULL)
        return;

    assert( A!= NULL );
    fprintf(stdout, "\n");
    for (size_t i = 0; i < GMPmat_Rows(A); ++i)
    {
        for (size_t j = 0; j < GMPmat_Cols(A); ++j)
        {
            mpq_out_str(stdout, 10, A->data[i*(A->n) + j]);
            fprintf(stdout, " ");
        }
        fprintf(stdout, "\n" );
    }
    fprintf(stdout, "\n" );
}
#endif /* DEBUG */

void GMPmat_destroy (struct GMPmat *A)
     {
       assert (A != NULL);
       size_t i, j;
       for (i = 0; i < GMPmat_Rows(A); ++i)
       {
        for (j = 0; j < GMPmat_Cols(A); ++j)
        {
            mpq_clear(A->data[i*(A->n) + j]);
        }
       }
       if (A->data != NULL) free (A->data);
       A->data = NULL;
       A->n   = 0;
       A->m   = 0;
       free (A);
     }

void GMPmat_getRow(mpz_t *ropN, mpz_t *ropD, struct GMPmat *A, size_t r)
{
    assert( r < A->m );
    assert( ropN != NULL && ropD != NULL );
    size_t i;
    for (i = 0; i < GMPmat_Cols(A); ++i)
    {
        mpz_set(ropN[i], mpq_numref(A->data[r*(A->n) + i]));
        mpz_set(ropD[i], mpq_denref(A->data[r*(A->n) + i]));
    }
}

#ifdef DEBUG
void GMPmat_printRow(const struct GMPmat *A, size_t r)
{
  assert ( r < A->m );
  size_t i, n;
  n = GMPmat_Cols(A);
  mpq_t curVal;
  mpq_init(curVal);
  fprintf(stdout, "\n" );
  for ( i = 0; i < n ; ++i){
    GMPmat_getValue(curVal, A, r, i);
    mpq_out_str(stdout, 10, curVal);
    fprintf(stdout, " ");
  }
  fprintf(stdout, "\n" );
  mpq_init(curVal);
}
#endif /* DEBUG */

struct GMPmat *GMPmat_dropCols(struct GMPmat *A, size_t d)
{
  assert ( A != NULL );
  assert ( d < GMPmat_Cols(A) );
  size_t m, n, i, j;
  m = GMPmat_Rows(A);
  n = GMPmat_Cols(A) - d;

  struct GMPmat *retVal;
  retVal = GMPmat_create( m, n, 0 );

  for (i=0; i!=m; ++i){
    for (j=0; j!=n; ++j){
      mpq_init ( retVal->data[i*n + j] );
      mpq_set ( retVal->data[i*n + j], A->data[i*(A->n) + j] );
    }
  }

  GMPmat_destroy(A);
  return retVal;
}

void GMPmat_everyNrows(struct GMPmat *A, size_t N, char *type)
{
  #ifndef NOINFO
  assert ( A != NULL );
  if ( GMPmat_Rows(A) % N == 0 )
  {
    printf( "So far %zu %s found.\n", (GMPmat_Rows(A)/N)*N, type);
  }
  #endif /* NOINFO */
}

void GMPmat_invertSignForFacetEnumeration(struct GMPmat *A)
{
  assert ( A != NULL );
  size_t i, j, m, n;
  m = GMPmat_Rows(A);
  n = GMPmat_Cols(A);
  mpq_t holder;
  mpq_init(holder);
  for (i = 0; i < m; i++)
  {
    for (j = 1; j < n; j++)
    {
      GMPmat_getValue(holder, A, i, j);
      mpq_neg(A->data[i*(A->n) + j], holder);
    }
  }
  mpq_clear(holder);
}

void mpz_row_clean(mpz_t *row, size_t m)
{
    assert(row != NULL);
    size_t i;
    for (i = 0; i < m; ++i)
    {
        mpz_clear(row[i]);
    }
    free(row);
}

void mpq_row_clean(mpq_t *row, size_t m)
{
    assert(row != NULL);
    size_t i;
    for (i = 0; i < m; ++i)
    {
        mpq_clear(row[i]);
    }
    free(row);
}

void mpz_row_init(mpz_t *row, size_t m)
{
    assert(row != NULL);
    size_t i;
    for (i = 0; i < m; ++i)
    {
        mpz_init(row[i]);
    }
}

void mpq_row_init(mpq_t *row, size_t m)
{
    assert(row != NULL);
    size_t i;
    for (i = 0; i < m; ++i)
    {
        mpq_init(row[i]);
    }
}

void mpz_norm(mpz_t norm, mpz_t *row, size_t m)
{
    assert( row != NULL );
    mpz_t help;
    mpz_init(help);
    size_t i;
    for (i = 0; i < m; ++i)
    {
        mpz_addmul(help, row[i], row[i]);
    }
    mpz_sqrt(norm, help);
    mpz_clear(help);
}

struct GMPmat *GMPmat_appendRow(struct GMPmat *A, mpq_t *row)
{
    assert( A != NULL && row != NULL );
    struct GMPmat *retVal;
    size_t m,n, i, j;
    m = GMPmat_Rows(A) + 1;
    n = GMPmat_Cols(A);
    retVal = GMPmat_create( m, n, 0);
    assert( m != 0 );

    for (i = 0; i < (m-1); ++i)
    {
        for (j = 0; j < n; ++j)
        {
            mpq_init(retVal->data[i*n + j]);
            mpq_set(retVal->data[i*n + j],A->data[i*n + j]);
        }
    }
    for (i = 0; i < n; ++i)
    {
        mpq_init(retVal->data[(m-1)*n + i]);
        mpq_set(retVal->data[(m-1)*n + i], row[i]);
    }
    GMPmat_destroy(A);
    return retVal;
}

void mpz_to_mpq(mpq_t *rop, mpz_t *op, size_t m)
{
    size_t i;
    assert( rop != NULL && op != NULL);
    if (!mpz_sgn(op[0]))
    {
        mpz_t norm;
        mpz_init(norm);
        mpz_norm(norm, op, m);
        for (i = 0; i < m; ++i)
        {
            mpq_set_num(rop[i], op[i]);
            mpq_set_den(rop[i], norm);
            mpq_canonicalize(rop[i]);
        }
        mpz_clear(norm);
    }else{
        for (i = 0; i < m; ++i)
        {
            mpq_set_num(rop[i], op[i]);
            mpq_set_den(rop[i], op[0]);
            mpq_canonicalize(rop[i]);
        }
    }
}

#ifdef DEBUG
void mpz_row_print(mpz_t *row, size_t n)
{
    assert ( row != NULL );
    size_t i;
    for (i = 0; i < n; ++i)
    {
        mpz_out_str(stdout, 10, row[i]);
        fprintf(stdout, " " );
    }
    fprintf(stdout, "\n" );
}

void mpq_row_print(mpq_t *row, size_t n)
{
    assert ( row != NULL );
    size_t i;
    for (i = 0; i < n; ++i)
    {
        mpq_out_str(stdout, 10, row[i]);
        fprintf(stdout, " " );
    }
    fprintf(stdout, "\n" );
}



void mpz_print(mpz_t op)
{
    assert( op != NULL );
    fprintf(stdout, "\n");
    mpz_out_str(stdout, 10, op);
    fprintf(stdout, "\n");
}

void mpq_print(mpq_t op)
{
    assert( op != NULL );
    fprintf(stdout, "\n");
    mpq_out_str(stdout, 10, op);
    fprintf(stdout, "\n");
}

void mpz_print_product(mpz_t numA, mpz_t denA, mpz_t numB, mpz_t denB)
{
    assert( numA != NULL && denA != NULL && numB != NULL && denB != NULL );
    mpq_t a,b, res;
    mpq_init(a);
    mpq_init(b);
    mpq_init(res);
    
    mpq_set_num(a, numA);
    mpq_set_den(a, denA);
    mpq_canonicalize(a);
    
    mpq_set_num(b, numB);
    mpq_set_den(b, denB);
    mpq_canonicalize(b);
    
    mpq_mul(res, a, b);

    mpq_print(res);
    
    mpq_clear(res);
    mpq_clear(b);
    mpq_clear(a);
}
#endif /* DEBUG */


mpq_t *mpq_row_extract(const struct GMPmat *A, size_t r)
{
    assert ( A != NULL && r < GMPmat_Rows(A) );
    mpq_t *retVal;
    size_t i, n; 
    n = GMPmat_Cols(A);
    retVal = malloc( n*sizeof(*retVal) );
    for ( i = 0; i < n; ++i )
    {
        mpq_init( retVal[i] );
        GMPmat_getValue( retVal[i], A, r, i);
    }
    return retVal;

}

void lprint(long *array, size_t l)
{
    fprintf(stdout, "\n");
    size_t i;
    for (i = 0; i < l; i++)
        fprintf(stdout, "%ld ", array[i]);
    fprintf(stdout, "\n");
}


#ifdef SIGNALS

/*
   If given a signal
   USR1            print current cobasis and continue
   TERM            print current cobasis and terminate
   INT (ctrl-C) ditto
   HUP                     ditto
 */
static void
setup_signals ()
{
  errcheck ("signal", signal (SIGTERM, die_gracefully));
  errcheck ("signal", signal (SIGALRM, timecheck));
  errcheck ("signal", signal (SIGHUP, die_gracefully));
  errcheck ("signal", signal (SIGINT, die_gracefully));
  errcheck ("signal", signal (SIGUSR1, checkpoint));
}

static void
timecheck ()
{
  lrs_dump_state ();
  errcheck ("signal", signal (SIGALRM, timecheck));
  alarm (lrs_checkpoint_seconds);
}

static void
checkpoint ()
{
  lrs_dump_state ();
  errcheck ("signal", signal (SIGUSR1, checkpoint));
}

static void
die_gracefully ()
{
  lrs_dump_state ();

  exit (1);
}

#endif