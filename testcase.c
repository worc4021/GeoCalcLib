/*
File: testcase.c
Description:
The interface used here was developed by implementing calls to the 
lrs library as they are outlined in chdemo.c, vedemo.c.

In this file a matrix of vertices is created, the vertices lie on the
unit circle, i.e. each vertex v_i = {sin(2*pi*i/N),cos(2*pi*i/N)}.
For these vertices a facet enumeration is can be trivially obtained.
For the purpose of illustrating the malfunction we compute a facet 
enumeration of the convex hull of {v_i+offset} in such a way that
the origin is not contained in the polytope. The facet enumeration
produces wrong results.
To verify that the function calls otherwise do what is expected we
calculate a facet enumeration of the convex hull of {v_i} (with no
offset), the facet enumeration produces the correct result, which
then is shifted.
conv{v_i} = {A*x<=b} --> conv{v_i+offset} = {Ax<=b-A*offset}
This produces vastly different results.

All matrices can be printed to the stream either as floating point 
numbers or in the GMP rational format (I avoided this because the
rational representation is rather untransparent to a human.
(You can use it by calling GMPmat_print instead of GMPmat_printAsDouble).
*/


#include "translation_functions.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

void GMPmat_printAsDouble(const struct GMPmat *A);
void GMPmat_shiftHRep(const struct GMPmat *A, const struct GMPmat *y);
// void GMPmat_setRow(const struct GMPmat *A, const mpq_t *row, size_t r);

int main(int argc, char const *argv[])
{
	size_t i,N = 30; /* N = number of vertices */
	double t;
	struct GMPmat *Vertices, *Offset;
	
	Vertices = GMPmat_create(N,3,1);
	
	double offset[] = {5.,5.};
	Offset = GMPmat_create(2,1,1);
	GMPmat_setValue(Offset,0,0,offset[0]);
	GMPmat_setValue(Offset,1,0,offset[1]);

	/* Creating the vertices with an offset */
	for (i = 0; i < N; ++i)
	{
		t = ((double)i/(double)N)*2.*M_PI;
		GMPmat_setValue(Vertices,i,0,1.);
		GMPmat_setValue(Vertices,i,1,sin(t)+offset[0]);
		GMPmat_setValue(Vertices,i,2,cos(t)+offset[1]);
	}

	// GMPmat_print(Vertices);
	Vertices = GMPmat_lift(Vertices);
	// GMPmat_print(Vertices);
	Vertices = H2V(Vertices);
	GMPmat_print(Vertices);
	Vertices = GMPmat_unlift(Vertices);
	GMPmat_printAsDouble(Vertices);
	// fprintf(stdout,"Initial Vertices:\n");
	// GMPmat_printAsDouble(Vertices);

	// /* Reducing (non-existent) redundant vertices using the LRS checkindex function */
	// Vertices = reducevertices(Vertices);

	// fprintf(stdout,"Reduced Vertices:\n");
	// GMPmat_printAsDouble(Vertices);

	// /* Calculating the facet enumeration using analogue procedures as those in chdemo.c */
	// Vertices = V2H(Vertices);
	// GMPmat_invertSignForFacetEnumeration(Vertices);

	// fprintf(stdout,"Facet description:\n");
	// GMPmat_printAsDouble(Vertices);

	// /* Reducing redundant inequalities using analogue procedures as those in redund_main in lrslib.c */
	// Vertices = reducemat(Vertices);
	// fprintf(stdout,"Reduced H-representation:\n");
	// GMPmat_printAsDouble(Vertices);

	// GMPmat_destroy(Vertices);

	// /* Same all over again using no offset to verify that the function calls are correct otherwise */
	// Vertices = GMPmat_create(N,3,1);

	// for (i = 0; i < N; ++i)
	// {
	// 	t = ((double)i/(double)N)*2.*M_PI;
	// 	GMPmat_setValue(Vertices,i,0,1.);
	// 	GMPmat_setValue(Vertices,i,1,sin(t));
	// 	GMPmat_setValue(Vertices,i,2,cos(t));
	// }

	// fprintf(stdout,"Vertices without offset:\n");
	// GMPmat_printAsDouble(Vertices);
	
	// Vertices = V2H(Vertices);
	// GMPmat_invertSignForFacetEnumeration(Vertices);

	// fprintf(stdout,"H-representation for unshifted Vertices:\n");
	// GMPmat_printAsDouble(Vertices);

	// Vertices = reducemat(Vertices);
	// fprintf(stdout,"Reduced H-representation:\n");
	// GMPmat_printAsDouble(Vertices);

	// GMPmat_shiftHRep(Vertices,Offset);
	
	// fprintf(stdout, "(b-A*offset,-A) \n");
	// GMPmat_printAsDouble(Vertices);

	// Vertices = H2V(Vertices);
	// fprintf(stdout, "Vertex enumeration again:\n");
	// GMPmat_printAsDouble(Vertices);	

	// GMPmat_destroy(Vertices);
	// GMPmat_destroy(Offset);
	return 0;
}


void GMPmat_printAsDouble(const struct GMPmat *A)
{
    assert( A!= NULL );
    double curVal;
    size_t i, j;
    fprintf(stdout, "\n");
    fprintf(stdout, "Matrix has %zu rows and %zu columns.\n", GMPmat_Rows(A), GMPmat_Cols(A));
    for (i = 0; i < GMPmat_Rows(A); ++i)
    {
        for (j = 0; j < GMPmat_Cols(A); ++j)
        {
        	curVal = mpq_get_d(A->data[i*(A->n) + j]);
            fprintf(stdout, "%5.2f ",curVal);
        }
        fprintf(stdout, "\n" );
    }
    fprintf(stdout, "\n" );
}

void GMPmat_shiftHRep(const struct GMPmat *A, const struct GMPmat *y)
{
	mpq_t sum, prod, curElemA, curElemY;
	mpq_init(sum);
	mpq_init(prod);
	mpq_init(curElemA);
	mpq_init(curElemY);

	mpq_t *currentRow;
	currentRow = malloc(GMPmat_Cols(A)*sizeof(*currentRow));
	mpq_row_init(currentRow, GMPmat_Cols(A) );

	size_t i,j;
	for (i = 0; i < GMPmat_Rows(A); ++i)
	{
		mpq_set_si(sum, 0L, 1L);
		currentRow = mpq_row_extract(A,i);
		for (j = 1; j < GMPmat_Cols(A); ++j)
		{
			GMPmat_getValue(curElemA, A, i, j);
			GMPmat_getValue(curElemY, y, j-1, 0);
			mpq_mul(prod, curElemA, curElemY);
			mpq_sub(sum, sum, prod);
		}
		mpq_add(currentRow[0],currentRow[0],sum);
		GMPmat_setRow(A, currentRow, i);
	}

	mpq_clear(sum);
	mpq_clear(prod);
	mpq_clear(curElemA);
	mpq_clear(curElemY);
}


// void GMPmat_setRow(const struct GMPmat *A, const mpq_t *row, size_t r)
// {
// 	assert( r < GMPmat_Rows(A) );
// 	size_t i;
// 	for (i = 0; i<GMPmat_Cols(A); ++i)
// 		mpq_set(A->data[r*(A->n) + i], row[i]);
// }