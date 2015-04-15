/*	File: mainFunctions.h
	
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

#include <gmp.h>
#include "mex.h"
#include "matrix.h"

struct GMPmat{
    mpq_t *data;
    size_t m, n;
};

/* Interface between mxArray and GMPmat structures */

struct GMPmat *GMPmat_fromMXArray (const mxArray *pm);
mxArray *MXArray_fromGMPmat(const struct GMPmat *A);
size_t MXArray_to_integer(const mxArray *pm);
mxArray *VertConcat(const mxArray *A, const mxArray *b);
mxArray *VertBreakdown(const mxArray *res);

/* Interface to the LRS library */

int my_lrs_init();
struct GMPmat *H2V(struct GMPmat *inp);
struct GMPmat *V2H(struct GMPmat *inp);
struct GMPmat *reducemat(struct GMPmat *inp);
struct GMPmat *reducevertices(struct GMPmat *inp);

/* Functions to work with GMPmat structures */

void GMPmat_getRow(mpz_t *ropN, mpz_t *ropD, struct GMPmat *A, size_t r);
void GMPmat_destroy (struct GMPmat *A);

#ifdef DEBUG
void GMPmat_print(const struct GMPmat *A);
#endif /* DEBUG */

struct GMPmat *GMPmat_create (size_t m, size_t n, int init);
size_t GMPmat_Rows (const struct GMPmat *A);
size_t GMPmat_Cols (const struct GMPmat *A);
void GMPmat_setValue (const struct GMPmat *A, size_t r, size_t c, double val);
void GMPmat_getValue (mpq_t rop, const struct GMPmat *A, size_t r, size_t c);
void GMPmat_validate_indices (const struct GMPmat *A, size_t r, size_t c);

#ifdef DEBUG
void GMPmat_printRow(const struct GMPmat *A, size_t r);
#endif /* DEBUG */

struct GMPmat *GMPmat_dropCols(struct GMPmat *A, size_t d);
void GMPmat_everyNrows(struct GMPmat *A, size_t N, char *type);
void GMPmat_invertSignForFacetEnumeration(struct GMPmat *A);

/* Functions to work with GMP data types */

void mpz_to_mpq(mpq_t *rop, mpz_t *op, size_t m);
struct GMPmat *GMPmat_appendRow(struct GMPmat *A, mpq_t *row);
void mpz_norm(mpz_t norm, mpz_t *row, size_t m);
void mpq_row_init(mpq_t *row, size_t m);
void mpz_row_init(mpz_t *row, size_t m);
void mpq_row_clean(mpq_t *row, size_t m);
void mpz_row_clean(mpz_t *row, size_t m);

#ifdef DEBUG
void mpz_row_print(mpz_t *row, size_t n);
void mpq_row_print(mpq_t *row, size_t n);
void mpq_print(mpq_t op);
void mpz_print(mpz_t op);
void mpz_print_product(mpz_t numA, mpz_t denA, mpz_t numB, mpz_t denB);
#endif /* DEBUG */

mpq_t *mpq_row_extract(const struct GMPmat *A, size_t r);
void lprint(long *array, size_t l);