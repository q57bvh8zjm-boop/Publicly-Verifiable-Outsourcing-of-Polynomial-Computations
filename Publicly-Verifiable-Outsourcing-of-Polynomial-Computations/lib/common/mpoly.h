#ifndef MPOLY_H
#define MPOLY_H

#include <time.h>
#include "fmpz.h"
#include "fq.h"
#include "fmpz_poly.h"
#include "fq_poly.h"
#include "stdio.h"
#include "gmp.h"
#include "stdlib.h"
#include "string.h"
#include "fq_vec.h"
#include "fq_mat.h"
#include "flint.h"
#include <sodium.h>

// ===========  ===========

// 
typedef struct {
    int m;          // 
    int d;          // 
    int t;          // 
    int k;          // 
    int iN;         // 
    int L;          // Protocol 8: 
    int ***I;       // 
    fmpz_t p;       // 
    fq_ctx_t Fp;    // 
    fq_t *SID;      // ID
} pubpar;

// ===========  ===========

void Eval(fq_t y, fq_mat_t F, fq_mat_t X, pubpar *par);

// Number of coefficients in a dense m-variate polynomial of total degree <= d.
// Coefficients use graded lexicographic order, represented by nondecreasing
// variable-index tuples within each degree. Returns -1 on invalid input or
// overflow.
slong mpoly_num_coefficients(int m, int d);

void IntPoly(fq_poly_t f, fq_t *Y, pubpar *par);

void fq_mat_scal_mul(fq_mat_t op, fq_t c, fq_ctx_t Fp);

// ===========  ===========

// fqfmpz
void fq2fmpz(fmpz_t out, fq_t in, fq_ctx_t ctx);

// unsigned charfmpz
void chars2fmpz(fmpz_t zo, unsigned char r[]);

// fmpzunsigned char
void fmpz2chars(unsigned char r[], fmpz_t zo);

#endif // MPOLY_H 
