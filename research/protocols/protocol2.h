#ifndef PROTOCOL2_H
#define PROTOCOL2_H

#include "../lib/common/mpoly.h"

// =========== Protocol 2  ===========

// Protocol 2 （F_{q^2}）
typedef struct {
    int m;          // 
    int d;          // 
    int t;          // 
    int k;          // 
    int iN;         // 
    fmpz_t p;       // 
    fq_ctx_t Fq2;   // F_{q^2}
    fq_t *SID;      // ID
} pubpar_p2;

// =========== Protocol 2  ===========

// 
// : F - 
// : （）
void protocol2_keygen(fq_mat_t F);

// 
// : X - , par - 
// : alpha - , c - 
void protocol2_probgen(fq_t *alpha, fq_mat_t *c, fq_mat_t X, pubpar_p2 *par);

// 
// : F - , ci - , par - 
// : yi - 
void protocol2_compute(fq_t yi, fq_mat_t F, fq_mat_t ci, pubpar_p2* par);

// 
// : alpha - , y - , par - 
// : （1，0）
int protocol2_verify(fq_t alpha, fq_t *y, pubpar_p2 *par);

// =========== Protocol 2  ===========

// ：y = F(X)
void protocol2_eval(fq_t y, fq_mat_t F, fq_mat_t X, pubpar_p2 *par);

// ：ff(SID) = Y
void protocol2_interpolate(fq_poly_t f, fq_t *Y, pubpar_p2 *par);

#endif // PROTOCOL2_H 