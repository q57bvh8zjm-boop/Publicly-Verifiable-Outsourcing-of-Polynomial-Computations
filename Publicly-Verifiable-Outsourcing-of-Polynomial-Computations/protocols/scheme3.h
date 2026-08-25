#ifndef SCHEME3_H
#define SCHEME3_H

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sodium.h>
#include "fmpz.h"
#include "fq.h" 
#include "fmpz_poly.h"
#include "fq_poly.h"
#include "fq_vec.h"
#include "fq_mat.h"
#include "gmp.h"
#include "../lib/common/mpoly.h"

// ========== Protocol 9 (Paillier-based HSS)  ==========

// Paillier ( HSSElg.h Elgamal_PK)
typedef struct {
    fmpz_t f;      // f = n + 1
    fmpz_t g;      // generator
    fmpz_t h;      // h = g^sk mod n^2
    fmpz_t n;      // n = p * q
    fmpz_t n2;     // n^2
    fmpz_t mu;     // (For compatibility, though GroupHSS uses DDLog)
} paillier_pk_t;

// Paillier/HSS ( HSSElg.h HSS_CT)
// HSS_CT is array<Elgamal_CT, 2>, where Elgamal_CT is array<ZZ, 2>
// So we store 4 integers: c[0][0], c[0][1], c[1][0], c[1][1]
typedef struct {
    fmpz_t c[2][2]; 
} paillier_ct_t;

typedef struct {
    fmpz_t ek;     // 
} hss_ek_t;

// HSS_MV is array<ZZ, 2>
typedef struct {
    fmpz_t mv[2];     // 
} hss_mv_t;

// 
typedef struct {
    paillier_pk_t pk;      // HSS
    hss_ek_t ek[2];        // 2
    fmpz_t sk;
} scheme3_keygen_out;

//   
typedef struct {
    unsigned char **vk_x;
    
    paillier_ct_t *ct_inputs;
    
    fq_t *v_values;
} scheme3_probgen_out;

typedef struct {
    hss_mv_t *mv;
} scheme3_compute_out;

// 
typedef struct {
    fq_mat_t aux_R_prime;
    int valid;
} scheme3_verify_out;

// 
typedef struct {
    fq_mat_t result;
} scheme3_reconstruct_out;

// ========== Protocol 9  ==========

// scheme3_keygen: 
// : HSSpk, 2{ek_0, ek_1}
void scheme3_keygen(scheme3_keygen_out *out, fq_mat_t F, pubpar *par);

// scheme3_probgen: 
// : pk, x
void scheme3_probgen(scheme3_probgen_out *out, paillier_pk_t *pk, 
                       fq_mat_t x, pubpar *par);

// scheme3_compute: 
// : id b (01), ek_b, F, HSS{ct_j}
// : HSS{mv_{b,j}}
void scheme3_compute(scheme3_compute_out *out, int server_id, 
                       hss_ek_t *ek_b, paillier_pk_t *pk, fq_mat_t F,
                       paillier_ct_t *ct_inputs, pubpar *par);

// scheme3_verify: 
// : vk_x, , v_values, F, pk
int scheme3_verify(scheme3_verify_out *out, unsigned char **vk_x,
                     scheme3_compute_out *compute_results, fq_t *v_values,
                     fq_mat_t F, paillier_pk_t *pk, pubpar *par);

// scheme3_reconstruct: 
// : aux'_R=f(0)
// : F(x)
void scheme3_reconstruct(scheme3_reconstruct_out *out, fq_mat_t aux_R_prime, pubpar *par);

// ========== Protocol 9  ==========

// Paillier
void paillier_keygen(paillier_pk_t *pk, fmpz_t sk, int security_bits);

// Paillier
void paillier_encrypt(paillier_ct_t *ct, paillier_pk_t *pk, fmpz_t m);

// Paillier
void paillier_decrypt(fmpz_t m, paillier_pk_t *pk, fmpz_t sk, paillier_ct_t *ct);

void paillier_add(paillier_ct_t *result, paillier_pk_t *pk, 
                  paillier_ct_t *ct1, paillier_ct_t *ct2);

void paillier_scalar_mul(paillier_ct_t *result, paillier_pk_t *pk,
                         paillier_ct_t *ct, fmpz_t scalar);

void hss_gen(paillier_pk_t *pk, fmpz_t sk, hss_ek_t *ek0, hss_ek_t *ek1, int security_bits);

void hss_input(paillier_ct_t *ct, paillier_pk_t *pk, fq_t x, pubpar *par);

// ct_inputsm
void hss_evaluate(hss_mv_t *result, int server_id, hss_ek_t *ek_b,
                  paillier_pk_t *pk, fq_mat_t F, paillier_ct_t *ct_inputs, pubpar *par);

void hss_decrypt(fq_t result, paillier_pk_t *pk, hss_mv_t *mv0, hss_mv_t *mv1, pubpar *par);

void hss_mul_paillier(hss_mv_t *result, int server_id, paillier_pk_t *pk, paillier_ct_t *Ix, hss_mv_t *My, int *prf_key);

void scheme3_compute_lambda_prime(fq_t *lambda_prime, fq_t *v_values, int d, pubpar *par);

// ========== Protocol 9  ==========

// Paillier
void paillier_pk_clear(paillier_pk_t *pk);

// Paillier
void paillier_ct_clear(paillier_ct_t *ct);

// HSS
void hss_ek_clear(hss_ek_t *ek);

// HSS
void hss_mv_clear(hss_mv_t *mv);

// keygen
void scheme3_keygen_out_clear(scheme3_keygen_out *out);

// probgen
void scheme3_probgen_out_clear(scheme3_probgen_out *out, pubpar *par);

// compute
void scheme3_compute_out_clear(scheme3_compute_out *out, pubpar *par);

// verify
void scheme3_verify_out_clear(scheme3_verify_out *out, pubpar *par);

// reconstruct
void scheme3_reconstruct_out_clear(scheme3_reconstruct_out *out, pubpar *par);

#endif // SCHEME3_H
