#ifndef SCHEME1_H
#define SCHEME1_H

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
#include "../lib/common/extended_field.h"
#include "scheme4.h"  // scheme4

// ========== Protocol 3  ==========

// 
typedef struct {
    fq_t theta;  // 
    fq_t phi;    //   
    fq_t omega;  // 
} scheme1_field_ext_param;

// 
typedef struct {
    // pk_F =  ()
    // vk_F =  () 
    // rho_i = F for all i ()
    // aux_R0 =  ()
} scheme1_keygen_out;

//   
typedef struct {
    // vk_x:  {g^{m_i}, g^{n_i} | i ∈ [dt]}
    unsigned char **vk_x_m;  // g^{m_i}
    unsigned char **vk_x_n;  // g^{n_i}
    // sigma_i: （）
    ext_field_mat_t *sigma; // σ_i = c(i) for i ∈ [k] ()
    // aux_R:  {m_i, n_i | i ∈ [dt]}
    fq_t *aux_R_m;   // m_i
    fq_t *aux_R_n;   // n_i
} scheme1_probgen_out;

// 
typedef struct {
    // aux'_R:  (a_i, b_i) for i = 0,1,2,...,dt
    fq_t *aux_R_prime_a;  // a_i
    fq_t *aux_R_prime_b;  // b_i
    int valid;            // 
} scheme1_verify_out;

// 
typedef struct {
    fq_t result;  // F(x)
} scheme1_reconstruct_out;

// ========== Protocol 3  ==========

// scheme1_keygen: 、、、
// : F
// : pk_F = ⊥, vk_F = ⊥, {ρ_i}_{i=1}^k = F, aux_R0 = ⊥
void scheme1_keygen(fq_mat_t F, pubpar *par);

// scheme1_probgen: 、、
// : X
// : vk_x, {σ_i}_{i=1}^k (), aux_R
void scheme1_probgen(unsigned char *ga, ext_field_mat_t *c, ext_field_elem_t *b, fq_mat_t X, pubpar *par);

// scheme1_compute: 
// : F, ci, bi
// : yi ()  zi
void scheme1_compute(ext_field_elem_t *yi, ext_field_elem_t *zi, fq_mat_t F, ext_field_mat_t *ci, ext_field_elem_t *bi, pubpar *par);

// scheme1_verify: 
// : ga, y ()  z, 
// : 1，0
int scheme1_verify(unsigned char *ga, ext_field_elem_t *y, ext_field_elem_t *z, pubpar *par);

// scheme1_reconstruct: 
// : ga, y ()  z, 
// : result
void scheme1_reconstruct(fq_t result, unsigned char *ga, ext_field_elem_t *y, ext_field_elem_t *z, pubpar *par);

// ==========  ==========

// ：F_{q^2}
void scheme1_ext_field_interpolate(ext_field_elem_t *phi_coeffs, ext_field_elem_t *y_values, pubpar *par);

// 
void scheme1_eval_extended(ext_field_elem_t *result, fq_mat_t F, ext_field_mat_t *X_ext, ext_field_ctx_t *ctx, pubpar *par);

// 
ext_field_ctx_t* scheme1_get_ext_ctx(void);

// 
void scheme1_init_ext_ctx(pubpar *par);

#endif 