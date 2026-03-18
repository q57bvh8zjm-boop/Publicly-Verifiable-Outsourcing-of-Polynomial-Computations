#ifndef SCHEME2_H
#define SCHEME2_H

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

// ========== Protocol 4 (Π₂)  ==========

// 
typedef struct {
    // pk_F = ⊥ (scheme2.md)
    // vk_F = ⊥ (scheme2.md) 
    // rho_i = F for all i (，F)
    // aux_R0 = ⊥ (scheme2.md)
} scheme2_keygen_out;

//   
typedef struct {
    // vk_x:  {g^{λ'_j} | j ∈ {0,1,...,d}}
    unsigned char **vk_x;  // g^{λ'_j}，d+1
    
    // sigma_i: ，id+1
    // sigma_i = {σ_{i_j} | j ∈ {0,1,...,d}}
    fq_mat_t **sigma; // σ_{i_j} = c_{v_j}(i) for i ∈ [k], j ∈ [d]; σ_{i_0} = c(i)
    
    // aux_R = ⊥ (scheme2.md)
    
    // v_j，
    fq_t *v_values;  // v_1, v_2, ..., v_d
} scheme2_probgen_out;

// 
typedef struct {
    // π_i = {π_{i_j} | j ∈ {0,1,...,d}} for server i
    fq_mat_t *pi; // π_{i_j} = F(σ_{i_j})，d+1，1x1
} scheme2_compute_out;

// 
typedef struct {
    // aux'_R: f(0) if verification succeeds
    fq_mat_t aux_R_prime;  // f(0)，1x1
    int valid;             // ：1，0
} scheme2_verify_out;

// 
typedef struct {
    fq_mat_t result;  // F(x)，1x1
} scheme2_reconstruct_out;

// ========== Protocol 4  ==========

// scheme2_keygen: 
// : λ, F
// : pk_F=⊥, vk_F=⊥, {ρ_i}_{i=1}^k=F, aux_R0=⊥
void scheme2_keygen(fq_mat_t F, pubpar *par);

// scheme2_probgen: 
// : pk_F=⊥, x
// : vk_x={g^{λ'_j}}, {σ_i}, aux_R=⊥
void scheme2_probgen(scheme2_probgen_out *out, fq_mat_t x, pubpar *par);

// scheme2_compute: 
// : id i, ρ_i=F, σ_i
// : π_i = {π_{i_j} | j ∈ {0,1,...,d}}
void scheme2_compute(scheme2_compute_out *out, int server_id, fq_mat_t F, 
                       fq_mat_t *sigma_i, pubpar *par);

// scheme2_verify: 
// : vk_F=⊥, vk_x, {π_i}_{i=1}^k, v_values, F
// : aux'_R = f(0)  ⊥
int scheme2_verify(scheme2_verify_out *out, unsigned char **vk_x, 
                     scheme2_compute_out *pi_results, fq_t *v_values, 
                     fq_mat_t F, pubpar *par);

// scheme2_reconstruct: 
// : aux_R0=⊥, aux_R=⊥, aux'_R=f(0)
// : F(x)
void scheme2_reconstruct(scheme2_reconstruct_out *out, fq_mat_t aux_R_prime, pubpar *par);

// ========== Protocol 4  ==========

// λ'_j
void scheme2_compute_lambda_prime(fq_t *lambda_prime, fq_t *v_values, int d, pubpar *par);

// ：d+1f, f_{v_1}, ..., f_{v_d}
void scheme2_interpolate_polynomials(fq_poly_t *f_polys, scheme2_compute_out *pi_results, 
                                       fq_t *v_values, pubpar *par);

// ========== Protocol 4  ==========

// probgen
void scheme2_probgen_out_clear(scheme2_probgen_out *out, pubpar *par);

// compute
void scheme2_compute_out_clear(scheme2_compute_out *out, pubpar *par);

// verify
void scheme2_verify_out_clear(scheme2_verify_out *out, pubpar *par);

// reconstruct
void scheme2_reconstruct_out_clear(scheme2_reconstruct_out *out, pubpar *par);

#endif // SCHEME2_H 