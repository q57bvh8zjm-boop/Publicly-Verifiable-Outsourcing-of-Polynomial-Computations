#ifndef SCHEME4_H
#define SCHEME4_H

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

// ========== Protocol 8  ==========

// CNF 
//  i，F_i = {F_{i,1}, ..., F_{i,L}} 
typedef struct {
    fq_mat_t *F_i;  // ， L
    int L;          // 
} cnf_function_share_t;

// 
typedef struct {
    // pk_F = ⊥
    // vk_F = ⊥
    // rho_i = F_i (CNF )
    cnf_function_share_t *rho;  // ， k
    // aux_R0 = ⊥
} scheme4_keygen_out;

// 
typedef struct {
    // vk_x = g^{α^d} ()
    unsigned char *vk_x;  // 32
    
    // sigma_i = {s_i, tilde{s}_i} ()
    // s_i = (s_i^{(1)}, ..., s_i^{(m)})
    // s_i^{(j)} = {x_{j,u} | u_i = 0}
    // tilde{s}_i = (tilde{s}_i^{(1)}, ..., tilde{s}_i^{(m)})
    // tilde{s}_i^{(j)} = {r_{j,u} | u_i = 0}
    fq_t ***s_i;        // s_i[j][u_idx]  i， j，CNF 
    fq_t ***tilde_s_i;  // tilde{s}_i[j][u_idx]
    int *num_shares_per_coord;  //  CNF 
    
    // aux_D = x ( Decode)
    fq_mat_t aux_D;
    
    // ：α （）
    fq_t alpha;
    
    //  r = αx（）
    fq_mat_t r;
} scheme4_probgen_out;

// 
typedef struct {
    // Y_i = (y_{i,1}, ..., y_{i,L})  y_{i,ell'} = F_{i,ell'}(s_i)
    fq_t *Y_i;
    
    // tilde{Y}_i = (tilde{y}_{i,1}, ..., tilde{y}_{i,L})  tilde{y}_{i,ell'} = F_{i,ell'}(tilde{s}_i)
    fq_t *tilde_Y_i;
    
    // r_i^{(j)} = {g^{x_{j,u}} | u_i = 0}
    // tilde{r}_i^{(j)} = {g^{r_{j,u}} | u_i = 0}
    unsigned char ***r_i;      // r_i[j][u_idx] 
    unsigned char ***tilde_r_i; // tilde{r}_i[j][u_idx] 
    
    int L;  // 
} scheme4_compute_out;

// 
typedef struct {
    // aux'_D = {Y_i | i ∈ [k]} if verification succeeds
    fq_t **aux_D_prime;  // Y_i ， k
    int *aux_D_prime_L;  //  Y_i ， k
    int valid;           // ：1，0
} scheme4_verify_out;

// 
typedef struct {
    fq_t result;  // F(x) 
} scheme4_decode_out;

// ========== Protocol 8  ==========

// scheme4_keygen: 
// : λ, F
// : pk_F=⊥, vk_F=⊥, {ρ_i}_{i=1}^k=F_i, aux_R0=⊥
void scheme4_keygen(scheme4_keygen_out *out, fq_mat_t F, pubpar *par);

// scheme4_probgen: 
// : pk_F=⊥, x
// : vk_x=g^{α^d}, {σ_i}, aux_D=x
void scheme4_probgen(scheme4_probgen_out *out, fq_mat_t x, pubpar *par);

// scheme4_compute: 
// : id i, ρ_i=F_i, σ_i={s_i, tilde{s}_i}
// : π_i = {Y_i, tilde{Y}_i, {r_i^{(j)}, tilde{r}_i^{(j)} | j ∈ [m]}}
void scheme4_compute(scheme4_compute_out *out, int server_id, 
                      cnf_function_share_t *rho_i, scheme4_probgen_out *probgen_out, 
                      pubpar *par);

// scheme4_verify: 
// : vk_F=⊥, vk_x=g^{α^d}, {π_i}_{i=1}^k, probgen_out
// : aux'_D = {Y_i | i ∈ [k]}  ⊥
int scheme4_verify(scheme4_verify_out *out, unsigned char *vk_x,
                    scheme4_compute_out *pi_results, int k,
                    scheme4_probgen_out *probgen_out, pubpar *par);

// scheme4_decode: 
// : aux_D=x, aux'_D={Y_i | i ∈ [k]}
// : F(x)
void scheme4_decode(scheme4_decode_out *out, fq_mat_t aux_D, 
                     fq_t **aux_D_prime, scheme4_keygen_out *keygen_out, 
                     pubpar *par);

// ========== Protocol 8  ==========

// CNF secret sharing:  F  F_i
// ， CNF 
void scheme4_cnf_split(cnf_function_share_t *F_shares, fq_mat_t F, pubpar *par);

// CNF reconstruction: （ Verify，l=1 ）
// ， g^
void scheme4_cnf_reconstruct_verify(fq_t *result, unsigned char **g_coeffs, 
                                      int num_shares, pubpar *par);

// CNF reconstruction: （ Decode，）
void scheme4_cnf_reconstruct_decode(fq_t *result, fq_t *shares, 
                                     int num_shares, pubpar *par);

//  H(u, 0) = t
int scheme4_hamming_weight(int u, int k);

//  H(u, 0) = t  u ∈ {0,1}^k
int scheme4_generate_cnf_vectors(int **u_vectors, int k, int t);

// ========== Protocol 8  ==========

void scheme4_keygen_out_clear(scheme4_keygen_out *out, pubpar *par);
void scheme4_probgen_out_clear(scheme4_probgen_out *out, pubpar *par);
void scheme4_compute_out_clear(scheme4_compute_out *out, pubpar *par);
void scheme4_verify_out_clear(scheme4_verify_out *out, int k, pubpar *par);
void scheme4_decode_out_clear(scheme4_decode_out *out, pubpar *par);

#endif // SCHEME4_H
