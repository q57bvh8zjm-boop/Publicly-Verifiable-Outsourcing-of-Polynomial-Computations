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


// 
typedef struct {
} scheme2_keygen_out;

//   
typedef struct {
    unsigned char **vk_x;
    
    fq_mat_t **sigma;
    
    
    fq_t *v_values;  // v_1, v_2, ..., v_d
} scheme2_probgen_out;

// 
typedef struct {
    fq_mat_t *pi;
} scheme2_compute_out;

// 
typedef struct {
    // aux'_R: f(0) if verification succeeds
    fq_mat_t aux_R_prime;
    int valid;
} scheme2_verify_out;

// 
typedef struct {
    fq_mat_t result;
} scheme2_reconstruct_out;

// ========== Protocol 4  ==========

// scheme2_keygen: 
void scheme2_keygen(fq_mat_t F, pubpar *par);

// scheme2_probgen: 
void scheme2_probgen(scheme2_probgen_out *out, fq_mat_t x, pubpar *par);

// scheme2_compute: 
void scheme2_compute(scheme2_compute_out *out, int server_id, fq_mat_t F, 
                       fq_mat_t *sigma_i, pubpar *par);

// scheme2_verify: 
int scheme2_verify(scheme2_verify_out *out, unsigned char **vk_x, 
                     scheme2_compute_out *pi_results, fq_t *v_values, 
                     fq_mat_t F, pubpar *par);

// scheme2_reconstruct: 
// : F(x)
void scheme2_reconstruct(scheme2_reconstruct_out *out, fq_mat_t aux_R_prime, pubpar *par);

// ========== Protocol 4  ==========

void scheme2_compute_lambda_prime(fq_t *lambda_prime, fq_t *v_values, int d, pubpar *par);

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