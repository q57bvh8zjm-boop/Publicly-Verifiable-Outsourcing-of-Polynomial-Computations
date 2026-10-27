#ifndef SCHEME2_H
#define SCHEME2_H
#include "../lib/common/crypto_util.h"

/* Each algorithm initializes its output; release it with the matching clear. */
typedef struct {
    unsigned char **vk_x;
    fq_mat_t **sigma;             /* sigma[server][column], each m x 1 */
    fq_t *v_values;               /* client-only scaling points v_1,...,v_d */
} scheme2_probgen_out;
typedef struct { fq_mat_t *pi; } scheme2_compute_out;
typedef struct { fq_mat_t aux_R_prime; int valid; } scheme2_verify_out;
typedef struct { fq_mat_t result; } scheme2_reconstruct_out;

void scheme2_keygen(fq_mat_t F, pubpar *par);
int scheme2_probgen(scheme2_probgen_out *out, fq_mat_t x, pubpar *par);
void scheme2_compute(scheme2_compute_out *out, int server_id, fq_mat_t F,
                     fq_mat_t *sigma_i, pubpar *par);
/* Public verification takes neither the input nor the secret scaling points. */
int scheme2_verify(scheme2_verify_out *out, unsigned char **vk_x,
                   scheme2_compute_out *responses, pubpar *par);
void scheme2_reconstruct(scheme2_reconstruct_out *out, fq_mat_t eta,
                         const fq_t constant, pubpar *par);
void scheme2_probgen_out_clear(scheme2_probgen_out *out, pubpar *par);
void scheme2_compute_out_clear(scheme2_compute_out *out, pubpar *par);
void scheme2_verify_out_clear(scheme2_verify_out *out, pubpar *par);
void scheme2_reconstruct_out_clear(scheme2_reconstruct_out *out, pubpar *par);
#endif
