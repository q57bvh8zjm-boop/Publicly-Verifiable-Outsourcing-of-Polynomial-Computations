#ifndef GAMMA_OV_SH_H
#define GAMMA_OV_SH_H
#include "../lib/common/crypto_util.h"

/* Each algorithm initializes its output; release it with the matching clear. */
typedef struct {
    unsigned char **vk_x;
    fq_mat_t **sigma;             /* sigma[server][column], each m x 1 */
    fq_t *v_values;               /* client-only scaling points v_1,...,v_d */
} gamma_ov_sh_probgen_out;
typedef struct { fq_mat_t *pi; } gamma_ov_sh_compute_out;
typedef struct { fq_mat_t aux_R_prime; int valid; } gamma_ov_sh_verify_out;
typedef struct { fq_mat_t result; } gamma_ov_sh_reconstruct_out;

void gamma_ov_sh_keygen(fq_mat_t F, pubpar *par);
int gamma_ov_sh_probgen(gamma_ov_sh_probgen_out *out, fq_mat_t x, pubpar *par);
void gamma_ov_sh_compute(gamma_ov_sh_compute_out *out, int server_id, fq_mat_t F,
                     fq_mat_t *sigma_i, pubpar *par);
/* Public verification takes neither the input nor the secret scaling points. */
int gamma_ov_sh_verify(gamma_ov_sh_verify_out *out, unsigned char **vk_x,
                   gamma_ov_sh_compute_out *responses, pubpar *par);
void gamma_ov_sh_reconstruct(gamma_ov_sh_reconstruct_out *out, fq_mat_t eta,
                         const fq_t constant, pubpar *par);
void gamma_ov_sh_probgen_out_clear(gamma_ov_sh_probgen_out *out, pubpar *par);
void gamma_ov_sh_compute_out_clear(gamma_ov_sh_compute_out *out, pubpar *par);
void gamma_ov_sh_verify_out_clear(gamma_ov_sh_verify_out *out, pubpar *par);
void gamma_ov_sh_reconstruct_out_clear(gamma_ov_sh_reconstruct_out *out, pubpar *par);
#endif
