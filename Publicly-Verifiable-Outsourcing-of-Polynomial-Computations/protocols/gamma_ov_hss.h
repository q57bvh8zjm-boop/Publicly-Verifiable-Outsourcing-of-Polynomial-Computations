#ifndef GAMMA_OV_HSS_H
#define GAMMA_OV_HSS_H
#include "../lib/common/hss_paillier.h"

typedef struct { hss_public_key pk; hss_evaluation_key ek[2]; } gamma_ov_hss_keygen_out;
typedef struct { unsigned char **vk_x; hss_ciphertext *ct_inputs; } gamma_ov_hss_probgen_out;
typedef struct { fq_t *values; int valid; } gamma_ov_hss_compute_out;
typedef struct { fq_mat_t aux_R_prime; int valid; } gamma_ov_hss_verify_out;
typedef struct { fq_mat_t result; } gamma_ov_hss_reconstruct_out;

int gamma_ov_hss_keygen(gamma_ov_hss_keygen_out *out, pubpar *par);
int gamma_ov_hss_probgen(gamma_ov_hss_probgen_out *out, const hss_public_key *pk,
                    fq_mat_t input, pubpar *par);
int gamma_ov_hss_compute(gamma_ov_hss_compute_out *out, int server,
                     const hss_evaluation_key *ek, const hss_public_key *pk,
                     fq_mat_t polynomial, hss_ciphertext *inputs, pubpar *par);
/* Only the public OV key and field output shares are used for verification. */
int gamma_ov_hss_verify(gamma_ov_hss_verify_out *out, unsigned char **vk,
                    gamma_ov_hss_compute_out responses[2], pubpar *par);
void gamma_ov_hss_reconstruct(gamma_ov_hss_reconstruct_out *out, fq_mat_t eta,
                         const fq_t constant, pubpar *par);
void gamma_ov_hss_keygen_out_clear(gamma_ov_hss_keygen_out *out);
void gamma_ov_hss_probgen_out_clear(gamma_ov_hss_probgen_out *out, pubpar *par);
void gamma_ov_hss_compute_out_clear(gamma_ov_hss_compute_out *out, pubpar *par);
void gamma_ov_hss_verify_out_clear(gamma_ov_hss_verify_out *out, pubpar *par);
void gamma_ov_hss_reconstruct_out_clear(gamma_ov_hss_reconstruct_out *out, pubpar *par);
#endif
