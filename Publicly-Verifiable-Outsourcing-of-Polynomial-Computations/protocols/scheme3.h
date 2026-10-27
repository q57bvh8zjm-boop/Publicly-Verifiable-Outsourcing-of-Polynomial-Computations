#ifndef SCHEME3_H
#define SCHEME3_H
#include "../lib/common/hss_paillier.h"

typedef struct { hss_public_key pk; hss_evaluation_key ek[2]; } scheme3_keygen_out;
typedef struct { unsigned char **vk_x; hss_ciphertext *ct_inputs; } scheme3_probgen_out;
typedef struct { fq_t *values; int valid; } scheme3_compute_out;
typedef struct { fq_mat_t aux_R_prime; int valid; } scheme3_verify_out;
typedef struct { fq_mat_t result; } scheme3_reconstruct_out;

int scheme3_keygen(scheme3_keygen_out *out, pubpar *par);
int scheme3_probgen(scheme3_probgen_out *out, const hss_public_key *pk,
                    fq_mat_t input, pubpar *par);
int scheme3_compute(scheme3_compute_out *out, int server,
                     const hss_evaluation_key *ek, const hss_public_key *pk,
                     fq_mat_t polynomial, hss_ciphertext *inputs, pubpar *par);
/* Only the public OV key and field output shares are used for verification. */
int scheme3_verify(scheme3_verify_out *out, unsigned char **vk,
                    scheme3_compute_out responses[2], pubpar *par);
void scheme3_reconstruct(scheme3_reconstruct_out *out, fq_mat_t eta,
                         const fq_t constant, pubpar *par);
void scheme3_keygen_out_clear(scheme3_keygen_out *out);
void scheme3_probgen_out_clear(scheme3_probgen_out *out, pubpar *par);
void scheme3_compute_out_clear(scheme3_compute_out *out, pubpar *par);
void scheme3_verify_out_clear(scheme3_verify_out *out, pubpar *par);
void scheme3_reconstruct_out_clear(scheme3_reconstruct_out *out, pubpar *par);
#endif
