#ifndef GAMMA_FE_SH_H
#define GAMMA_FE_SH_H
#include "../lib/common/crypto_util.h"
#include "../lib/common/extended_field.h"

typedef struct {
    ext_field_ctx_t field;
    unsigned char **m, **n;       /* encoded powers; safe to give to verifier */
    int count;
} gamma_fe_sh_public_key;
typedef struct {
    gamma_fe_sh_public_key vk;
    ext_field_elem_t *powers;     /* client-only decoding key */
    ext_field_mat_t *shares;
    ext_field_elem_t *masks;
} gamma_fe_sh_probgen_out;
typedef struct { ext_field_elem_t *coefficients; int valid; } gamma_fe_sh_verify_out;

int gamma_fe_sh_probgen(gamma_fe_sh_probgen_out *out, fq_mat_t input, pubpar *par);
void gamma_fe_sh_compute(ext_field_elem_t *out, fq_mat_t F, ext_field_mat_t *share,
                     const ext_field_elem_t *mask, ext_field_ctx_t *field, pubpar *par);
int gamma_fe_sh_verify(gamma_fe_sh_verify_out *out, const gamma_fe_sh_public_key *vk,
                    ext_field_elem_t *responses, const pubpar *par);
int gamma_fe_sh_decode(fq_t out, const gamma_fe_sh_probgen_out *request,
                   const gamma_fe_sh_verify_out *verified, const pubpar *par);
void gamma_fe_sh_probgen_out_clear(gamma_fe_sh_probgen_out *out, pubpar *par);
void gamma_fe_sh_verify_out_clear(gamma_fe_sh_verify_out *out, const pubpar *par);
void gamma_fe_sh_eval_extended(ext_field_elem_t *result, fq_mat_t F, ext_field_mat_t *X,
                           ext_field_ctx_t *field, pubpar *par);
#endif
