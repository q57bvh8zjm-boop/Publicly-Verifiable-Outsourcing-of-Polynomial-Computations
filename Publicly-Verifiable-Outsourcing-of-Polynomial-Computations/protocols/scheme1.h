#ifndef SCHEME1_H
#define SCHEME1_H
#include "../lib/common/crypto_util.h"
#include "../lib/common/extended_field.h"

typedef struct {
    ext_field_ctx_t field;
    unsigned char **m, **n;       /* encoded powers; safe to give to verifier */
    int count;
} scheme1_public_key;
typedef struct {
    scheme1_public_key vk;
    ext_field_elem_t *powers;     /* client-only decoding key */
    ext_field_mat_t *shares;
    ext_field_elem_t *masks;
} scheme1_probgen_out;
typedef struct { ext_field_elem_t *coefficients; int valid; } scheme1_verify_out;

int scheme1_probgen(scheme1_probgen_out *out, fq_mat_t input, pubpar *par);
void scheme1_compute(ext_field_elem_t *out, fq_mat_t F, ext_field_mat_t *share,
                     const ext_field_elem_t *mask, ext_field_ctx_t *field, pubpar *par);
int scheme1_verify(scheme1_verify_out *out, const scheme1_public_key *vk,
                    ext_field_elem_t *responses, const pubpar *par);
int scheme1_decode(fq_t out, const scheme1_probgen_out *request,
                   const scheme1_verify_out *verified, const pubpar *par);
void scheme1_probgen_out_clear(scheme1_probgen_out *out, pubpar *par);
void scheme1_verify_out_clear(scheme1_verify_out *out, const pubpar *par);
void scheme1_eval_extended(ext_field_elem_t *result, fq_mat_t F, ext_field_mat_t *X,
                           ext_field_ctx_t *field, pubpar *par);
#endif
