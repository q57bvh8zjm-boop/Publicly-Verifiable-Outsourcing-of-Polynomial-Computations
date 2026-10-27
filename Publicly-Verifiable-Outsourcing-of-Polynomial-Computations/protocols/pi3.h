#ifndef PV_PI3_H
#define PV_PI3_H
#include "gamma_fe_sh.h"
/* Private extension-field baseline; public-verifier keys are not generated. */
typedef struct {
    ext_field_ctx_t field;
    ext_field_elem_t alpha;
    ext_field_mat_t *shares;
    int count;
} pi3_probgen_out;
typedef struct {
    ext_field_elem_t value;
    int valid;
} pi3_verify_out;
int pi3_probgen(pi3_probgen_out *out, fq_mat_t x, pubpar *p);
void pi3_compute(ext_field_elem_t *out, fq_mat_t F, ext_field_mat_t *share, ext_field_ctx_t *field,
                 pubpar *p);
int pi3_verify(pi3_verify_out *out, const pi3_probgen_out *req, ext_field_elem_t *responses,
               const pubpar *p);
int pi3_decode(fq_t out, const pi3_verify_out *verified, const pubpar *p);
void pi3_probgen_out_clear(pi3_probgen_out *out);
void pi3_verify_out_clear(pi3_verify_out *out);
#endif
