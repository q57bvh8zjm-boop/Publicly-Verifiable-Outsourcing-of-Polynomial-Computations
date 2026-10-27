#ifndef PV_CRYPTO_UTIL_H
#define PV_CRYPTO_UTIL_H

#include "mpoly.h"
#include <stdint.h>

#define PV_POINT_BYTES 32
#define PV_STAT_BITS 128

void *pv_alloc(size_t count, size_t size);
void pv_random_below(fmpz_t out, const fmpz_t bound);
void pv_random_bits(fmpz_t out, unsigned bits);
void pv_random_field(fq_t out, const fq_ctx_t ctx);
void pv_group_order(fmpz_t order);
int pv_field_matches_group(const pubpar *par);
int pv_params_init(pubpar *par, int m, int d, int t);
void pv_params_clear(pubpar *par);
int pv_group_base(unsigned char out[32], const fmpz_t scalar);
int pv_group_add_scaled(unsigned char sum[32], const unsigned char point[32],
                        const fmpz_t scalar);
int pv_group_is_identity(const unsigned char point[32]);
void pv_sample_points(fq_t *v, int d, const fq_ctx_t ctx);
void pv_orthogonal_weights(fq_t *weights, fq_t *v, int d, const fq_ctx_t ctx);
int pv_ov_verify(unsigned char **vk, fq_t *values, int d, const pubpar *par);
void pv_interpolate_zero(fq_t result, fq_t *values, const pubpar *par);

#endif
