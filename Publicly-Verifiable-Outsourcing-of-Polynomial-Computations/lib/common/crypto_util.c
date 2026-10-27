#include "crypto_util.h"
#include <limits.h>

void *pv_alloc(size_t count, size_t size) {
    if (size && count > SIZE_MAX / size) abort();
    void *p = calloc(count ? count : 1, size ? size : 1);
    if (!p) abort();
    return p;
}

void pv_random_bits(fmpz_t out, unsigned bits) {
    if (sodium_init() < 0) abort();
    size_t n = ((size_t)bits + 7) / 8;
    unsigned char *bytes = pv_alloc(n, 1);
    randombytes_buf(bytes, n);
    if (bits % 8) bytes[n - 1] &= (1U << (bits % 8)) - 1;
    mpz_t z;
    mpz_init(z);
    mpz_import(z, n, -1, 1, 0, 0, bytes);
    fmpz_set_mpz(out, z);
    mpz_clear(z);
    sodium_memzero(bytes, n);
    free(bytes);
}

void pv_random_below(fmpz_t out, const fmpz_t bound) {
    if (fmpz_sgn(bound) <= 0) abort();
    unsigned bits = (unsigned)fmpz_bits(bound);
    do { pv_random_bits(out, bits); } while (fmpz_cmp(out, bound) >= 0);
}

void pv_random_field(fq_t out, const fq_ctx_t ctx) {
    fmpz_poly_t polynomial;
    fmpz_t value;
    fmpz_poly_init(polynomial);
    fmpz_init(value);
    for (slong i = 0; i < fq_ctx_degree(ctx); i++) {
        pv_random_below(value, fq_ctx_prime(ctx));
        fmpz_poly_set_coeff_fmpz(polynomial, i, value);
    }
    fq_set_fmpz_poly(out, polynomial, ctx);
    fmpz_clear(value);
    fmpz_poly_clear(polynomial);
}

void pv_group_order(fmpz_t order) {
    fmpz_t tail;
    fmpz_init(tail);
    fmpz_one(order);
    fmpz_mul_2exp(order, order, 252);
    fmpz_set_str(tail, "27742317777372353535851937790883648493", 10);
    fmpz_add(order, order, tail);
    fmpz_clear(tail);
}

int pv_field_matches_group(const pubpar *par) {
    fmpz_t order;
    fmpz_init(order);
    pv_group_order(order);
    int ok = fq_ctx_degree(par->Fp) == 1 && fmpz_equal(par->p, order)
             && fmpz_equal(fq_ctx_prime(par->Fp), order);
    fmpz_clear(order);
    return ok;
}

int pv_params_init(pubpar *par, int m, int d, int t) {
    if (!par || m < 1 || d < 1 || t < 1 || d > (INT_MAX - 1) / t)
        return 0;
    memset(par, 0, sizeof(*par));
    par->m = m; par->d = d; par->t = t; par->k = d * t + 1;
    fmpz_init(par->p);
    pv_group_order(par->p);
    fq_ctx_init(par->Fp, par->p, 1, "x");
    par->SID = pv_alloc((size_t)par->k, sizeof(fq_t));
    for (int i = 0; i < par->k; i++) {
        fq_init(par->SID[i], par->Fp);
        fq_set_ui(par->SID[i], (ulong)i + 1, par->Fp);
    }
    return sodium_init() >= 0;
}

void pv_params_clear(pubpar *par) {
    for (int i = 0; i < par->k; i++) fq_clear(par->SID[i], par->Fp);
    free(par->SID);
    fq_ctx_clear(par->Fp);
    fmpz_clear(par->p);
    memset(par, 0, sizeof(*par));
}

static int encode_scalar(unsigned char out[32], const fmpz_t scalar) {
    fmpz_t order;
    fmpz_init(order); pv_group_order(order);
    int ok = fmpz_sgn(scalar) >= 0 && fmpz_cmp(scalar, order) < 0;
    fmpz_clear(order);
    if (!ok) return 0;
    mpz_t z;
    mpz_init(z); fmpz_get_mpz(z, scalar);
    memset(out, 0, 32);
    mpz_export(out, NULL, -1, 1, 0, 0, z);
    mpz_clear(z);
    return 1;
}

int pv_group_base(unsigned char out[32], const fmpz_t scalar) {
    unsigned char s[32];
    if (!encode_scalar(s, scalar)) return 0;
    if (fmpz_is_zero(scalar)) { memset(out, 0, 32); return 1; }
    return crypto_scalarmult_ristretto255_base(out, s) == 0;
}

int pv_group_is_identity(const unsigned char point[32]) {
    return sodium_is_zero(point, 32);
}

int pv_group_add_scaled(unsigned char sum[32], const unsigned char point[32],
                        const fmpz_t scalar) {
    unsigned char s[32], term[32];
    if (!crypto_core_ristretto255_is_valid_point(point) ||
        !encode_scalar(s, scalar)) return 0;
    if (fmpz_is_zero(scalar) || pv_group_is_identity(point)) return 1;
    if (crypto_scalarmult_ristretto255(term, s, point) != 0) return 0;
    return crypto_core_ristretto255_add(sum, sum, term) == 0;
}

void pv_sample_points(fq_t *v, int d, const fq_ctx_t ctx) {
    fq_one(v[0], ctx);
    for (int j = 1; j <= d; j++) {
        int duplicate;
        do {
            pv_random_field(v[j], ctx);
            duplicate = fq_is_zero(v[j], ctx);
            for (int i = 0; i < j; i++) duplicate |= fq_equal(v[j], v[i], ctx);
        } while (duplicate);
    }
}

void pv_orthogonal_weights(fq_t *weights, fq_t *v, int d, const fq_ctx_t ctx) {
    fq_t common, numerator, denominator, difference;
    fq_init(common, ctx); fq_init(numerator, ctx);
    fq_init(denominator, ctx); fq_init(difference, ctx);
    fq_one(common, ctx);
    for (int i = 0; i <= d; i++)
        for (int j = i + 1; j <= d; j++) {
            fq_sub(difference, v[j], v[i], ctx);
            fq_mul(common, common, difference, ctx);
        }
    for (int j = 0; j <= d; j++) {
        fq_set(numerator, common, ctx); fq_one(denominator, ctx);
        for (int i = 0; i <= d; i++) if (i != j) {
            fq_mul(numerator, numerator, v[i], ctx);
            fq_sub(difference, v[j], v[i], ctx);
            fq_mul(denominator, denominator, difference, ctx);
        }
        if (fq_is_zero(denominator, ctx)) abort();
        fq_inv(denominator, denominator, ctx);
        fq_mul(weights[j], numerator, denominator, ctx);
    }
    fq_clear(common, ctx); fq_clear(numerator, ctx);
    fq_clear(denominator, ctx); fq_clear(difference, ctx);
}

int pv_ov_verify(unsigned char **vk, fq_t *values, int d, const pubpar *par) {
    if (!pv_field_matches_group(par) || !vk) return 0;
    unsigned char sum[32] = {0};
    fmpz_t exponent;
    fmpz_init(exponent);
    int ok = 1;
    for (int j = 0; j <= d && ok; j++) {
        fq2fmpz(exponent, values[j], par->Fp);
        ok = vk[j] && !pv_group_is_identity(vk[j]) &&
             pv_group_add_scaled(sum, vk[j], exponent);
    }
    fmpz_clear(exponent);
    return ok && pv_group_is_identity(sum);
}

void pv_interpolate_zero(fq_t result, fq_t *values, const pubpar *par) {
    fq_t weight, diff, term;
    fq_init(weight, par->Fp); fq_init(diff, par->Fp); fq_init(term, par->Fp);
    fq_zero(result, par->Fp);
    for (int i = 0; i < par->k; i++) {
        fq_one(weight, par->Fp);
        for (int j = 0; j < par->k; j++) if (i != j) {
            fq_sub(diff, par->SID[j], par->SID[i], par->Fp);
            fq_inv(diff, diff, par->Fp);
            fq_mul(diff, diff, par->SID[j], par->Fp);
            fq_mul(weight, weight, diff, par->Fp);
        }
        fq_mul(term, weight, values[i], par->Fp);
        fq_add(result, result, term, par->Fp);
    }
    fq_clear(weight, par->Fp); fq_clear(diff, par->Fp); fq_clear(term, par->Fp);
}
