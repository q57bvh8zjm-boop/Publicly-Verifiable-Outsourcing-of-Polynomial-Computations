#ifndef EXTENDED_FIELD_H
#define EXTENDED_FIELD_H

#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_mat.h>
#include <flint/fq_poly.h>
#include <flint/flint.h>

// ===========  ===========

/**
 * 
 *  a + b*ω ， a, b ∈ F_p
 */
typedef struct {
    fmpz_t a;              // 
    fmpz_t b;              // ω
} ext_field_elem_t;

/**
 * 
 * 
 */
typedef struct {
    fmpz_t p;              // F_p
    fmpz_t c;              // g(ω) = ω^2 + c
    flint_rand_t state;    // 
} ext_field_ctx_t;

// ===========  ===========

/**
 * 
 * @param ctx 
 * @param p 
 * @param c 
 */
void ext_field_ctx_init(ext_field_ctx_t *ctx, const fmpz_t p, const fmpz_t c);

/**
 * 
 * @param ctx 
 */
void ext_field_ctx_clear(ext_field_ctx_t *ctx);

/**
 * 
 * @param elem 
 */
void ext_field_elem_init(ext_field_elem_t *elem);

/**
 * 
 * @param elem 
 */
void ext_field_elem_clear(ext_field_elem_t *elem);

/**
 * 
 * @param elem 
 */
void ext_field_elem_zero(ext_field_elem_t *elem);

/**
 * 
 * @param elem 
 */
void ext_field_elem_one(ext_field_elem_t *elem);

/**
 * 
 * @param elem 
 * @param a 
 * @param b ω
 * @param ctx 
 */
void ext_field_elem_set(ext_field_elem_t *elem, const fmpz_t a, const fmpz_t b, const ext_field_ctx_t *ctx);

/**
 * 
 * @param dest 
 * @param src 
 */
void ext_field_elem_copy(ext_field_elem_t *dest, const ext_field_elem_t *src);

/**
 * 
 * @param elem 
 * @param ctx 
 */
void ext_field_elem_random(ext_field_elem_t *elem, ext_field_ctx_t *ctx);

/**
 * 
 * @param result 
 * @param a 
 * @param b 
 * @param ctx 
 */
void ext_field_elem_add(ext_field_elem_t *result, const ext_field_elem_t *a, const ext_field_elem_t *b, const ext_field_ctx_t *ctx);

/**
 * 
 * @param result 
 * @param a 
 * @param b 
 * @param ctx 
 */
void ext_field_elem_sub(ext_field_elem_t *result, const ext_field_elem_t *a, const ext_field_elem_t *b, const ext_field_ctx_t *ctx);

/**
 * 
 * @param result 
 * @param a 
 * @param b 
 * @param ctx 
 */
void ext_field_elem_mul(ext_field_elem_t *result, const ext_field_elem_t *a, const ext_field_elem_t *b, const ext_field_ctx_t *ctx);

/**
 * 
 * @param a 
 * @param b 
 * @return 1 if equal, 0 otherwise
 */
int ext_field_elem_equal(const ext_field_elem_t *a, const ext_field_elem_t *b);

/**
 * 
 * @param result 
 * @param a 
 * @param ctx 
 */
void ext_field_elem_inv(ext_field_elem_t *result, const ext_field_elem_t *a, const ext_field_ctx_t *ctx);

/**
 * 
 * @param elem 
 */
void ext_field_elem_print(const ext_field_elem_t *elem);

// ===========  ===========
typedef struct {
    ext_field_elem_t *entries;
    slong r;
    slong c;
    ext_field_elem_t **rows;
} ext_field_mat_t;

// ===========  ===========
typedef struct {
    ext_field_elem_t *coeffs;
    int degree;
} ext_field_poly_t;

void ext_field_poly_init(ext_field_poly_t *poly, int degree);
void ext_field_poly_clear(ext_field_poly_t *poly);
void ext_field_poly_set_coeff(ext_field_poly_t *poly, int d, const ext_field_elem_t *val);
void ext_field_poly_add(ext_field_poly_t *res, const ext_field_poly_t *a, const ext_field_poly_t *b, const ext_field_ctx_t *ctx);
void ext_field_poly_sub(ext_field_poly_t *res, const ext_field_poly_t *a, const ext_field_poly_t *b, const ext_field_ctx_t *ctx);
void ext_field_poly_mul(ext_field_poly_t *res, const ext_field_poly_t *a, const ext_field_poly_t *b, const ext_field_ctx_t *ctx);
void ext_field_poly_pow(ext_field_poly_t *res, const ext_field_poly_t *base, int exp, const ext_field_ctx_t *ctx);
void ext_field_poly_evaluate(ext_field_elem_t *result, const ext_field_poly_t *poly, const ext_field_elem_t *x, const ext_field_ctx_t *ctx);

// ===========  ===========
void ext_field_mat_init(ext_field_mat_t *mat, slong r, slong c);
void ext_field_mat_clear(ext_field_mat_t *mat);
ext_field_elem_t* ext_field_mat_entry(ext_field_mat_t *mat, slong i, slong j);
void ext_field_mat_entry_set(ext_field_mat_t *mat, slong i, slong j, const ext_field_elem_t *val);

// ========== Operation Counters ==========

/**
 * Get extension field multiplication count
 */
long long ext_field_get_mul_count(void);

/**
 * Get extension field addition count
 */
long long ext_field_get_add_count(void);

/**
 * Get estimated base field multiplication count
 */
long long ext_field_get_base_mul_count(void);

/**
 * Reset all operation counters
 */
void ext_field_reset_counters(void);

#endif // EXTENDED_FIELD_H