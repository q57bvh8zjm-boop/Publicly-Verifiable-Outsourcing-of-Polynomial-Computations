/**
 * 
 * 
 */

#include "extended_field.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ========== Operation Counters ==========
static long long g_ext_mul_count = 0;
static long long g_ext_add_count = 0;
static long long g_base_mul_count = 0;

// ===========  ===========

void ext_field_ctx_init(ext_field_ctx_t *ctx, const fmpz_t p, const fmpz_t c) {
    fmpz_init(ctx->p);
    fmpz_init(ctx->c);
    fmpz_set(ctx->p, p);
    fmpz_set(ctx->c, c);
    flint_randinit(ctx->state);
}

void ext_field_ctx_clear(ext_field_ctx_t *ctx) {
    fmpz_clear(ctx->p);
    fmpz_clear(ctx->c);
    flint_randclear(ctx->state);
}

// ===========  ===========

void ext_field_elem_init(ext_field_elem_t *elem) {
    fmpz_init(elem->a);
    fmpz_init(elem->b);
}

void ext_field_elem_clear(ext_field_elem_t *elem) {
    fmpz_clear(elem->a);
    fmpz_clear(elem->b);
}

void ext_field_elem_zero(ext_field_elem_t *elem) {
    fmpz_zero(elem->a);
    fmpz_zero(elem->b);
}

void ext_field_elem_one(ext_field_elem_t *elem) {
    fmpz_one(elem->a);
    fmpz_zero(elem->b);
}

void ext_field_elem_set(ext_field_elem_t *elem, const fmpz_t a, const fmpz_t b, const ext_field_ctx_t *ctx) {
    fmpz_mod(elem->a, a, ctx->p);
    fmpz_mod(elem->b, b, ctx->p);
}

void ext_field_elem_copy(ext_field_elem_t *dest, const ext_field_elem_t *src) {
    fmpz_set(dest->a, src->a);
    fmpz_set(dest->b, src->b);
}

void ext_field_elem_random(ext_field_elem_t *elem, ext_field_ctx_t *ctx) {
    fmpz_randm(elem->a, ctx->state, ctx->p);
    fmpz_randm(elem->b, ctx->state, ctx->p);
}

// ===========  ===========

void ext_field_elem_add(ext_field_elem_t *result, const ext_field_elem_t *a, const ext_field_elem_t *b, const ext_field_ctx_t *ctx) {
    g_ext_add_count++;
    fmpz_add(result->a, a->a, b->a);
    fmpz_add(result->b, a->b, b->b);
    
    // p
    fmpz_mod(result->a, result->a, ctx->p);
    fmpz_mod(result->b, result->b, ctx->p);
}

void ext_field_elem_sub(ext_field_elem_t *result, const ext_field_elem_t *a, const ext_field_elem_t *b, const ext_field_ctx_t *ctx) {
    fmpz_sub(result->a, a->a, b->a);
    fmpz_sub(result->b, a->b, b->b);
    
    // p
    fmpz_mod(result->a, result->a, ctx->p);
    fmpz_mod(result->b, result->b, ctx->p);
}

void ext_field_elem_mul(ext_field_elem_t *result, const ext_field_elem_t *a, const ext_field_elem_t *b, const ext_field_ctx_t *ctx) {
    
    g_ext_mul_count++;
    g_base_mul_count += 5; // 5 base field multiplications
    
    fmpz_t temp1, temp2, temp3, temp4;
    fmpz_init(temp1);
    fmpz_init(temp2);
    fmpz_init(temp3);
    fmpz_init(temp4);
    
    // : a1*a2 - b1*b2*c
    fmpz_mul(temp1, a->a, b->a);    // a1*a2
    fmpz_mul(temp2, a->b, b->b);    // b1*b2
    fmpz_mul(temp2, temp2, ctx->c); // b1*b2*c
    fmpz_sub(result->a, temp1, temp2); // a1*a2 - b1*b2*c
    
    fmpz_mul(temp3, a->a, b->b);    // a1*b2
    fmpz_mul(temp4, a->b, b->a);    // a2*b1
    fmpz_add(result->b, temp3, temp4); // a1*b2 + a2*b1
    
    // p
    fmpz_mod(result->a, result->a, ctx->p);
    fmpz_mod(result->b, result->b, ctx->p);
    
    fmpz_clear(temp1);
    fmpz_clear(temp2);
    fmpz_clear(temp3);
    fmpz_clear(temp4);
}

int ext_field_elem_equal(const ext_field_elem_t *a, const ext_field_elem_t *b) {
    return fmpz_equal(a->a, b->a) && fmpz_equal(a->b, b->b);
}

void ext_field_elem_inv(ext_field_elem_t *result, const ext_field_elem_t *a, const ext_field_ctx_t *ctx) {
    
    fmpz_t norm, temp1, temp2;
    fmpz_init(norm);
    fmpz_init(temp1);
    fmpz_init(temp2);
    
    //  norm = a^2 + b^2*c
    fmpz_mul(temp1, a->a, a->a);    // a^2
    fmpz_mul(temp2, a->b, a->b);    // b^2
    fmpz_mul(temp2, temp2, ctx->c); // b^2*c
    fmpz_add(norm, temp1, temp2);   // a^2 + b^2*c
    
    //  norm 
    fmpz_invmod(norm, norm, ctx->p);
    
    // 
    fmpz_mul(result->a, a->a, norm);  // a / norm
    fmpz_neg(temp1, a->b);            // -b
    fmpz_mul(result->b, temp1, norm); // -b / norm
    
    // p
    fmpz_mod(result->a, result->a, ctx->p);
    fmpz_mod(result->b, result->b, ctx->p);
    
    fmpz_clear(norm);
    fmpz_clear(temp1);
    fmpz_clear(temp2);
}

// ========== Operation Counter Functions ==========

long long ext_field_get_mul_count(void) {
    return g_ext_mul_count;
}

long long ext_field_get_add_count(void) {
    return g_ext_add_count;
}

long long ext_field_get_base_mul_count(void) {
    return g_base_mul_count;
}

void ext_field_reset_counters(void) {
    g_ext_mul_count = 0;
    g_ext_add_count = 0;
    g_base_mul_count = 0;
}

void ext_field_elem_print(const ext_field_elem_t *elem) {
    printf("(");
    fmpz_print(elem->a);
    printf(" + ");
    fmpz_print(elem->b);
    printf("*w)");
} 

// ===========  ===========
void ext_field_poly_init(ext_field_poly_t *poly, int degree) {
    poly->degree = degree;
    poly->coeffs = (ext_field_elem_t*)malloc((degree + 1) * sizeof(ext_field_elem_t));
    for (int i = 0; i <= degree; i++) {
        ext_field_elem_init(&poly->coeffs[i]);
        ext_field_elem_zero(&poly->coeffs[i]);
    }
}

void ext_field_poly_clear(ext_field_poly_t *poly) {
    for (int i = 0; i <= poly->degree; i++) {
        ext_field_elem_clear(&poly->coeffs[i]);
    }
    free(poly->coeffs);
}

void ext_field_poly_set_coeff(ext_field_poly_t *poly, int d, const ext_field_elem_t *val) {
    if (d < 0 || d > poly->degree) return;
    ext_field_elem_copy(&poly->coeffs[d], val);
}

void ext_field_poly_add(ext_field_poly_t *res, const ext_field_poly_t *a, const ext_field_poly_t *b, const ext_field_ctx_t *ctx) {
    int deg = (a->degree > b->degree) ? a->degree : b->degree;
    ext_field_poly_init(res, deg);
    for (int i = 0; i <= deg; i++) {
        if (i <= a->degree && i <= b->degree)
            ext_field_elem_add(&res->coeffs[i], &a->coeffs[i], &b->coeffs[i], ctx);
        else if (i <= a->degree)
            ext_field_elem_copy(&res->coeffs[i], &a->coeffs[i]);
        else if (i <= b->degree)
            ext_field_elem_copy(&res->coeffs[i], &b->coeffs[i]);
    }
}

void ext_field_poly_sub(ext_field_poly_t *res, const ext_field_poly_t *a, const ext_field_poly_t *b, const ext_field_ctx_t *ctx) {
    int deg = (a->degree > b->degree) ? a->degree : b->degree;
    ext_field_poly_init(res, deg);
    for (int i = 0; i <= deg; i++) {
        if (i <= a->degree && i <= b->degree)
            ext_field_elem_sub(&res->coeffs[i], &a->coeffs[i], &b->coeffs[i], ctx);
        else if (i <= a->degree)
            ext_field_elem_copy(&res->coeffs[i], &a->coeffs[i]);
        else if (i <= b->degree) {
            ext_field_elem_t zero; ext_field_elem_init(&zero); ext_field_elem_zero(&zero);
            ext_field_elem_sub(&res->coeffs[i], &zero, &b->coeffs[i], ctx);
            ext_field_elem_clear(&zero);
        }
    }
}

void ext_field_poly_mul(ext_field_poly_t *res, const ext_field_poly_t *a, const ext_field_poly_t *b, const ext_field_ctx_t *ctx) {
    int deg = a->degree + b->degree;
    ext_field_poly_init(res, deg);
    for (int i = 0; i <= a->degree; i++) {
        for (int j = 0; j <= b->degree; j++) {
            ext_field_elem_t temp;
            ext_field_elem_init(&temp);
            ext_field_elem_mul(&temp, &a->coeffs[i], &b->coeffs[j], ctx);
            ext_field_elem_add(&res->coeffs[i+j], &res->coeffs[i+j], &temp, ctx);
            ext_field_elem_clear(&temp);
        }
    }
}

void ext_field_poly_pow(ext_field_poly_t *res, const ext_field_poly_t *base, int exp, const ext_field_ctx_t *ctx) {
    ext_field_poly_t tmp, acc;
    ext_field_poly_init(&acc, 0);
    ext_field_elem_one(&acc.coeffs[0]);
    if (exp == 0) {
        ext_field_poly_init(res, 0);
        ext_field_elem_one(&res->coeffs[0]);
        ext_field_poly_clear(&acc);
        return;
    }
    ext_field_poly_init(&tmp, base->degree);
    for (int i = 0; i <= base->degree; i++)
        ext_field_elem_copy(&tmp.coeffs[i], &base->coeffs[i]);
    for (int e = 0; e < exp; e++) {
        ext_field_poly_t prod;
        ext_field_poly_mul(&prod, &acc, &tmp, ctx);
        ext_field_poly_clear(&acc);
        acc = prod;
    }
    ext_field_poly_init(res, acc.degree);
    for (int i = 0; i <= acc.degree; i++)
        ext_field_elem_copy(&res->coeffs[i], &acc.coeffs[i]);
    ext_field_poly_clear(&acc);
    ext_field_poly_clear(&tmp);
}

void ext_field_poly_evaluate(ext_field_elem_t *result, const ext_field_poly_t *poly, const ext_field_elem_t *x, const ext_field_ctx_t *ctx) {
    ext_field_elem_zero(result);
    ext_field_elem_t xpow, temp;
    ext_field_elem_init(&xpow);
    ext_field_elem_init(&temp);
    ext_field_elem_one(&xpow);
    for (int i = 0; i <= poly->degree; i++) {
        ext_field_elem_mul(&temp, &poly->coeffs[i], &xpow, ctx);
        ext_field_elem_add(result, result, &temp, ctx);
        if (i < poly->degree)
            ext_field_elem_mul(&xpow, &xpow, x, ctx);
    }
    ext_field_elem_clear(&xpow);
    ext_field_elem_clear(&temp);
}

// ===========  ===========

void ext_field_mat_init(ext_field_mat_t *mat, slong r, slong c) {
    mat->r = r;
    mat->c = c;
    mat->entries = (ext_field_elem_t*)malloc(r * c * sizeof(ext_field_elem_t));
    mat->rows = (ext_field_elem_t**)malloc(r * sizeof(ext_field_elem_t*));
    
    for (slong i = 0; i < r; i++) {
        mat->rows[i] = mat->entries + i * c;
        for (slong j = 0; j < c; j++) {
            ext_field_elem_init(&mat->rows[i][j]);
            ext_field_elem_zero(&mat->rows[i][j]);
        }
    }
}

void ext_field_mat_clear(ext_field_mat_t *mat) {
    for (slong i = 0; i < mat->r; i++) {
        for (slong j = 0; j < mat->c; j++) {
            ext_field_elem_clear(&mat->rows[i][j]);
        }
    }
    free(mat->entries);
    free(mat->rows);
}

ext_field_elem_t* ext_field_mat_entry(ext_field_mat_t *mat, slong i, slong j) {
    if (i < 0 || i >= mat->r || j < 0 || j >= mat->c) return NULL;
    return &mat->rows[i][j];
}

void ext_field_mat_entry_set(ext_field_mat_t *mat, slong i, slong j, const ext_field_elem_t *val) {
    if (i < 0 || i >= mat->r || j < 0 || j >= mat->c) return;
    ext_field_elem_copy(&mat->rows[i][j], val);
} 
