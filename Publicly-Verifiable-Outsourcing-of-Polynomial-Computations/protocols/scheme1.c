#include "scheme1.h"

static void extension_init(ext_field_ctx_t *field, const fmpz_t prime) {
    fmpz_t c, minus_c;
    fmpz_init(c); fmpz_init(minus_c); fmpz_one(c);
    do {
        fmpz_neg(minus_c, c);
        if (fmpz_jacobi(minus_c, prime) == -1) break;
        fmpz_add_ui(c, c, 1);
    } while (1);
    ext_field_ctx_init(field, prime, c);
    fmpz_clear(c); fmpz_clear(minus_c);
}

/* Lagrange interpolation at 1,...,k; all basis coefficients are in F_q.
 * Build the product once and use synthetic division: O(k^2) field operations.
 */
static void interpolate(ext_field_elem_t *coefficients, ext_field_elem_t *values,
                        int k, const ext_field_ctx_t *field) {
    fmpz_t *product = pv_alloc((size_t)k + 1, sizeof(fmpz_t));
    fmpz_t *basis = pv_alloc((size_t)k, sizeof(fmpz_t));
    for (int r = 0; r <= k; r++) fmpz_init(product[r]);
    for (int r = 0; r < k; r++) { fmpz_init(basis[r]); ext_field_elem_zero(&coefficients[r]); }
    fmpz_one(product[0]);
    fmpz_t temp, denominator, weight;
    fmpz_init(temp); fmpz_init(denominator); fmpz_init(weight);
    for (int j = 1; j <= k; j++) {
        for (int r = j; r > 0; r--) {
            fmpz_mul_ui(temp, product[r], (ulong)j);
            fmpz_sub(product[r], product[r-1], temp);
            fmpz_mod(product[r], product[r], field->p);
        }
        fmpz_mul_si(product[0], product[0], -j);
        fmpz_mod(product[0], product[0], field->p);
    }
    for (int j = 1; j <= k; j++) {
        fmpz_one(basis[k-1]);
        for (int r = k-2; r >= 0; r--) {
            fmpz_mul_ui(temp, basis[r+1], (ulong)j);
            fmpz_add(basis[r], product[r+1], temp);
            fmpz_mod(basis[r], basis[r], field->p);
        }
        fmpz_one(denominator);
        for (int i = 1; i <= k; i++) if (i != j)
            fmpz_mul_si(denominator, denominator, j-i);
        fmpz_mod(denominator, denominator, field->p);
        if (!fmpz_invmod(denominator, denominator, field->p)) abort();
        for (int r = 0; r < k; r++) {
            fmpz_mul(weight, basis[r], denominator);
            fmpz_mod(weight, weight, field->p);
            fmpz_addmul(coefficients[r].a, values[j-1].a, weight);
            fmpz_addmul(coefficients[r].b, values[j-1].b, weight);
            fmpz_mod(coefficients[r].a, coefficients[r].a, field->p);
            fmpz_mod(coefficients[r].b, coefficients[r].b, field->p);
        }
    }
    for (int r = 0; r <= k; r++) fmpz_clear(product[r]);
    for (int r = 0; r < k; r++) fmpz_clear(basis[r]);
    free(product); free(basis);
    fmpz_clear(temp); fmpz_clear(denominator); fmpz_clear(weight);
}
typedef struct {
    ext_field_elem_t *result;
    fq_mat_struct *F;
    ext_field_mat_t *X;
    ext_field_ctx_t *ctx;
    pubpar *par;
    slong next_coeff;
    slong coeff_count;
} scheme1_eval_state_t;

static void scheme1_eval_degree(scheme1_eval_state_t *state, int start_var,
                                int remaining,
                                const ext_field_elem_t *monomial) {
    if (state->next_coeff >= state->coeff_count)
        return;

    if (remaining == 0) {
        fq_t coeff_fq;
        fmpz_t coeff_int, zero;
        ext_field_elem_t coeff, term;
        fq_init(coeff_fq, state->par->Fp);
        fmpz_init(coeff_int);
        fmpz_init(zero);
        ext_field_elem_init(&coeff);
        ext_field_elem_init(&term);
        fmpz_zero(zero);

        fq_set(coeff_fq,
               fq_mat_entry(state->F, state->next_coeff, 0),
               state->par->Fp);
        fq2fmpz(coeff_int, coeff_fq, state->par->Fp);
        ext_field_elem_set(&coeff, coeff_int, zero, state->ctx);
        ext_field_elem_mul(&term, &coeff, monomial, state->ctx);
        ext_field_elem_add(state->result, state->result, &term, state->ctx);
        state->next_coeff++;

        ext_field_elem_clear(&term);
        ext_field_elem_clear(&coeff);
        fmpz_clear(zero);
        fmpz_clear(coeff_int);
        fq_clear(coeff_fq, state->par->Fp);
        return;
    }

    for (int i = start_var; i < state->par->m; i++) {
        ext_field_elem_t next_monomial;
        ext_field_elem_init(&next_monomial);
        ext_field_elem_mul(&next_monomial, monomial,
                           ext_field_mat_entry(state->X, i, 0), state->ctx);
        scheme1_eval_degree(state, i, remaining - 1, &next_monomial);
        ext_field_elem_clear(&next_monomial);
        if (state->next_coeff >= state->coeff_count)
            return;
    }
}

void scheme1_eval_extended(ext_field_elem_t *result, fq_mat_t F,
                           ext_field_mat_t *X_ext, ext_field_ctx_t *ctx,
                           pubpar *par) {
    slong rows = fq_mat_nrows(F, par->Fp);
    ext_field_elem_zero(result);
    if (rows == 0 || par->d < 0)
        return;

    fq_t constant_fq;
    fmpz_t constant, zero;
    fq_init(constant_fq, par->Fp);
    fmpz_init(constant);
    fmpz_init(zero);
    fmpz_zero(zero);
    fq_set(constant_fq, fq_mat_entry(F, 0, 0), par->Fp);
    fq2fmpz(constant, constant_fq, par->Fp);
    ext_field_elem_set(result, constant, zero, ctx);

    ext_field_elem_t one;
    ext_field_elem_init(&one);
    ext_field_elem_one(&one);
    scheme1_eval_state_t state = { result, F, X_ext, ctx, par, 1, rows };
    for (int degree = 1; degree <= par->d && state.next_coeff < rows; degree++)
        scheme1_eval_degree(&state, 0, degree, &one);

    ext_field_elem_clear(&one);
    fmpz_clear(zero);
    fmpz_clear(constant);
    fq_clear(constant_fq, par->Fp);
}



int scheme1_probgen(scheme1_probgen_out *out, fq_mat_t input, pubpar *par) {
    memset(out, 0, sizeof(*out));
    int k = par->k, t = par->t, m = par->m;
    if (!pv_field_matches_group(par) || par->d < 1 || t < 1 ||
        (long long)par->d * t + 1 != k || fq_mat_nrows(input, par->Fp) != m ||
        fq_mat_ncols(input, par->Fp) != 1) return 0;
    extension_init(&out->vk.field, par->p); out->vk.count = k;
    ext_field_ctx_t *field = &out->vk.field;
    ext_field_elem_t alpha, u, difference, product, term;
    ext_field_elem_init(&alpha); ext_field_elem_init(&u);
    ext_field_elem_init(&difference); ext_field_elem_init(&product); ext_field_elem_init(&term);
    do { ext_field_elem_random(&alpha, field); }
    while (fmpz_is_zero(alpha.b) && fmpz_cmp_ui(alpha.a, (ulong)k) <= 0);
    out->powers = pv_alloc((size_t)k, sizeof(*out->powers));
    out->vk.m = pv_alloc((size_t)k, sizeof(*out->vk.m));
    out->vk.n = pv_alloc((size_t)k, sizeof(*out->vk.n));
    for (int i = 0; i < k; i++) {
        ext_field_elem_init(&out->powers[i]);
        if (!i) ext_field_elem_one(&out->powers[i]);
        else ext_field_elem_mul(&out->powers[i], &out->powers[i-1], &alpha, field);
        out->vk.m[i] = pv_alloc(32, 1); out->vk.n[i] = pv_alloc(32, 1);
        if (!pv_group_base(out->vk.m[i], out->powers[i].a) ||
            !pv_group_base(out->vk.n[i], out->powers[i].b)) abort();
    }
    out->shares = pv_alloc((size_t)k, sizeof(*out->shares));
    out->masks = pv_alloc((size_t)k, sizeof(*out->masks));
    for (int i = 0; i < k; i++) {
        ext_field_mat_init(&out->shares[i], m, 1);
        ext_field_elem_init(&out->masks[i]);
    }
    ext_field_elem_t *random = pv_alloc((size_t)t, sizeof(*random));
    for (int s = 0; s < t; s++) ext_field_elem_init(&random[s]);
    for (int row = 0; row < m; row++) {
        for (int s = 0; s < t; s++) ext_field_elem_random(&random[s], field);
        for (int i = 0; i < k; i++) {
            ext_field_elem_t *share = ext_field_mat_entry(&out->shares[i], row, 0);
            fq2fmpz(share->a, fq_mat_entry(input, row, 0), par->Fp);
            fmpz_set_ui(u.a, (ulong)i+1); fmpz_zero(u.b);
            ext_field_elem_one(&product);
            for (int s = 1; s <= t; s++) {
                ext_field_elem_mul(&product, &product, &u, field);
                ext_field_elem_sub(&difference, &product, &out->powers[s], field);
                ext_field_elem_mul(&term, &random[s-1], &difference, field);
                ext_field_elem_add(share, share, &term, field);
            }
        }
    }
    for (int s = 0; s < t; s++) ext_field_elem_clear(&random[s]);
    free(random);
    ext_field_elem_t *h = pv_alloc((size_t)k, sizeof(*h));
    for (int s = 0; s < k; s++) {
        ext_field_elem_init(&h[s]);
        if (s) ext_field_elem_random(&h[s], field);
    }
    for (int i = 0; i < k; i++) {
        fmpz_set_ui(u.a, (ulong)i+1); fmpz_zero(u.b);
        ext_field_elem_sub(&difference, &u, &alpha, field);
        ext_field_elem_one(&product); ext_field_elem_zero(&out->masks[i]);
        for (int s = 1; s < k; s++) {
            ext_field_elem_mul(&product, &product, &difference, field);
            ext_field_elem_mul(&term, &h[s], &product, field);
            ext_field_elem_add(&out->masks[i], &out->masks[i], &term, field);
        }
    }
    for (int s = 0; s < k; s++) ext_field_elem_clear(&h[s]);
    free(h);
    ext_field_elem_clear(&alpha); ext_field_elem_clear(&u);
    ext_field_elem_clear(&difference); ext_field_elem_clear(&product); ext_field_elem_clear(&term);
    return 1;
}

void scheme1_compute(ext_field_elem_t *out, fq_mat_t F, ext_field_mat_t *share,
                     const ext_field_elem_t *mask, ext_field_ctx_t *field, pubpar *par) {
    scheme1_eval_extended(out, F, share, field, par);
    ext_field_elem_add(out, out, mask, field);
}

int scheme1_verify(scheme1_verify_out *out, const scheme1_public_key *vk,
                    ext_field_elem_t *responses, const pubpar *par) {
    memset(out, 0, sizeof(*out));
    if (!pv_field_matches_group(par) || !vk || vk->count != par->k || !responses ||
        !fmpz_equal(vk->field.p, par->p)) return 0;
    out->coefficients = pv_alloc((size_t)par->k, sizeof(*out->coefficients));
    for (int i = 0; i < par->k; i++) ext_field_elem_init(&out->coefficients[i]);
    interpolate(out->coefficients, responses, par->k, &vk->field);
    unsigned char sum[32] = {0};
    int ok = 1;
    for (int i = 0; i < par->k && ok; i++)
        ok = pv_group_add_scaled(sum, vk->n[i], out->coefficients[i].a) &&
             pv_group_add_scaled(sum, vk->m[i], out->coefficients[i].b);
    out->valid = ok && pv_group_is_identity(sum);
    return out->valid;
}

int scheme1_decode(fq_t out, const scheme1_probgen_out *request,
                   const scheme1_verify_out *verified, const pubpar *par) {
    if (!verified->valid) return 0;
    ext_field_elem_t sum, term;
    ext_field_elem_init(&sum); ext_field_elem_init(&term); ext_field_elem_zero(&sum);
    for (int i = 0; i < par->k; i++) {
        ext_field_elem_mul(&term, &verified->coefficients[i], &request->powers[i], &request->vk.field);
        ext_field_elem_add(&sum, &sum, &term, &request->vk.field);
    }
    int ok = fmpz_is_zero(sum.b);
    if (ok) fq_set_fmpz(out, sum.a, par->Fp);
    ext_field_elem_clear(&sum); ext_field_elem_clear(&term);
    return ok;
}

void scheme1_probgen_out_clear(scheme1_probgen_out *out, pubpar *par) {
    (void)par;
    for (int i = 0; i < out->vk.count; i++) {
        ext_field_elem_clear(&out->powers[i]); ext_field_elem_clear(&out->masks[i]);
        ext_field_mat_clear(&out->shares[i]); free(out->vk.m[i]); free(out->vk.n[i]);
    }
    free(out->powers); free(out->masks); free(out->shares); free(out->vk.m); free(out->vk.n);
    if (out->vk.count) ext_field_ctx_clear(&out->vk.field);
    memset(out, 0, sizeof(*out));
}
void scheme1_verify_out_clear(scheme1_verify_out *out, const pubpar *par) {
    if (out->coefficients) for (int i = 0; i < par->k; i++) ext_field_elem_clear(&out->coefficients[i]);
    free(out->coefficients); memset(out, 0, sizeof(*out));
}
