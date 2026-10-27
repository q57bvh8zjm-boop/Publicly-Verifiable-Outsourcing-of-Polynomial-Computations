#include "pi3.h"
int pi3_probgen(pi3_probgen_out *out, fq_mat_t x, pubpar *p) {
    memset(out, 0, sizeof(*out));
    if (!pv_field_matches_group(p) || p->d < 1 || p->t < 1 || p->k != p->d * p->t + 1 ||
        fq_mat_nrows(x, p->Fp) != p->m || fq_mat_ncols(x, p->Fp) != 1)
        return 0;
    ext_field_ctx_init_default(&out->field, p->p);
    ext_field_elem_init(&out->alpha);
    out->count = p->k;
    do {
        ext_field_elem_random(&out->alpha, &out->field);
    } while (fmpz_is_zero(out->alpha.b) && fmpz_cmp_ui(out->alpha.a, (ulong)p->k) <= 0);
    ext_field_elem_t *powers = pv_alloc(p->t + 1, sizeof(*powers));
    for (int s = 0; s <= p->t; s++) {
        ext_field_elem_init(&powers[s]);
        if (!s)
            ext_field_elem_one(&powers[s]);
        else
            ext_field_elem_mul(&powers[s], &powers[s - 1], &out->alpha, &out->field);
    }
    out->shares = pv_alloc(p->k, sizeof(*out->shares));
    for (int i = 0; i < p->k; i++)
        ext_field_mat_init(&out->shares[i], p->m, 1);
    ext_field_elem_t *random = pv_alloc(p->t, sizeof(*random));
    for (int s = 0; s < p->t; s++)
        ext_field_elem_init(&random[s]);
    ext_field_elem_t u, power, difference, term;
    ext_field_elem_init(&u);
    ext_field_elem_init(&power);
    ext_field_elem_init(&difference);
    ext_field_elem_init(&term);
    for (int row = 0; row < p->m; row++) {
        for (int s = 0; s < p->t; s++)
            ext_field_elem_random(&random[s], &out->field);
        for (int i = 0; i < p->k; i++) {
            ext_field_elem_t *share = ext_field_mat_entry(&out->shares[i], row, 0);
            fq2fmpz(share->a, fq_mat_entry(x, row, 0), p->Fp);
            fmpz_set_ui(u.a, i + 1);
            fmpz_zero(u.b);
            ext_field_elem_one(&power);
            for (int s = 1; s <= p->t; s++) {
                ext_field_elem_mul(&power, &power, &u, &out->field);
                ext_field_elem_sub(&difference, &power, &powers[s], &out->field);
                ext_field_elem_mul(&term, &random[s - 1], &difference, &out->field);
                ext_field_elem_add(share, share, &term, &out->field);
            }
        }
    }
    for (int s = 0; s < p->t; s++)
        ext_field_elem_clear(&random[s]);
    free(random);
    for (int s = 0; s <= p->t; s++)
        ext_field_elem_clear(&powers[s]);
    free(powers);
    ext_field_elem_clear(&u);
    ext_field_elem_clear(&power);
    ext_field_elem_clear(&difference);
    ext_field_elem_clear(&term);
    return 1;
}
void pi3_compute(ext_field_elem_t *out, fq_mat_t F, ext_field_mat_t *share, ext_field_ctx_t *field,
                 pubpar *p) {
    gamma_fe_sh_eval_extended(out, F, share, field, p);
}
int pi3_verify(pi3_verify_out *out, const pi3_probgen_out *req, ext_field_elem_t *responses,
               const pubpar *p) {
    ext_field_elem_init(&out->value);
    out->valid = 0;
    if (!responses || req->count != p->k || !pv_field_matches_group(p))
        return 0;
    ext_field_elem_t *coeff = pv_alloc(p->k, sizeof(*coeff));
    for (int i = 0; i < p->k; i++)
        ext_field_elem_init(&coeff[i]);
    ext_field_interpolate(coeff, responses, p->k, &req->field);
    /* Horner evaluation at the private point, followed by base-field membership. */
    ext_field_elem_copy(&out->value, &coeff[p->k - 1]);
    for (int i = p->k - 2; i >= 0; i--) {
        ext_field_elem_mul(&out->value, &out->value, &req->alpha, &req->field);
        ext_field_elem_add(&out->value, &out->value, &coeff[i], &req->field);
    }
    out->valid = fmpz_is_zero(out->value.b);
    for (int i = 0; i < p->k; i++)
        ext_field_elem_clear(&coeff[i]);
    free(coeff);
    return out->valid;
}
int pi3_decode(fq_t out, const pi3_verify_out *v, const pubpar *p) {
    if (!v->valid)
        return 0;
    fq_set_fmpz(out, v->value.a, p->Fp);
    return 1;
}
void pi3_probgen_out_clear(pi3_probgen_out *out) {
    for (int i = 0; i < out->count; i++)
        ext_field_mat_clear(&out->shares[i]);
    free(out->shares);
    if (out->count) {
        ext_field_elem_clear(&out->alpha);
        ext_field_ctx_clear(&out->field);
    }
    memset(out, 0, sizeof(*out));
}
void pi3_verify_out_clear(pi3_verify_out *out) {
    ext_field_elem_clear(&out->value);
    out->valid = 0;
}
