#include "gamma_ov_sh.h"

void gamma_ov_sh_keygen(fq_mat_t F, pubpar *par) { (void)F; (void)par; }

int gamma_ov_sh_probgen(gamma_ov_sh_probgen_out *out, fq_mat_t x, pubpar *par) {
    memset(out, 0, sizeof(*out));
    int d = par->d, k = par->k, t = par->t, m = par->m;
    if (!pv_field_matches_group(par) || d < 1 || t < 1 ||
        (long long)d * t + 1 != k || fq_mat_nrows(x, par->Fp) != m ||
        fq_mat_ncols(x, par->Fp) != 1) return 0;
    fq_t *v = pv_alloc((size_t)d + 1, sizeof(fq_t));
    fq_t *weights = pv_alloc((size_t)d + 1, sizeof(fq_t));
    for (int j = 0; j <= d; j++) {
        fq_init(v[j], par->Fp); fq_init(weights[j], par->Fp);
    }
    pv_sample_points(v, d, par->Fp);
    pv_orthogonal_weights(weights, v, d, par->Fp);
    out->vk_x = pv_alloc((size_t)d + 1, sizeof(*out->vk_x));
    out->v_values = pv_alloc((size_t)d, sizeof(fq_t));
    fmpz_t scalar;
    fmpz_init(scalar);
    for (int j = 0; j <= d; j++) {
        out->vk_x[j] = pv_alloc(32, 1);
        fq2fmpz(scalar, weights[j], par->Fp);
        if (!pv_group_base(out->vk_x[j], scalar)) abort();
        if (j) {
            fq_init(out->v_values[j-1], par->Fp);
            fq_set(out->v_values[j-1], v[j], par->Fp);
        }
    }
    fmpz_clear(scalar);
    out->sigma = pv_alloc((size_t)k, sizeof(*out->sigma));
    for (int i = 0; i < k; i++) {
        out->sigma[i] = pv_alloc((size_t)d + 1, sizeof(fq_mat_t));
        for (int j = 0; j <= d; j++) fq_mat_init(out->sigma[i][j], m, 1, par->Fp);
    }
    fq_t *coeff = pv_alloc((size_t)t, sizeof(fq_t));
    for (int s = 0; s < t; s++) fq_init(coeff[s], par->Fp);
    fq_t base, value;
    fq_init(base, par->Fp); fq_init(value, par->Fp);
    for (int j = 0; j <= d; j++) for (int row = 0; row < m; row++) {
        /* Independent, uniform coefficients for every coordinate and column. */
        for (int s = 0; s < t; s++) pv_random_field(coeff[s], par->Fp);
        fq_mul(base, v[j], fq_mat_entry(x, row, 0), par->Fp);
        for (int i = 0; i < k; i++) {
            fq_set(value, coeff[t-1], par->Fp);
            for (int s = t-2; s >= 0; s--) {
                fq_mul(value, value, par->SID[i], par->Fp);
                fq_add(value, value, coeff[s], par->Fp);
            }
            fq_mul(value, value, par->SID[i], par->Fp);
            fq_add(fq_mat_entry(out->sigma[i][j], row, 0), value, base, par->Fp);
        }
    }
    for (int s = 0; s < t; s++) fq_clear(coeff[s], par->Fp);
    free(coeff); fq_clear(base, par->Fp); fq_clear(value, par->Fp);
    for (int j = 0; j <= d; j++) {
        fq_clear(v[j], par->Fp); fq_clear(weights[j], par->Fp);
    }
    free(v); free(weights);
    return 1;
}

void gamma_ov_sh_compute(gamma_ov_sh_compute_out *out, int server_id, fq_mat_t F,
                     fq_mat_t *sigma_i, pubpar *par) {
    (void)server_id;
    out->pi = pv_alloc((size_t)par->d + 1, sizeof(fq_mat_t));
    for (int j = 0; j <= par->d; j++) {
        fq_mat_init(out->pi[j], 1, 1, par->Fp);
        fq_t result;
        Eval(result, F, sigma_i[j], par);
        fq_sub(fq_mat_entry(out->pi[j], 0, 0), result,
               fq_mat_entry(F, 0, 0), par->Fp);
        fq_clear(result, par->Fp);
    }
}

int gamma_ov_sh_verify(gamma_ov_sh_verify_out *out, unsigned char **vk_x,
                   gamma_ov_sh_compute_out *responses, pubpar *par) {
    fq_mat_init(out->aux_R_prime, 1, 1, par->Fp);
    out->valid = 0;
    if (!pv_field_matches_group(par) || !responses || !vk_x) return 0;
    fq_t *values = pv_alloc((size_t)par->d + 1, sizeof(fq_t));
    fq_t *column = pv_alloc((size_t)par->k, sizeof(fq_t));
    for (int i = 0; i < par->k; i++) fq_init(column[i], par->Fp);
    for (int j = 0; j <= par->d; j++) {
        fq_init(values[j], par->Fp);
        for (int i = 0; i < par->k; i++)
            fq_set(column[i], fq_mat_entry(responses[i].pi[j], 0, 0), par->Fp);
        pv_interpolate_zero(values[j], column, par);
    }
    out->valid = pv_ov_verify(vk_x, values, par->d, par);
    if (out->valid) fq_set(fq_mat_entry(out->aux_R_prime, 0, 0), values[0], par->Fp);
    for (int i = 0; i < par->k; i++) fq_clear(column[i], par->Fp);
    for (int j = 0; j <= par->d; j++) fq_clear(values[j], par->Fp);
    free(column); free(values);
    return out->valid;
}

void gamma_ov_sh_reconstruct(gamma_ov_sh_reconstruct_out *out, fq_mat_t eta,
                         const fq_t constant, pubpar *par) {
    fq_mat_init(out->result, 1, 1, par->Fp);
    fq_add(fq_mat_entry(out->result, 0, 0), fq_mat_entry(eta, 0, 0), constant, par->Fp);
}

void gamma_ov_sh_probgen_out_clear(gamma_ov_sh_probgen_out *out, pubpar *par) {
    if (out->sigma) for (int i = 0; i < par->k; i++) {
        for (int j = 0; j <= par->d; j++) fq_mat_clear(out->sigma[i][j], par->Fp);
        free(out->sigma[i]);
    }
    if (out->vk_x) for (int j = 0; j <= par->d; j++) free(out->vk_x[j]);
    if (out->v_values) for (int j = 0; j < par->d; j++) fq_clear(out->v_values[j], par->Fp);
    free(out->sigma); free(out->vk_x); free(out->v_values); memset(out, 0, sizeof(*out));
}
void gamma_ov_sh_compute_out_clear(gamma_ov_sh_compute_out *out, pubpar *par) {
    if (out->pi) for (int j = 0; j <= par->d; j++) fq_mat_clear(out->pi[j], par->Fp);
    free(out->pi); out->pi = NULL;
}
void gamma_ov_sh_verify_out_clear(gamma_ov_sh_verify_out *out, pubpar *par) {
    fq_mat_clear(out->aux_R_prime, par->Fp);
}
void gamma_ov_sh_reconstruct_out_clear(gamma_ov_sh_reconstruct_out *out, pubpar *par) {
    fq_mat_clear(out->result, par->Fp);
}
