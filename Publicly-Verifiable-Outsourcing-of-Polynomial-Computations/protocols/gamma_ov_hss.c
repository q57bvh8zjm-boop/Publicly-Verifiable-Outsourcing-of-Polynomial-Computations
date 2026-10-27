#include "gamma_ov_hss.h"

int gamma_ov_hss_keygen(gamma_ov_hss_keygen_out *out, pubpar *par) {
    memset(out, 0, sizeof(*out));
    if (!pv_field_matches_group(par)) return 0;
    return hss_setup(&out->pk, out->ek, par);
}

int gamma_ov_hss_probgen(gamma_ov_hss_probgen_out *out, const hss_public_key *pk,
                    fq_mat_t input, pubpar *par) {
    memset(out, 0, sizeof(*out));
    if (!pv_field_matches_group(par) || par->m > pk->max_m || par->d > pk->max_d ||
        fq_mat_nrows(input, par->Fp) != par->m || fq_mat_ncols(input, par->Fp) != 1)
        return 0;
    fq_t *points = pv_alloc((size_t)par->d + 1, sizeof(fq_t));
    fq_t *weights = pv_alloc((size_t)par->d + 1, sizeof(fq_t));
    for (int j = 0; j <= par->d; j++) {
        fq_init(points[j], par->Fp); fq_init(weights[j], par->Fp);
    }
    pv_sample_points(points, par->d, par->Fp);
    pv_orthogonal_weights(weights, points, par->d, par->Fp);
    out->vk_x = pv_alloc((size_t)par->d + 1, sizeof(*out->vk_x));
    out->ct_inputs = pv_alloc(((size_t)par->d + 1) * par->m, sizeof(hss_ciphertext));
    fmpz_t scalar; fmpz_init(scalar);
    fq_t value; fq_init(value, par->Fp);
    for (int j = 0; j <= par->d; j++) {
        out->vk_x[j] = pv_alloc(32, 1);
        fq2fmpz(scalar, weights[j], par->Fp);
        if (!pv_group_base(out->vk_x[j], scalar)) abort();
        for (int row = 0; row < par->m; row++) {
            fq_mul(value, points[j], fq_mat_entry(input, row, 0), par->Fp);
            hss_input(&out->ct_inputs[j*par->m+row], pk, value, par->Fp);
        }
        fq_clear(points[j], par->Fp); fq_clear(weights[j], par->Fp);
    }
    fq_clear(value, par->Fp); fmpz_clear(scalar); free(points); free(weights);
    return 1;
}

int gamma_ov_hss_compute(gamma_ov_hss_compute_out *out, int server,
                     const hss_evaluation_key *ek, const hss_public_key *pk,
                     fq_mat_t polynomial, hss_ciphertext *inputs, pubpar *par) {
    memset(out, 0, sizeof(*out));
    out->values = pv_alloc((size_t)par->d + 1, sizeof(fq_t));
    out->valid = 1;
    for (int j = 0; j <= par->d; j++) {
        fq_init(out->values[j], par->Fp);
        out->valid &= hss_eval(out->values[j], server, pk, ek,
                               &inputs[j*par->m], polynomial, par);
        /* Server 1 carries the constant; OV checks the polynomial minus F(0). */
        if (server == 1 && out->valid)
            fq_sub(out->values[j], out->values[j], fq_mat_entry(polynomial, 0, 0), par->Fp);
    }
    return out->valid;
}

int gamma_ov_hss_verify(gamma_ov_hss_verify_out *out, unsigned char **vk,
                    gamma_ov_hss_compute_out responses[2], pubpar *par) {
    fq_mat_init(out->aux_R_prime, 1, 1, par->Fp); out->valid = 0;
    if (!responses[0].valid || !responses[1].valid || !pv_field_matches_group(par)) return 0;
    fq_t *values = pv_alloc((size_t)par->d + 1, sizeof(fq_t));
    for (int j = 0; j <= par->d; j++) {
        fq_init(values[j], par->Fp);
        fq_add(values[j], responses[0].values[j], responses[1].values[j], par->Fp);
    }
    out->valid = pv_ov_verify(vk, values, par->d, par);
    if (out->valid) fq_set(fq_mat_entry(out->aux_R_prime, 0, 0), values[0], par->Fp);
    for (int j = 0; j <= par->d; j++) fq_clear(values[j], par->Fp);
    free(values); return out->valid;
}
void gamma_ov_hss_reconstruct(gamma_ov_hss_reconstruct_out *out, fq_mat_t eta,
                         const fq_t constant, pubpar *par) {
    fq_mat_init(out->result, 1, 1, par->Fp);
    fq_add(fq_mat_entry(out->result, 0, 0), fq_mat_entry(eta, 0, 0), constant, par->Fp);
}
void gamma_ov_hss_keygen_out_clear(gamma_ov_hss_keygen_out *out) {
    hss_evaluation_key_clear(&out->ek[0]); hss_evaluation_key_clear(&out->ek[1]);
    hss_public_key_clear(&out->pk);
}
void gamma_ov_hss_probgen_out_clear(gamma_ov_hss_probgen_out *out, pubpar *par) {
    if (out->vk_x) for (int j = 0; j <= par->d; j++) free(out->vk_x[j]);
    if (out->ct_inputs) for (size_t i = 0; i < ((size_t)par->d+1)*par->m; i++)
        hss_ciphertext_clear(&out->ct_inputs[i]);
    free(out->vk_x); free(out->ct_inputs); memset(out, 0, sizeof(*out));
}
void gamma_ov_hss_compute_out_clear(gamma_ov_hss_compute_out *out, pubpar *par) {
    if (out->values) for (int j = 0; j <= par->d; j++) fq_clear(out->values[j], par->Fp);
    free(out->values); memset(out, 0, sizeof(*out));
}
void gamma_ov_hss_verify_out_clear(gamma_ov_hss_verify_out *out, pubpar *par) {
    fq_mat_clear(out->aux_R_prime, par->Fp);
}
void gamma_ov_hss_reconstruct_out_clear(gamma_ov_hss_reconstruct_out *out, pubpar *par) {
    fq_mat_clear(out->result, par->Fp);
}
