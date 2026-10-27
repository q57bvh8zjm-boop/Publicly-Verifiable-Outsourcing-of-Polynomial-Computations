#define _POSIX_C_SOURCE 200809L
#include "wire.h"
#include <signal.h>
#include <unistd.h>

int main(int argc, char **argv) {
    if (argc != 3) {
        fputs("Usage: vm_worker SERVER_ID MAX_COEFFICIENTS (binary TCP/stdin protocol)\n", stderr);
        return 2;
    }
    signal(SIGPIPE, SIG_IGN);
    int server_id = vm_number(argv[1], 0, VM_MAX_SERVERS - 1);
    uint64_t limit = (uint64_t)vm_number(argv[2], 1, INT_MAX);
    wire_buf b = {0}, reply = {0};
    vm_expect_frame(STDIN_FILENO, VM_SETUP, &b);
    vm_check(wire_get_u32(&b) == VM_WIRE_VERSION, "wire version mismatch");
    int scheme = (int)wire_get_u32(&b), m = (int)wire_get_u32(&b), d = (int)wire_get_u32(&b),
        t = (int)wire_get_u32(&b);
    int sid = (int)wire_get_u32(&b);
    unsigned char job[32];
    wire_get(&b, job, 32);
    vm_check(scheme >= 1 && scheme <= 4 && m > 0 && d > 0 && d <= 32 && t > 0 && t <= 32 &&
                 d * t + 1 <= VM_MAX_SERVERS,
             "invalid setup dimensions");
    vm_check(sid == server_id && sid < (scheme == 3 ? 2 : d * t + 1) && (scheme != 3 || t == 1),
             "wrong server or HSS threshold");
    slong count = vm_coefficient_count(m, d, limit);
    pubpar p;
    vm_check(pv_params_init(&p, m, d, t), "parameters");
    hss_public_key hpk;
    hss_evaluation_key hek;
    if (scheme == 3)
        wire_get_hss_setup(&b, &hpk, &hek, &p);
    wire_end(&b);
    wire_clear(&b);
    fq_mat_t F;
    fq_mat_init(F, count, 1, p.Fp);
    vm_fill_polynomial(F, &p);
    ext_field_ctx_t field;
    if (scheme == 1 || scheme == 4)
        ext_field_ctx_init_default(&field, p.p);
    wire_u32(&reply, VM_WIRE_VERSION);
    wire_u32(&reply, (uint32_t)sid);
    wire_put(&reply, job, 32);
    char info[512];
    int n = snprintf(info, sizeof(info), "FLINT=%s GMP=%s sodium=%s compiler=%s", FLINT_VERSION,
                     gmp_version, sodium_version_string(), __VERSION__);
    wire_u32(&reply, (uint32_t)n);
    wire_put(&reply, info, (size_t)n);
    wire_frame_write(STDOUT_FILENO, VM_READY, &reply);
    wire_clear(&reply);
    uint32_t expected_trial = 0;
    for (;;) {
        uint32_t type;
        if (!wire_frame_read(STDIN_FILENO, &type, &b))
            break;
        if (type == VM_BYE) {
            wire_end(&b);
            wire_clear(&b);
            break;
        }
        vm_check(type == VM_QUERY, "expected query");
        unsigned char query_job[32];
        wire_get(&b, query_job, 32);
        uint32_t trial = wire_get_u32(&b);
        vm_check(!memcmp(job, query_job, 32) && trial == expected_trial++,
                 "stale or reordered query");
        uint32_t key_count = wire_get_u32(&b);
        uint32_t required = scheme == 1 ? 2 * (uint32_t)p.k : scheme == 4 ? 0 : (uint32_t)p.d + 1;
        vm_check(key_count == required, "verification-key size mismatch");
        for (uint32_t i = 0; i < key_count; i++) {
            unsigned char point[32];
            wire_get(&b, point, 32);
            vm_check(crypto_core_ristretto255_is_valid_point(point),
                     "invalid public verification key");
        }
        ext_field_mat_t share;
        ext_field_elem_t mask;
        fq_mat_t *columns = NULL;
        hss_ciphertext *ct = NULL;
        if (scheme == 1 || scheme == 4) {
            ext_field_mat_init(&share, m, 1);
            for (int i = 0; i < m; i++)
                wire_get_ext(&b, ext_field_mat_entry(&share, i, 0), &p);
            if (scheme == 1) {
                ext_field_elem_init(&mask);
                wire_get_ext(&b, &mask, &p);
            }
        } else if (scheme == 2) {
            columns = pv_alloc(d + 1, sizeof(fq_mat_t));
            for (int j = 0; j <= d; j++) {
                fq_mat_init(columns[j], m, 1, p.Fp);
                for (int i = 0; i < m; i++)
                    wire_get_fq(&b, fq_mat_entry(columns[j], i, 0), &p);
            }
        } else
            ct = wire_get_hss_inputs(&b, &hpk, &p);
        wire_end(&b);
        wire_clear(&b);
        wire_u32(&reply, trial);
        wire_frame_write(STDOUT_FILENO, VM_LOADED, &reply);
        wire_clear(&reply);
        vm_expect_frame(STDIN_FILENO, VM_GO, &b);
        vm_check(wire_get_u32(&b) == trial, "start token mismatch");
        wire_end(&b);
        wire_clear(&b);
        ext_field_elem_t ext_response;
        gamma_ov_sh_compute_out ov_response;
        gamma_ov_hss_compute_out hss_response;
        if (scheme == 1 || scheme == 4)
            ext_field_elem_init(&ext_response);
        uint64_t begin = vm_ns();
        if (scheme == 1)
            gamma_fe_sh_compute(&ext_response, F, &share, &mask, &field, &p);
        else if (scheme == 2)
            gamma_ov_sh_compute(&ov_response, sid, F, columns, &p);
        else if (scheme == 3)
            vm_check(gamma_ov_hss_compute(&hss_response, sid, &hek, &hpk, F, ct, &p),
                     "HSS Compute failed");
        else
            pi3_compute(&ext_response, F, &share, &field, &p);
        uint64_t elapsed = vm_ns() - begin;
        wire_put(&reply, job, 32);
        wire_u32(&reply, trial);
        wire_u32(&reply, (uint32_t)sid);
        wire_u64(&reply, elapsed);
        if (scheme == 1 || scheme == 4)
            wire_ext(&reply, &ext_response, &p);
        else
            for (int j = 0; j <= d; j++)
                wire_fq(&reply,
                        scheme == 2 ? fq_mat_entry(ov_response.pi[j], 0, 0)
                                    : hss_response.values[j],
                        &p);
        wire_frame_write(STDOUT_FILENO, VM_RESULT, &reply);
        wire_clear(&reply);
        if (scheme == 1 || scheme == 4) {
            ext_field_elem_clear(&ext_response);
            ext_field_mat_clear(&share);
            if (scheme == 1)
                ext_field_elem_clear(&mask);
        } else if (scheme == 2) {
            gamma_ov_sh_compute_out_clear(&ov_response, &p);
            for (int j = 0; j <= d; j++)
                fq_mat_clear(columns[j], p.Fp);
            free(columns);
        } else {
            gamma_ov_hss_compute_out_clear(&hss_response, &p);
            for (size_t i = 0; i < (size_t)m * (d + 1); i++)
                hss_ciphertext_clear(&ct[i]);
            free(ct);
        }
    }
    if (scheme == 3) {
        hss_evaluation_key_clear(&hek);
        hss_public_key_clear(&hpk);
    }
    if (scheme == 1 || scheme == 4)
        ext_field_ctx_clear(&field);
    fq_mat_clear(F, p.Fp);
    pv_params_clear(&p);
    return 0;
}
