#define _POSIX_C_SOURCE 200809L
#include "wire.h"
#include <signal.h>
#include <unistd.h>

static void batch_add(wire_buf *batch, int index, wire_buf *payload) {
    wire_u32(batch, (uint32_t)index);
    wire_u64(batch, payload->len);
    wire_put(batch, payload->data, payload->len);
}
static wire_buf batch_next(wire_buf *batch, int expected) {
    vm_check(wire_get_u32(batch) == (uint32_t)expected, "server response order");
    uint64_t n = wire_get_u64(batch);
    vm_check(n <= batch->len - batch->pos, "truncated server response");
    wire_buf b = {batch->data + batch->pos, (size_t)n, (size_t)n, 0};
    batch->pos += (size_t)n;
    return b;
}
static void send_json(const char *json) {
    wire_buf out = {0};
    wire_put(&out, json, strlen(json));
    wire_frame_write(STDOUT_FILENO, VM_ROW, &out);
    wire_clear(&out);
}
int main(int argc, char **argv) {
    if (argc != 9) {
        fputs("Usage: vm_client SCHEME M D T TRIALS WARMUPS NATIVE_REPEATS MAX_COEFFICIENTS\n",
              stderr);
        return 2;
    }
    signal(SIGPIPE, SIG_IGN);
    fprintf(stderr, "client: FLINT=%s GMP=%s sodium=%s compiler=%s\n", FLINT_VERSION,
            gmp_version, sodium_version_string(), __VERSION__);
    int scheme = vm_scheme(argv[1]), m = vm_number(argv[2], 1, 1000000),
        d = vm_number(argv[3], 1, 32), t = vm_number(argv[4], 1, 32);
    int trials = vm_number(argv[5], 1, 100000), warmups = vm_number(argv[6], 0, 100),
        native_repeats = vm_number(argv[7], 0, 100000);
    uint64_t limit = (uint64_t)vm_number(argv[8], 1, INT_MAX);
    int servers = scheme == 3 ? 2 : d * t + 1;
    vm_check(servers <= VM_MAX_SERVERS && (scheme != 3 || t == 1),
             "unsupported server count or threshold");
    slong count = vm_coefficient_count(m, d, limit);
    pubpar p;
    vm_check(pv_params_init(&p, m, d, t), "parameters");
    fq_mat_t F, x;
    fq_mat_init(F, count, 1, p.Fp);
    fq_mat_init(x, m, 1, p.Fp);
    vm_fill_polynomial(F, &p);
    unsigned char job[32];
    randombytes_buf(job, sizeof(job));
    gamma_ov_hss_keygen_out hkey;
    uint64_t keygen_ns = 0, begin = vm_ns();
    if (scheme == 3) {
        vm_check(gamma_ov_hss_keygen(&hkey, &p), "HSS setup");
        keygen_ns = vm_ns() - begin;
    }
    wire_buf batch = {0}, incoming = {0};
    wire_u32(&batch, (uint32_t)servers);
    for (int i = 0; i < servers; i++) {
        wire_buf setup = {0};
        wire_u32(&setup, VM_WIRE_VERSION);
        wire_u32(&setup, (uint32_t)scheme);
        wire_u32(&setup, (uint32_t)m);
        wire_u32(&setup, (uint32_t)d);
        wire_u32(&setup, (uint32_t)t);
        wire_u32(&setup, (uint32_t)i);
        wire_put(&setup, job, sizeof(job));
        if (scheme == 3)
            wire_hss_setup(&setup, &hkey, i);
        batch_add(&batch, i, &setup);
        wire_clear(&setup);
    }
    wire_frame_write(STDOUT_FILENO, VM_BATCH_SETUP, &batch);
    wire_clear(&batch);
    vm_expect_frame(STDIN_FILENO, VM_BATCH_READY, &incoming);
    uint64_t setup_ns = wire_get_u64(&incoming);
    vm_check(wire_get_u32(&incoming) == (uint32_t)servers, "setup reply count");
    for (int i = 0; i < servers; i++) {
        wire_buf response = batch_next(&incoming, i);
        vm_check(wire_get_u32(&response) == VM_WIRE_VERSION &&
                     wire_get_u32(&response) == (uint32_t)i,
                 "worker version or identity");
        unsigned char received[32];
        wire_get(&response, received, 32);
        vm_check(!memcmp(received, job, 32), "setup session mismatch");
        uint32_t len = wire_get_u32(&response);
        vm_check(len < 1024 && len <= response.len - response.pos, "worker metadata length");
        fprintf(stderr, "server %d: %.*s\n", i, (int)len, (char *)response.data + response.pos);
        response.pos += len;
        wire_end(&response);
    }
    wire_end(&incoming);
    wire_clear(&incoming);
    for (int sequence = 0; sequence < trials + warmups; sequence++) {
        int trial = sequence - warmups;
        for (int i = 0; i < m; i++)
            vm_fixture(fq_mat_entry(x, i, 0), "input", m, d, i, trial, &p);
        begin = vm_ns();
        fq_t expected;
        Eval(expected, F, x, &p);
        uint64_t single = vm_ns() - begin;
        int repeats = native_repeats;
        if (!repeats) {
            repeats = single ? (int)(100000000ULL / single) : 10000;
            if (repeats < 1)
                repeats = 1;
            if (repeats > 10000)
                repeats = 10000;
        }
        begin = vm_ns();
        for (int i = 0; i < repeats; i++) {
            fq_t y;
            Eval(y, F, x, &p);
            fq_clear(y, p.Fp);
        }
        double native_ns = (double)(vm_ns() - begin) / repeats;
        gamma_fe_sh_probgen_out fe;
        gamma_ov_sh_probgen_out ov;
        gamma_ov_hss_probgen_out hs;
        pi3_probgen_out priv;
        uint64_t online_begin = vm_ns();
        begin = vm_ns();
        if (scheme == 1)
            vm_check(gamma_fe_sh_probgen(&fe, x, &p), "FE ProbGen");
        else if (scheme == 2)
            vm_check(gamma_ov_sh_probgen(&ov, x, &p), "OV ProbGen");
        else if (scheme == 3)
            vm_check(gamma_ov_hss_probgen(&hs, &hkey.pk, x, &p), "HSS ProbGen");
        else
            vm_check(pi3_probgen(&priv, x, &p), "Pi3 ProbGen");
        uint64_t pg_ns = vm_ns() - begin;
        begin = vm_ns();
        wire_u32(&batch, (uint32_t)servers);
        for (int i = 0; i < servers; i++) {
            wire_buf query = {0};
            wire_put(&query, job, 32);
            wire_u32(&query, (uint32_t)sequence);
            wire_u32(&query, scheme == 1 ? 2 * (uint32_t)p.k : scheme == 4 ? 0 : (uint32_t)d + 1);
            if (scheme == 1)
                for (int j = 0; j < p.k; j++) {
                    wire_put(&query, fe.vk.m[j], 32);
                    wire_put(&query, fe.vk.n[j], 32);
                }
            else if (scheme == 2 || scheme == 3)
                for (int j = 0; j <= d; j++)
                    wire_put(&query, scheme == 2 ? ov.vk_x[j] : hs.vk_x[j], 32);
            if (scheme == 1 || scheme == 4) {
                ext_field_mat_t *share = scheme == 1 ? &fe.shares[i] : &priv.shares[i];
                for (int j = 0; j < m; j++)
                    wire_ext(&query, ext_field_mat_entry(share, j, 0), &p);
                if (scheme == 1)
                    wire_ext(&query, &fe.masks[i], &p);
            } else if (scheme == 2)
                for (int j = 0; j <= d; j++)
                    for (int r = 0; r < m; r++)
                        wire_fq(&query, fq_mat_entry(ov.sigma[i][j], r, 0), &p);
            else
                wire_hss_inputs(&query, hs.ct_inputs, &hkey.pk, &p);
            batch_add(&batch, i, &query);
            wire_clear(&query);
        }
        uint64_t encode_ns = vm_ns() - begin;
        wire_frame_write(STDOUT_FILENO, VM_BATCH_QUERY, &batch);
        wire_clear(&batch);
        vm_expect_frame(STDIN_FILENO, VM_BATCH_RESULT, &incoming);
        uint64_t upload_ns = wire_get_u64(&incoming), compute_download_ns = wire_get_u64(&incoming),
                 roundtrip_ns = wire_get_u64(&incoming);
        uint64_t request_bytes = wire_get_u64(&incoming), response_bytes = wire_get_u64(&incoming);
        vm_check(wire_get_u32(&incoming) == (uint32_t)servers, "result count");
        ext_field_elem_t *er = NULL;
        gamma_ov_sh_compute_out *orr = NULL;
        gamma_ov_hss_compute_out hr[2];
        if (scheme == 1 || scheme == 4) {
            er = pv_alloc(servers, sizeof(*er));
            for (int i = 0; i < servers; i++)
                ext_field_elem_init(&er[i]);
        } else if (scheme == 2) {
            orr = pv_alloc(servers, sizeof(*orr));
            for (int i = 0; i < servers; i++) {
                orr[i].pi = pv_alloc(d + 1, sizeof(fq_mat_t));
                for (int j = 0; j <= d; j++)
                    fq_mat_init(orr[i].pi[j], 1, 1, p.Fp);
            }
        } else
            for (int i = 0; i < 2; i++) {
                hr[i].valid = 1;
                hr[i].values = pv_alloc(d + 1, sizeof(fq_t));
                for (int j = 0; j <= d; j++)
                    fq_init(hr[i].values[j], p.Fp);
            }
        uint64_t server_ns[VM_MAX_SERVERS], server_max = 0, server_sum = 0, server_min = UINT64_MAX;
        begin = vm_ns();
        for (int i = 0; i < servers; i++) {
            wire_buf response = batch_next(&incoming, i);
            unsigned char got[32];
            wire_get(&response, got, 32);
            vm_check(!memcmp(got, job, 32) && wire_get_u32(&response) == (uint32_t)sequence &&
                         wire_get_u32(&response) == (uint32_t)i,
                     "response session/trial/server mismatch");
            server_ns[i] = wire_get_u64(&response);
            vm_check(server_ns[i] <= UINT64_MAX - server_sum, "server time overflow");
            server_sum += server_ns[i];
            if (server_ns[i] > server_max)
                server_max = server_ns[i];
            if (server_ns[i] < server_min)
                server_min = server_ns[i];
            if (scheme == 1 || scheme == 4)
                wire_get_ext(&response, &er[i], &p);
            else
                for (int j = 0; j <= d; j++)
                    wire_get_fq(&response,
                                scheme == 2 ? fq_mat_entry(orr[i].pi[j], 0, 0) : hr[i].values[j],
                                &p);
            wire_end(&response);
        }
        wire_end(&incoming);
        wire_clear(&incoming);
        uint64_t parse_ns = vm_ns() - begin;
        gamma_fe_sh_verify_out fv;
        gamma_ov_sh_verify_out ovv;
        gamma_ov_hss_verify_out hv;
        pi3_verify_out pv;
        begin = vm_ns();
        int ok;
        if (scheme == 1)
            ok = gamma_fe_sh_verify(&fv, &fe.vk, er, &p);
        else if (scheme == 2)
            ok = gamma_ov_sh_verify(&ovv, ov.vk_x, orr, &p);
        else if (scheme == 3)
            ok = gamma_ov_hss_verify(&hv, hs.vk_x, hr, &p);
        else
            ok = pi3_verify(&pv, &priv, er, &p);
        uint64_t verify_ns = vm_ns() - begin;
        vm_check(ok, "verification failed; no timing row emitted");
        fq_t actual;
        fq_init(actual, p.Fp);
        gamma_ov_sh_reconstruct_out od;
        gamma_ov_hss_reconstruct_out hd;
        begin = vm_ns();
        if (scheme == 1)
            ok = gamma_fe_sh_decode(actual, &fe, &fv, &p);
        else if (scheme == 2)
            gamma_ov_sh_reconstruct(&od, ovv.aux_R_prime, fq_mat_entry(F, 0, 0), &p);
        else if (scheme == 3)
            gamma_ov_hss_reconstruct(&hd, hv.aux_R_prime, fq_mat_entry(F, 0, 0), &p);
        else
            ok = pi3_decode(actual, &pv, &p);
        uint64_t decode_ns = vm_ns() - begin, online_ns = vm_ns() - online_begin;
        if (scheme == 2)
            fq_set(actual, fq_mat_entry(od.result, 0, 0), p.Fp);
        if (scheme == 3)
            fq_set(actual, fq_mat_entry(hd.result, 0, 0), p.Fp);
        vm_check(ok && fq_equal(actual, expected, p.Fp),
                 "decoded result differs from native reference; no timing row emitted");
        if (trial >= 0) {
            char row[16384], times[8192];
            size_t pos = 0;
            times[pos++] = '[';
            for (int i = 0; i < servers; i++)
                pos += (size_t)snprintf(times + pos, sizeof(times) - pos, "%s%.9f", i ? "," : "",
                                        server_ns[i] / 1e6);
            times[pos++] = ']';
            times[pos] = 0;
            double client_ns = (double)pg_ns + verify_ns + decode_ns;
            int n = snprintf(
                row, sizeof(row),
                "{\"mode\":\"tcp_full_protocol\",\"scheme\":\"%s\",\"m\":%d,\"d\":%d,\"t\":%d,"
                "\"servers\":%d,\"coefficients\":%ld,\"trial\":%d,\"native_inner_runs\":%d,"
                "\"native_ms\":%.9f,\"keygen_ms\":%.9f,\"remote_setup_ms\":%.9f,\"probgen_ms\":%."
                "9f,\"server_max_ms\":%.9f,\"server_sum_ms\":%.9f,\"server_min_ms\":%.9f,\"server_"
                "ms\":%s,\"verify_ms\":%.9f,\"decode_ms\":%.9f,\"client_ms\":%.9f,\"client_ratio\":"
                "%.9f,\"request_encode_ms\":%.9f,\"response_parse_ms\":%.9f,\"upload_wall_ms\":%."
                "9f,\"compute_download_wall_ms\":%.9f,\"tcp_roundtrip_ms\":%.9f,\"online_wall_ms\":"
                "%.9f,\"request_wire_bytes\":%llu,\"response_wire_bytes\":%llu,\"paillier_bits\":%"
                "u,\"correct\":1}",
                vm_scheme_name(scheme), m, d, t, servers, (long)count, trial, repeats,
                native_ns / 1e6, keygen_ns / 1e6, setup_ns / 1e6, pg_ns / 1e6, server_max / 1e6,
                server_sum / 1e6, server_min / 1e6, times, verify_ns / 1e6, decode_ns / 1e6,
                client_ns / 1e6, client_ns / native_ns, encode_ns / 1e6, parse_ns / 1e6,
                upload_ns / 1e6, compute_download_ns / 1e6, roundtrip_ns / 1e6, online_ns / 1e6,
                (unsigned long long)request_bytes, (unsigned long long)response_bytes,
                scheme == 3 ? hkey.pk.modulus_bits : 0);
            vm_check(n > 0 && (size_t)n < sizeof(row), "row overflow");
            send_json(row);
        }
        if (scheme == 1 || scheme == 4) {
            for (int i = 0; i < servers; i++)
                ext_field_elem_clear(&er[i]);
            free(er);
            if (scheme == 1) {
                gamma_fe_sh_verify_out_clear(&fv, &p);
                gamma_fe_sh_probgen_out_clear(&fe, &p);
            } else {
                pi3_verify_out_clear(&pv);
                pi3_probgen_out_clear(&priv);
            }
        } else if (scheme == 2) {
            for (int i = 0; i < servers; i++)
                gamma_ov_sh_compute_out_clear(&orr[i], &p);
            free(orr);
            gamma_ov_sh_verify_out_clear(&ovv, &p);
            gamma_ov_sh_reconstruct_out_clear(&od, &p);
            gamma_ov_sh_probgen_out_clear(&ov, &p);
        } else {
            for (int i = 0; i < 2; i++)
                gamma_ov_hss_compute_out_clear(&hr[i], &p);
            gamma_ov_hss_verify_out_clear(&hv, &p);
            gamma_ov_hss_reconstruct_out_clear(&hd, &p);
            gamma_ov_hss_probgen_out_clear(&hs, &p);
        }
        fq_clear(actual, p.Fp);
        fq_clear(expected, p.Fp);
    }
    wire_frame_write(STDOUT_FILENO, VM_DONE, NULL);
    if (scheme == 3)
        gamma_ov_hss_keygen_out_clear(&hkey);
    fq_mat_clear(F, p.Fp);
    fq_mat_clear(x, p.Fp);
    pv_params_clear(&p);
    return 0;
}
