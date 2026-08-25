#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "scheme3.h"

// Function to calculate nCr
double nCr(int n, int r) {
    if (r < 0 || r > n) return 0;
    if (r == 0 || r == n) return 1;
    if (r > n / 2) r = n - r;
    
    double res = 1;
    for (int i = 1; i <= r; i++) {
        res = res * (n - r + i) / i;
    }
    return res;
}

// Find minimal m such that mCr >= N
int find_m(int r, long long N) {
    int m = r;
    while (nCr(m, r) < N) {
        m++;
    }
    return m;
}

int main() {
    long long n_db = 1LL << 20; // 2^20
    printf("Database size n = %lld\n", n_db);
    printf("Benchmarking Scheme 3 client operations...\n\n");

    // Initialize required structures for benchmarking
    pubpar par;
    par.m = 100; // dummy
    par.d = 2;   // dummy
    fmpz_init(par.p);
    fmpz_set_str(par.p, "2410312426921032588580116606028314112912093247945688951359675039065257391591803200669085024107346049663448766280888004787862416978794958324969612987890774651455213339381625224770782077917681499676845543137387820057597345857904599109461387122099507964997815641342300677629473355281617428411794163967785870370368969109221591943054232011562758450080579587850900993714892283476646631181515063804873375182260506246992837898705971012525843324401232986857004760339321639", 10);
    fq_ctx_init(par.Fp, par.p, 1, "a");

    paillier_pk_t pk;
    fmpz_t sk;
    fmpz_init(sk);
    hss_ek_t ek0, ek1;
    hss_gen(&pk, sk, &ek0, &ek1, 1024);

    fq_t x;
    fq_init(x, par.Fp);
    fq_set_ui(x, 123, par.Fp);
    paillier_ct_t ct;

    // Benchmark hss_input
    clock_t start = clock();
    int iterations = 100;
    for (int i = 0; i < iterations; i++) {
        hss_input(&ct, &pk, x, &par);
    }
    clock_t end = clock();
    double t_hss_input = (double)(end - start) / CLOCKS_PER_SEC / iterations;

    // Benchmark Ristretto Scalarmult
    unsigned char vk[32];
    unsigned char scalar[32] = {1};
    unsigned char point[32] = {0}; // identity or random
    crypto_scalarmult_ristretto255_base(point, scalar); // make a valid point

    start = clock();
    for (int i = 0; i < iterations; i++) {
        if (crypto_scalarmult_ristretto255(vk, scalar, point) != 0)
            memset(vk, 0, sizeof(vk));
    }
    end = clock();
    double t_ristretto = (double)(end - start) / CLOCKS_PER_SEC / iterations;

    // Benchmark hss_mul_paillier (Server-side)
    hss_mv_t mv_x, mv_y, mv_res;
    paillier_ct_t ct_x;
    hss_input(&ct_x, &pk, x, &par);
    fmpz_init(mv_x.mv[0]); fmpz_init(mv_x.mv[1]);
    fmpz_init(mv_y.mv[0]); fmpz_init(mv_y.mv[1]);
    fmpz_set_ui(mv_x.mv[0], 1); fmpz_set_ui(mv_x.mv[1], 1);
    fmpz_set_ui(mv_y.mv[0], 1); fmpz_set_ui(mv_y.mv[1], 1);

    int prf_key = 0;
    start = clock();
    int s_iterations = 20; // Paillier mul is slow
    for (int i = 0; i < s_iterations; i++) {
        hss_mul_paillier(&mv_res, 0, &pk, &ct_x, &mv_y, &prf_key);
    }
    end = clock();
    double t_hss_mul = (double)(end - start) / CLOCKS_PER_SEC / s_iterations;

    printf("Avg hss_input time: %.6f s\n", t_hss_input);
    printf("Avg hss_mul time (Server): %.6f s\n", t_hss_mul);
    printf("Avg Ristretto ScalarMult time: %.6f s\n", t_ristretto);
    printf("\n%-3s | %-3s | %-5s | %-12s | %-15s | %-15s\n", "t", "d", "m", "m choose d", "Client Time (s)", "Server Time (h)");
    printf("----------------------------------------------------------------------------------------\n");

    for (int t = 1; t <= 5; t++) {
        double min_time = 1e18;
        int best_d = 0;
        int best_m = 0;

        for (int d = 1; d <= 20; d++) {
            int m = find_m(d, n_db);
            // Scaled cost: (d+t)*m * t_hss_input + (d+t) * t_ristretto
            double client_time = (d + t) * m * t_hss_input + (d + t) * t_ristretto;
            
            if (client_time < min_time) {
                min_time = client_time;
                best_d = d;
                best_m = m;
            }
            
            // Print top few or just the best if t > 1? 
            // Better print best for each t to answer the user query clearly.
        }
        double best_server_h = (double)n_db * (best_d - 1) * t_hss_mul / 3600.0;
        printf("%-3d | %-3d | %-5d | %-12.0f | %-15.6f | %-15.2f (BEST)\n", t, best_d, best_m, nCr(best_m, best_d), min_time, best_server_h);
        
        // Also print details for the user's focus
        for (int d = best_d - 2; d <= best_d + 2; d++) {
            if (d < 1 || d > 20) continue;
            int m = find_m(d, n_db);
            double client_time = (d + t) * m * t_hss_input + (d + t) * t_ristretto;
            if (d == best_d) continue;
            printf("%-3d | %-3d | %-5d | %-12.0f | %-15.6f\n", t, d, m, nCr(m, d), client_time);
        }
        printf("------------------------------------------------------------------\n");
    }

    // Cleanup
    fq_clear(x, par.Fp);
    paillier_pk_clear(&pk);
    fmpz_clear(sk);
    hss_ek_clear(&ek0);
    hss_ek_clear(&ek1);
    paillier_ct_clear(&ct);
    paillier_ct_clear(&ct_x);
    hss_mv_clear(&mv_x); hss_mv_clear(&mv_y); hss_mv_clear(&mv_res);
    fq_ctx_clear(par.Fp);
    fmpz_clear(par.p);

    return 0;
}
