#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "protocol2.h"

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
    printf("Benchmarking Protocol 2 client operations...\n\n");

    // Initialize required structures for benchmarking
    pubpar_p2 par;
    fmpz_init(par.p);
    fmpz_set_str(par.p, "2410312426921032588580116606028314112912093247945688951359675039065257391591803200669085024107346049663448766280888004787862416978794958324969612987890774651455213339381625224770782077917681499676845543137387820057597345857904599109461387122099507964997815641342300677629473355281617428411794163967785870370368969109221591943054232011562758450080579587850900993714892283476646631181515063804873375182260506246992837898705971012525843324401232986857004760339321639", 10);
    fq_ctx_init(par.Fq2, par.p, 2, "a");

    par.m = 100;
    par.t = 1;
    par.k = 2;
    par.SID = malloc(sizeof(fq_t) * par.k);
    for(int i=0; i<par.k; i++) {
        fq_init(par.SID[i], par.Fq2);
        fq_set_ui(par.SID[i], i+1, par.Fq2);
    }

    fq_mat_t X;
    fq_mat_init(X, par.m, 1, par.Fq2);
    fq_mat_zero(X, par.Fq2);

    fq_t alpha;
    fq_mat_t *c = malloc(sizeof(fq_mat_t) * par.k);

    // Benchmark core operation in ProbGen: m * k * t field multiplications
    // Actually we can just benchmark protocol2_probgen directly or sum of field ops
    clock_t start = clock();
    int iterations = 100;
    for (int i = 0; i < iterations; i++) {
        protocol2_probgen(&alpha, c, X, &par);
        for(int j=0; j<par.k; j++) fq_mat_clear(c[j], par.Fq2);
        fq_clear(alpha, par.Fq2);
    }
    clock_t end = clock();
    double t_probgen_base = (double)(end - start) / CLOCKS_PER_SEC / iterations;
    // Normalize to m*k*t = 100*2*1
    double t_per_unit = t_probgen_base / (par.m * par.k * par.t);

    // Benchmark fq_mul (Server-side) over Fq2
    flint_rand_t s_state;
    flint_randinit(s_state);
    fq_t a, b, res;
    fq_init(a, par.Fq2); fq_init(b, par.Fq2); fq_init(res, par.Fq2);
    fq_randtest(a, s_state, par.Fq2); fq_randtest(b, s_state, par.Fq2);
    
    start = clock();
    int f_iterations = 1000000;
    for (int i = 0; i < f_iterations; i++) {
        fq_mul(res, a, b, par.Fq2);
    }
    end = clock();
    double t_fq_mul = (double)(end - start) / CLOCKS_PER_SEC / f_iterations;

    flint_randclear(s_state);

    printf("Avg field op unit time: %.9f s\n", t_per_unit);
    printf("Avg server field mul time: %.9f s\n", t_fq_mul);
    printf("\n%-3s | %-3s | %-5s | %-5s | %-12s | %-15s | %-15s\n", "t", "d", "m", "k", "m choose d", "Client Time (s)", "Server Time (s)");
    printf("----------------------------------------------------------------------------------------------------\n");

    for (int t = 1; t <= 5; t++) {
        double min_time = 1e18;
        int best_d = 0, best_m = 0, best_k = 0;

        for (int d = 1; d <= 20; d++) {
            int m = find_m(d, n_db);
            int k = d * t + 1; // Privacy t and degree d implies k > dt for IT-verifiability
            
            // Scaled cost: m * k * t * t_per_unit
            double client_time = (double)m * k * t * t_per_unit;
            
            if (client_time < min_time) {
                min_time = client_time;
                best_d = d;
                best_m = m;
                best_k = k;
            }
        }
        double best_server_s = (double)n_db * best_d * t_fq_mul;
        printf("%-3d | %-3d | %-5d | %-5d | %-12.0f | %-15.6f | %-15.6f (BEST)\n", t, best_d, best_m, best_k, nCr(best_m, best_d), min_time, best_server_s);
    }

    // Cleanup
    fq_mat_clear(X, par.Fq2);
    for(int i=0; i<par.k; i++) fq_clear(par.SID[i], par.Fq2);
    free(par.SID);
    free(c);
    fq_ctx_clear(par.Fq2);
    fmpz_clear(par.p);

    return 0;
}
