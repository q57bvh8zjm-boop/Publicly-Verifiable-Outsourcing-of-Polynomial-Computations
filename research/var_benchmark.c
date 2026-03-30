#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// Mock up schemes headers to use atomic operations
#include "scheme1.h"
#include "scheme2.h"
#include "scheme3.h"
#include "scheme4.h"
#include "sodium.h"

int main() {
    long long n_users = 1LL << 20; // 2^20 users
    printf("VAR Benchmark with %lld users\n\n", n_users);

    // Context setup
    pubpar par;
    par.m = 100;
    par.d = 1; // Sum is degree 1
    par.k = 3; // Number of servers
    par.t = 1; // Threshold
    fmpz_init(par.p);
    // Typical prime 
    fmpz_set_str(par.p, "2410312426921032588580116606028314112912093247945688951359675039065257391591803200669085024107346049663448766280888004787862416978794958324969612987890774651455213339381625224770782077917681499676845543137387820057597345857904599109461387122099507964997815641342300677629473355281617428411794163967785870370368969109221591943054232011562758450080579587850900993714892283476646631181515063804873375182260506246992837898705971012525843324401232986857004760339321639", 10);
    fq_ctx_init(par.Fp, par.p, 1, "a");
    scheme1_init_ext_ctx(&par);
    ext_field_ctx_t* s1_ctx = scheme1_get_ext_ctx();

    clock_t start, end;

    // ==========================================
    // Scheme 1 & 2 (Assuming F_q2 field ops)
    // ==========================================
    printf("Running Scheme 1/2 Issuance for %lld users...\n", n_users);
    ext_field_elem_t a, b, c;
    ext_field_elem_init(&a); ext_field_elem_init(&b); ext_field_elem_init(&c);
    ext_field_elem_one(&a); ext_field_elem_one(&b);
    
    start = clock();
    for (int i = 0; i < n_users * par.k; i++) {
        ext_field_elem_mul(&c, &a, &b, s1_ctx);
        ext_field_elem_add(&c, &c, &a, s1_ctx);
    }
    end = clock();
    double t_s1_issue = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Running Scheme 1/2 Server Audit for %lld users...\n", n_users);
    start = clock();
    for (int i = 0; i < n_users; i++) {
        ext_field_elem_add(&c, &c, &a, s1_ctx);
    }
    end = clock();
    double t_s1_server = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Running Scheme 1/2 Verification...\n");
    ext_field_elem_t *phi_coeffs = malloc(sizeof(ext_field_elem_t) * par.k);
    ext_field_elem_t *y_vals = malloc(sizeof(ext_field_elem_t) * par.k);
    for(int i=0; i<par.k; i++) { ext_field_elem_init(&y_vals[i]); ext_field_elem_one(&y_vals[i]); }
    
    start = clock();
    scheme1_ext_field_interpolate(phi_coeffs, y_vals, &par);
    end = clock();
    double t_s1_verify = (double)(end - start) / CLOCKS_PER_SEC;

    // ==========================================
    // Scheme 3 & 4 (Paillier)
    // ==========================================
    paillier_pk_t pk;
    fmpz_t sk; fmpz_init(sk);
    hss_ek_t ek0, ek1;
    hss_gen(&pk, sk, &ek0, &ek1, 1024);

    fq_t x; fq_init(x, par.Fp); fq_set_ui(x, 1, par.Fp);
    paillier_ct_t ct1, ct2;

    int n_sample_hss = 1000;
    printf("Running Scheme 3/4 Issuance (sampled %d encryptions to avoid 16h runtime)...\n", n_sample_hss);
    start = clock();
    for (int i = 0; i < n_sample_hss; i++) {
        hss_input(&ct1, &pk, x, &par);
    }
    end = clock();
    double t_s3_issue = ((double)(end - start) / CLOCKS_PER_SEC) / n_sample_hss * n_users; // scale up to 1 million for accurate comparison, otherwise wait 16 hours!

    fmpz_t n_sq; fmpz_init(n_sq);
    fmpz_mul(n_sq, pk.n, pk.n);

    printf("Running Scheme 3/4 Server Audit for %lld users...\n", n_users);
    start = clock();
    for (int i = 0; i < n_users; i++) {
        // Paillier addition
        fmpz_mul(ct1.c[0][0], ct1.c[0][0], ct1.c[0][0]);
        fmpz_mod(ct1.c[0][0], ct1.c[0][0], n_sq);
    }
    end = clock();
    double t_s3_server = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Running Scheme 3/4 Verification...\n");
    hss_mv_t mv0, mv1;
    fmpz_init(mv0.mv[0]); fmpz_init(mv0.mv[1]);
    fmpz_init(mv1.mv[0]); fmpz_init(mv1.mv[1]);
    fq_t out_x; fq_init(out_x, par.Fp);
    start = clock();
    for (int i = 0; i < 1; i++) {
        hss_decrypt(out_x, &pk, &mv0, &mv1, &par);
    }
    end = clock();
    double t_s3_verify = (double)(end - start) / CLOCKS_PER_SEC;

    printf("\nMarkdown Table:\n\n");
    printf("| Scheme | Issuance Time | Audit Proving Time | Audit Verification Time |\n");
    printf("|--------|---------------|--------------------|-------------------------|\n");
    printf("| P-VAR (Paper) | 1.52 s | 9.7 s | < 19 ms |\n");
    printf("| S-VAR (Paper, k=1) | 0.89 s | 34 s | < 100 µs |\n");
    printf("| Our Scheme 1/2 | %.4f s | %.4f s | %.4f ms |\n", t_s1_issue, t_s1_server, t_s1_verify * 1000.0);
    printf("| Our Scheme 3/4 | %.4f s* | %.4f s | %.4f ms |\n", t_s3_issue, t_s3_server, t_s3_verify * 1000.0);
    printf("(* Scheme 3 Issuance time reflects an extrapolated value due to an estimated 16h loop timeline for 1M Paillier encryptions.)\n");

    return 0;
}
