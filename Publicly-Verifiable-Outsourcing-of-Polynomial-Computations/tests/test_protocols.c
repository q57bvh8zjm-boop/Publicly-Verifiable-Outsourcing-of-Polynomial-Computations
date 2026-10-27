#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "gamma_fe_sh.h"
#include "gamma_ov_sh.h"
#include "gamma_ov_hss.h"
#include "pi3.h"

static void check_pi3(fq_mat_t F, fq_mat_t x, pubpar *p) {
    pi3_probgen_out request;
    assert(pi3_probgen(&request, x, p));
    ext_field_elem_t *responses = pv_alloc(p->k, sizeof(*responses));
    for (int i = 0; i < p->k; i++) {
        ext_field_elem_init(&responses[i]);
        pi3_compute(&responses[i], F, &request.shares[i], &request.field, p);
    }
    pi3_verify_out verified;
    assert(pi3_verify(&verified, &request, responses, p));
    fq_t actual, expected;
    fq_init(actual, p->Fp);
    Eval(expected, F, x, p);
    assert(pi3_decode(actual, &verified, p));
    assert(fq_equal(actual, expected, p->Fp));
    pi3_verify_out_clear(&verified);
    for (int i = 0; i < p->k; i++) {
        fmpz_add_ui(responses[i].b, responses[i].b, 1);
        fmpz_mod(responses[i].b, responses[i].b, p->p);
    }
    assert(!pi3_verify(&verified, &request, responses, p));
    assert(!pi3_decode(actual, &verified, p));
    pi3_verify_out_clear(&verified);
    for (int i = 0; i < p->k; i++) ext_field_elem_clear(&responses[i]);
    free(responses);
    fq_clear(actual, p->Fp); fq_clear(expected, p->Fp);
    pi3_probgen_out_clear(&request);
}

static void example(fq_mat_t F, fq_mat_t x, pubpar *p, int random) {
    for (slong i=0; i<fq_mat_nrows(F,p->Fp); i++) {
        if (random) pv_random_field(fq_mat_entry(F,i,0),p->Fp);
        else fq_set_ui(fq_mat_entry(F,i,0),(ulong)i+1,p->Fp);
    }
    for (int i=0;i<p->m;i++) {
        if (random) pv_random_field(fq_mat_entry(x,i,0),p->Fp);
        else fq_set_ui(fq_mat_entry(x,i,0),(ulong)i+2,p->Fp);
    }
}
static void check_fe(fq_mat_t F, fq_mat_t x, pubpar *p) {
    gamma_fe_sh_probgen_out req;
    assert(gamma_fe_sh_probgen(&req,x,p));
    ext_field_elem_t *r=pv_alloc(p->k,sizeof(*r));
    for (int i=0;i<p->k;i++) {
        ext_field_elem_init(&r[i]);
        gamma_fe_sh_compute(&r[i],F,&req.shares[i],&req.masks[i],&req.vk.field,p);
    }
    gamma_fe_sh_verify_out verified;
    assert(gamma_fe_sh_verify(&verified,&req.vk,r,p));
    fq_t actual, expected;
    fq_init(actual,p->Fp); Eval(expected,F,x,p);
    assert(gamma_fe_sh_decode(actual,&req,&verified,p));
    assert(fq_equal(actual,expected,p->Fp));
    gamma_fe_sh_verify_out_clear(&verified,p);
    /* Imaginary constant changes fail public verification deterministically. */
    for (int i=0;i<p->k;i++) {
        fmpz_add_ui(r[i].b,r[i].b,1); fmpz_mod(r[i].b,r[i].b,p->p);
    }
    assert(!gamma_fe_sh_verify(&verified,&req.vk,r,p));
    assert(!gamma_fe_sh_decode(actual,&req,&verified,p));
    gamma_fe_sh_verify_out_clear(&verified,p);
    for (int i=0;i<p->k;i++) {
        fmpz_sub_ui(r[i].b,r[i].b,1); fmpz_mod(r[i].b,r[i].b,p->p);
    }
    memset(req.vk.n[0],0xff,32);
    assert(!gamma_fe_sh_verify(&verified,&req.vk,r,p));
    gamma_fe_sh_verify_out_clear(&verified,p);
    for (int i=0;i<p->k;i++) ext_field_elem_clear(&r[i]);
    free(r); fq_clear(actual,p->Fp); fq_clear(expected,p->Fp);
    gamma_fe_sh_probgen_out_clear(&req,p);
}
static void check_ov(fq_mat_t F, fq_mat_t x, pubpar *p) {
    gamma_ov_sh_probgen_out req, next;
    assert(gamma_ov_sh_probgen(&req,x,p)); assert(gamma_ov_sh_probgen(&next,x,p));
    assert(!fq_mat_equal(req.sigma[0][0],next.sigma[0][0],p->Fp));
    gamma_ov_sh_probgen_out_clear(&next,p);
    gamma_ov_sh_compute_out *r=pv_alloc(p->k,sizeof(*r));
    for (int i=0;i<p->k;i++) gamma_ov_sh_compute(&r[i],i,F,req.sigma[i],p);
    gamma_ov_sh_verify_out v; gamma_ov_sh_reconstruct_out decoded;
    assert(gamma_ov_sh_verify(&v,req.vk_x,r,p));
    gamma_ov_sh_reconstruct(&decoded,v.aux_R_prime,fq_mat_entry(F,0,0),p);
    fq_t expected, one; Eval(expected,F,x,p); fq_init(one,p->Fp); fq_one(one,p->Fp);
    assert(fq_equal(fq_mat_entry(decoded.result,0,0),expected,p->Fp));
    gamma_ov_sh_verify_out_clear(&v,p); gamma_ov_sh_reconstruct_out_clear(&decoded,p);
    fq_add(fq_mat_entry(r[0].pi[0],0,0),fq_mat_entry(r[0].pi[0],0,0),one,p->Fp);
    assert(!gamma_ov_sh_verify(&v,req.vk_x,r,p)); gamma_ov_sh_verify_out_clear(&v,p);
    fq_sub(fq_mat_entry(r[0].pi[0],0,0),fq_mat_entry(r[0].pi[0],0,0),one,p->Fp);
    memset(req.vk_x[0],0xff,32);
    assert(!gamma_ov_sh_verify(&v,req.vk_x,r,p)); gamma_ov_sh_verify_out_clear(&v,p);
    for (int i=0;i<p->k;i++) gamma_ov_sh_compute_out_clear(&r[i],p);
    free(r); gamma_ov_sh_probgen_out_clear(&req,p); fq_clear(expected,p->Fp); fq_clear(one,p->Fp);
}
static void check_hss(fq_mat_t F, fq_mat_t x, pubpar *p, gamma_ov_hss_keygen_out *key) {
    gamma_ov_hss_probgen_out req;
    assert(gamma_ov_hss_probgen(&req,&key->pk,x,p));
    gamma_ov_hss_compute_out r[2], repeat;
    for (int b=0;b<2;b++) assert(gamma_ov_hss_compute(&r[b],b,&key->ek[b],&key->pk,F,req.ct_inputs,p));
    assert(gamma_ov_hss_compute(&repeat,0,&key->ek[0],&key->pk,F,req.ct_inputs,p));
    for (int j=0;j<=p->d;j++) assert(fq_equal(r[0].values[j],repeat.values[j],p->Fp));
    gamma_ov_hss_compute_out_clear(&repeat,p);
    gamma_ov_hss_verify_out v; gamma_ov_hss_reconstruct_out decoded;
    assert(gamma_ov_hss_verify(&v,req.vk_x,r,p));
    gamma_ov_hss_reconstruct(&decoded,v.aux_R_prime,fq_mat_entry(F,0,0),p);
    fq_t expected,one; Eval(expected,F,x,p); fq_init(one,p->Fp); fq_one(one,p->Fp);
    assert(fq_equal(fq_mat_entry(decoded.result,0,0),expected,p->Fp));
    gamma_ov_hss_verify_out_clear(&v,p); gamma_ov_hss_reconstruct_out_clear(&decoded,p);
    fq_add(r[0].values[0],r[0].values[0],one,p->Fp);
    assert(!gamma_ov_hss_verify(&v,req.vk_x,r,p)); gamma_ov_hss_verify_out_clear(&v,p);
    fq_sub(r[0].values[0],r[0].values[0],one,p->Fp);
    memset(req.vk_x[0],0xff,32);
    assert(!gamma_ov_hss_verify(&v,req.vk_x,r,p)); gamma_ov_hss_verify_out_clear(&v,p);
    fmpz_zero(req.ct_inputs[0].c[0]);
    assert(!gamma_ov_hss_compute(&repeat,0,&key->ek[0],&key->pk,F,req.ct_inputs,p));
    gamma_ov_hss_compute_out_clear(&repeat,p);
    fq_clear(expected,p->Fp); fq_clear(one,p->Fp);
    for (int b=0;b<2;b++) gamma_ov_hss_compute_out_clear(&r[b],p);
    gamma_ov_hss_probgen_out_clear(&req,p);
}
static void hss_tests(void) {
    pubpar p; assert(pv_params_init(&p,2,3,1));
    puts("Generating fresh >=3072-bit HSS key..."); fflush(stdout);
    gamma_ov_hss_keygen_out key; assert(gamma_ov_hss_keygen(&key,&p));
    fq_mat_t F,x;
    fq_mat_init(F,10,1,p.Fp); fq_mat_init(x,2,1,p.Fp);
    fq_set_ui(fq_mat_entry(F,3,0),1,p.Fp); /* x0^2 */
    fq_set_ui(fq_mat_entry(x,0,0),2,p.Fp);
    check_hss(F,x,&p,&key); puts("HSS x^2 regression passed."); fflush(stdout);
    example(F,x,&p,0); check_hss(F,x,&p,&key);
    puts("HSS dense mixed degree-3 polynomial passed."); fflush(stdout);
    example(F,x,&p,1); check_hss(F,x,&p,&key);
    puts("HSS full-field inputs and coefficients passed."); fflush(stdout);
    fq_mat_zero(F,p.Fp); check_hss(F,x,&p,&key);
    fq_set_ui(fq_mat_entry(F,0,0),7,p.Fp); check_hss(F,x,&p,&key);
    fq_mat_clear(F,p.Fp); fq_mat_clear(x,p.Fp);
    gamma_ov_hss_keygen_out_clear(&key); pv_params_clear(&p);
}
static void high_degree_hss(void) {
    pubpar p; assert(pv_params_init(&p,1,10,1));
    puts("Generating range-sized degree-10 HSS key..."); fflush(stdout);
    gamma_ov_hss_keygen_out key; assert(gamma_ov_hss_keygen(&key,&p));
    fq_mat_t F,x; fq_mat_init(F,11,1,p.Fp); fq_mat_init(x,1,1,p.Fp);
    fq_t expected,a,b,total; fq_init(a,p.Fp); fq_init(b,p.Fp); fq_init(total,p.Fp);
    fq_set_si(fq_mat_entry(x,0,0),-1,p.Fp);
    hss_ciphertext ct; hss_input(&ct,&key.pk,fq_mat_entry(x,0,0),p.Fp);
    for (int d=4;d<=10;d+=2) {
        fq_mat_zero(F,p.Fp); fq_set_ui(fq_mat_entry(F,d,0),1,p.Fp);
        Eval(expected,F,x,&p);
        assert(hss_eval(a,0,&key.pk,&key.ek[0],&ct,F,&p));
        assert(hss_eval(b,1,&key.pk,&key.ek[1],&ct,F,&p));
        fq_add(total,a,b,p.Fp); assert(fq_equal(total,expected,p.Fp));
        fq_clear(expected,p.Fp);
        printf("HSS degree %d at q-1 passed.\n",d); fflush(stdout);
    }
    hss_ciphertext_clear(&ct); fq_clear(a,p.Fp); fq_clear(b,p.Fp); fq_clear(total,p.Fp);
    fq_mat_clear(F,p.Fp); fq_mat_clear(x,p.Fp);
    gamma_ov_hss_keygen_out_clear(&key); pv_params_clear(&p);
}
int main(int argc, char **argv) {
    if (argc>1 && !strcmp(argv[1],"--hss")) { hss_tests(); return 0; }
    if (argc>1 && !strcmp(argv[1],"--hss-high")) { high_degree_hss(); return 0; }
    for (int d=1;d<=10;d++) for (int t=1;t<=3;t++) {
        pubpar p; assert(pv_params_init(&p,d<=3?2:1,d,t));
        fq_mat_t F,x; fq_mat_init(F,mpoly_num_coefficients(p.m,d),1,p.Fp);
        fq_mat_init(x,p.m,1,p.Fp);
        for (int trial=0;trial<2;trial++) { example(F,x,&p,trial); check_fe(F,x,&p); check_ov(F,x,&p); check_pi3(F,x,&p); }
        fq_mat_zero(F,p.Fp); check_fe(F,x,&p); check_ov(F,x,&p); check_pi3(F,x,&p);
        fq_set_ui(fq_mat_entry(F,0,0),7,p.Fp); check_fe(F,x,&p); check_ov(F,x,&p); check_pi3(F,x,&p);
        fq_mat_clear(F,p.Fp); fq_mat_clear(x,p.Fp); pv_params_clear(&p);
    }
    puts("FE/OV/Pi3 end-to-end, degree 1..10, threshold 1..3, tamper and malformed-key tests passed.");
    return 0;
}
