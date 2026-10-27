#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "scheme1.h"
#include "scheme2.h"
#include "scheme3.h"

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
    scheme1_probgen_out req;
    assert(scheme1_probgen(&req,x,p));
    ext_field_elem_t *r=pv_alloc(p->k,sizeof(*r));
    for (int i=0;i<p->k;i++) {
        ext_field_elem_init(&r[i]);
        scheme1_compute(&r[i],F,&req.shares[i],&req.masks[i],&req.vk.field,p);
    }
    scheme1_verify_out verified;
    assert(scheme1_verify(&verified,&req.vk,r,p));
    fq_t actual, expected;
    fq_init(actual,p->Fp); Eval(expected,F,x,p);
    assert(scheme1_decode(actual,&req,&verified,p));
    assert(fq_equal(actual,expected,p->Fp));
    scheme1_verify_out_clear(&verified,p);
    /* Imaginary constant changes fail public verification deterministically. */
    for (int i=0;i<p->k;i++) {
        fmpz_add_ui(r[i].b,r[i].b,1); fmpz_mod(r[i].b,r[i].b,p->p);
    }
    assert(!scheme1_verify(&verified,&req.vk,r,p));
    assert(!scheme1_decode(actual,&req,&verified,p));
    scheme1_verify_out_clear(&verified,p);
    for (int i=0;i<p->k;i++) {
        fmpz_sub_ui(r[i].b,r[i].b,1); fmpz_mod(r[i].b,r[i].b,p->p);
    }
    memset(req.vk.n[0],0xff,32);
    assert(!scheme1_verify(&verified,&req.vk,r,p));
    scheme1_verify_out_clear(&verified,p);
    for (int i=0;i<p->k;i++) ext_field_elem_clear(&r[i]);
    free(r); fq_clear(actual,p->Fp); fq_clear(expected,p->Fp);
    scheme1_probgen_out_clear(&req,p);
}
static void check_ov(fq_mat_t F, fq_mat_t x, pubpar *p) {
    scheme2_probgen_out req, next;
    assert(scheme2_probgen(&req,x,p)); assert(scheme2_probgen(&next,x,p));
    assert(!fq_mat_equal(req.sigma[0][0],next.sigma[0][0],p->Fp));
    scheme2_probgen_out_clear(&next,p);
    scheme2_compute_out *r=pv_alloc(p->k,sizeof(*r));
    for (int i=0;i<p->k;i++) scheme2_compute(&r[i],i,F,req.sigma[i],p);
    scheme2_verify_out v; scheme2_reconstruct_out decoded;
    assert(scheme2_verify(&v,req.vk_x,r,p));
    scheme2_reconstruct(&decoded,v.aux_R_prime,fq_mat_entry(F,0,0),p);
    fq_t expected, one; Eval(expected,F,x,p); fq_init(one,p->Fp); fq_one(one,p->Fp);
    assert(fq_equal(fq_mat_entry(decoded.result,0,0),expected,p->Fp));
    scheme2_verify_out_clear(&v,p); scheme2_reconstruct_out_clear(&decoded,p);
    fq_add(fq_mat_entry(r[0].pi[0],0,0),fq_mat_entry(r[0].pi[0],0,0),one,p->Fp);
    assert(!scheme2_verify(&v,req.vk_x,r,p)); scheme2_verify_out_clear(&v,p);
    fq_sub(fq_mat_entry(r[0].pi[0],0,0),fq_mat_entry(r[0].pi[0],0,0),one,p->Fp);
    memset(req.vk_x[0],0xff,32);
    assert(!scheme2_verify(&v,req.vk_x,r,p)); scheme2_verify_out_clear(&v,p);
    for (int i=0;i<p->k;i++) scheme2_compute_out_clear(&r[i],p);
    free(r); scheme2_probgen_out_clear(&req,p); fq_clear(expected,p->Fp); fq_clear(one,p->Fp);
}
static void check_hss(fq_mat_t F, fq_mat_t x, pubpar *p, scheme3_keygen_out *key) {
    scheme3_probgen_out req;
    assert(scheme3_probgen(&req,&key->pk,x,p));
    scheme3_compute_out r[2], repeat;
    for (int b=0;b<2;b++) assert(scheme3_compute(&r[b],b,&key->ek[b],&key->pk,F,req.ct_inputs,p));
    assert(scheme3_compute(&repeat,0,&key->ek[0],&key->pk,F,req.ct_inputs,p));
    for (int j=0;j<=p->d;j++) assert(fq_equal(r[0].values[j],repeat.values[j],p->Fp));
    scheme3_compute_out_clear(&repeat,p);
    scheme3_verify_out v; scheme3_reconstruct_out decoded;
    assert(scheme3_verify(&v,req.vk_x,r,p));
    scheme3_reconstruct(&decoded,v.aux_R_prime,fq_mat_entry(F,0,0),p);
    fq_t expected,one; Eval(expected,F,x,p); fq_init(one,p->Fp); fq_one(one,p->Fp);
    assert(fq_equal(fq_mat_entry(decoded.result,0,0),expected,p->Fp));
    scheme3_verify_out_clear(&v,p); scheme3_reconstruct_out_clear(&decoded,p);
    fq_add(r[0].values[0],r[0].values[0],one,p->Fp);
    assert(!scheme3_verify(&v,req.vk_x,r,p)); scheme3_verify_out_clear(&v,p);
    fq_sub(r[0].values[0],r[0].values[0],one,p->Fp);
    memset(req.vk_x[0],0xff,32);
    assert(!scheme3_verify(&v,req.vk_x,r,p)); scheme3_verify_out_clear(&v,p);
    fmpz_zero(req.ct_inputs[0].c[0]);
    assert(!scheme3_compute(&repeat,0,&key->ek[0],&key->pk,F,req.ct_inputs,p));
    scheme3_compute_out_clear(&repeat,p);
    fq_clear(expected,p->Fp); fq_clear(one,p->Fp);
    for (int b=0;b<2;b++) scheme3_compute_out_clear(&r[b],p);
    scheme3_probgen_out_clear(&req,p);
}
static void hss_tests(void) {
    pubpar p; assert(pv_params_init(&p,2,3,1));
    puts("Generating fresh >=3072-bit HSS key..."); fflush(stdout);
    scheme3_keygen_out key; assert(scheme3_keygen(&key,&p));
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
    scheme3_keygen_out_clear(&key); pv_params_clear(&p);
}
static void high_degree_hss(void) {
    pubpar p; assert(pv_params_init(&p,1,10,1));
    puts("Generating range-sized degree-10 HSS key..."); fflush(stdout);
    scheme3_keygen_out key; assert(scheme3_keygen(&key,&p));
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
    scheme3_keygen_out_clear(&key); pv_params_clear(&p);
}
int main(int argc, char **argv) {
    if (argc>1 && !strcmp(argv[1],"--hss")) { hss_tests(); return 0; }
    if (argc>1 && !strcmp(argv[1],"--hss-high")) { high_degree_hss(); return 0; }
    for (int d=1;d<=10;d++) for (int t=1;t<=3;t++) {
        pubpar p; assert(pv_params_init(&p,d<=3?2:1,d,t));
        fq_mat_t F,x; fq_mat_init(F,mpoly_num_coefficients(p.m,d),1,p.Fp);
        fq_mat_init(x,p.m,1,p.Fp);
        for (int trial=0;trial<2;trial++) { example(F,x,&p,trial); check_fe(F,x,&p); check_ov(F,x,&p); }
        fq_mat_zero(F,p.Fp); check_fe(F,x,&p); check_ov(F,x,&p);
        fq_set_ui(fq_mat_entry(F,0,0),7,p.Fp); check_fe(F,x,&p); check_ov(F,x,&p);
        fq_mat_clear(F,p.Fp); fq_mat_clear(x,p.Fp); pv_params_clear(&p);
    }
    puts("FE/OV end-to-end, degree 1..10, threshold 1..3, tamper and malformed-key tests passed.");
    return 0;
}
