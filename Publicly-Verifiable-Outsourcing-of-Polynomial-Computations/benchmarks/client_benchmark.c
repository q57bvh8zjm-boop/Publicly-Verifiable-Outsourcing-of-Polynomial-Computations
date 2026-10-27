#define _POSIX_C_SOURCE 200809L
#include <time.h>
#include "gamma_fe_sh.h"
#include "gamma_ov_sh.h"
#include "gamma_ov_hss.h"

/* Figure 5 CLIENT microbenchmark. Server responses are valid reference fixtures,
 * prepared outside every measured interval. No server or end-to-end time is
 * reported. Production ProbGen/Verify/Decode routines are called unchanged. */
static double tick(void) {
    struct timespec t; if (clock_gettime(CLOCK_MONOTONIC, &t)) abort();
    return (double)t.tv_sec * 1000.0 + (double)t.tv_nsec / 1000000.0;
}
static void check(int ok) { if (!ok) { fputs("client benchmark check failed\n",stderr); exit(1); } }
static void fixture(fq_t out, const char *kind, int m, int d, slong i, int trial, pubpar *p) {
    char label[128]; int n=snprintf(label,sizeof(label),"pv-benchmark-v1:%s:%d:%d:%ld:%d",kind,m,d,(long)i,trial);
    unsigned char digest[64]; crypto_hash_sha512(digest,(unsigned char *)label,n);
    mpz_t z; mpz_init(z); mpz_import(z,64,-1,1,0,0,digest);
    fmpz_t a; fmpz_init(a); fmpz_set_mpz(a,z); fq_set_fmpz(out,a,p->Fp);
    fmpz_clear(a); mpz_clear(z);
}
static void emit(const char *scheme,pubpar *p,int trial,slong count,int nr,double native,
                 double pg,double verify,double decode,double keygen,unsigned bits,
                 unsigned digits,const char *source) {
    double client=pg+verify+decode;
    printf("%s,%d,%d,%d,%ld,%d,%d,25,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%u,%u,%s,1\n",
        scheme,p->m,p->d,p->t,(long)count,trial,nr,native,pg,verify,decode,client,client/native,keygen,bits,digits,source);
    fflush(stdout);
}
static void scaled_eval(fq_t out,fq_mat_t F,fq_mat_t x,const fq_t v,pubpar *p) {
    fq_mat_t scaled; fq_mat_init(scaled,p->m,1,p->Fp);
    for(int i=0;i<p->m;i++) fq_mul(fq_mat_entry(scaled,i,0),fq_mat_entry(x,i,0),v,p->Fp);
    fq_t y; Eval(y,F,scaled,p); fq_sub(out,y,fq_mat_entry(F,0,0),p->Fp);
    fq_clear(y,p->Fp); fq_mat_clear(scaled,p->Fp);
}
static void fe(pubpar *p,fq_mat_t x,const fq_t expected,int trial,slong count,int nr,double native) {
    gamma_fe_sh_probgen_out req; double start=tick(); check(gamma_fe_sh_probgen(&req,x,p)); double pg=tick()-start;
    /* F(x) plus the random degree-dt masking polynomial: a full-field valid
     * response vector with the same public verification and decoding work. */
    ext_field_elem_t *responses=pv_alloc(p->k,sizeof(*responses)); fmpz_t y; fmpz_init(y); fq2fmpz(y,expected,p->Fp);
    for(int i=0;i<p->k;i++) { ext_field_elem_init(&responses[i]);
        fmpz_add(responses[i].a,req.masks[i].a,y); fmpz_mod(responses[i].a,responses[i].a,p->p);
        fmpz_set(responses[i].b,req.masks[i].b); }
    double verify=0,decode=0;
    for(int r=-1;r<25;r++) { gamma_fe_sh_verify_out v; start=tick(); int ok=gamma_fe_sh_verify(&v,&req.vk,responses,p); double elapsed=tick()-start; check(ok); if(r>=0)verify+=elapsed;
        fq_t actual; fq_init(actual,p->Fp); start=tick(); ok=gamma_fe_sh_decode(actual,&req,&v,p); elapsed=tick()-start;
        check(ok&&fq_equal(actual,expected,p->Fp)); if(r>=0)decode+=elapsed;
        fq_clear(actual,p->Fp); gamma_fe_sh_verify_out_clear(&v,p); }
    if(trial>=0)emit("Gamma_FE_Sh",p,trial,count,nr,native,pg,verify/25,decode/25,0,0,0,"masked_reference_value");
    for(int i=0;i<p->k;i++)ext_field_elem_clear(&responses[i]); free(responses); fmpz_clear(y); gamma_fe_sh_probgen_out_clear(&req,p);
}
static void ov(pubpar *p,fq_mat_t F,fq_mat_t x,const fq_t expected,int trial,slong count,int nr,double native) {
    gamma_ov_sh_probgen_out req; double start=tick(); check(gamma_ov_sh_probgen(&req,x,p)); double pg=tick()-start;
    gamma_ov_sh_compute_out *responses=pv_alloc(p->k,sizeof(*responses));
    for(int i=0;i<p->k;i++) { responses[i].pi=pv_alloc(p->d+1,sizeof(fq_mat_t)); for(int j=0;j<=p->d;j++)fq_mat_init(responses[i].pi[j],1,1,p->Fp); }
    fq_t y,value; fq_init(y,p->Fp); fq_init(value,p->Fp);
    fq_t *coeff=pv_alloc(p->k-1,sizeof(fq_t)); for(int s=0;s<p->k-1;s++)fq_init(coeff[s],p->Fp);
    for(int j=0;j<=p->d;j++) {
        if(!j)fq_sub(y,expected,fq_mat_entry(F,0,0),p->Fp); else scaled_eval(y,F,x,req.v_values[j-1],p);
        for(int s=0;s<p->k-1;s++)pv_random_field(coeff[s],p->Fp);
        for(int i=0;i<p->k;i++) { fq_set(value,coeff[p->k-2],p->Fp);
            for(int s=p->k-3;s>=0;s--) { fq_mul(value,value,p->SID[i],p->Fp); fq_add(value,value,coeff[s],p->Fp); }
            fq_mul(value,value,p->SID[i],p->Fp); fq_add(fq_mat_entry(responses[i].pi[j],0,0),value,y,p->Fp); }
    }
    double verify=0,decode=0;
    for(int r=-1;r<25;r++) { gamma_ov_sh_verify_out v; start=tick(); int ok=gamma_ov_sh_verify(&v,req.vk_x,responses,p); double elapsed=tick()-start; check(ok); if(r>=0)verify+=elapsed;
        gamma_ov_sh_reconstruct_out actual; start=tick(); gamma_ov_sh_reconstruct(&actual,v.aux_R_prime,fq_mat_entry(F,0,0),p); elapsed=tick()-start;
        check(fq_equal(fq_mat_entry(actual.result,0,0),expected,p->Fp)); if(r>=0)decode+=elapsed;
        gamma_ov_sh_reconstruct_out_clear(&actual,p); gamma_ov_sh_verify_out_clear(&v,p); }
    if(trial>=0)emit("Gamma_OV_Sh",p,trial,count,nr,native,pg,verify/25,decode/25,0,0,0,"reference_eval_shamir_split");
    for(int s=0;s<p->k-1;s++)fq_clear(coeff[s],p->Fp); free(coeff); fq_clear(y,p->Fp); fq_clear(value,p->Fp);
    for(int i=0;i<p->k;i++)gamma_ov_sh_compute_out_clear(&responses[i],p); free(responses); gamma_ov_sh_probgen_out_clear(&req,p);
}
static void hss(pubpar *p,fq_mat_t F,fq_mat_t x,const fq_t expected,gamma_ov_hss_keygen_out *key,
                const fmpz_t secret,int trial,slong count,int nr,double native,double keygen) {
    gamma_ov_hss_probgen_out req; double start=tick(); check(gamma_ov_hss_probgen(&req,&key->pk,x,p)); double pg=tick()-start;
    gamma_ov_hss_compute_out responses[2];
    for(int b=0;b<2;b++) { responses[b].valid=1; responses[b].values=pv_alloc(p->d+1,sizeof(fq_t)); for(int j=0;j<=p->d;j++)fq_init(responses[b].values[j],p->Fp); }
    /* BENCHMARK ONLY: the harness holds both evaluation keys, recovers the
     * Paillier secret, and decrypts the first coordinate to recover each v_j.
     * This builds reference responses, never a production server shortcut. */
    fmpz_t plaintext; fmpz_init(plaintext); fq_t y,v,inverse; fq_init(y,p->Fp); fq_init(v,p->Fp); fq_init(inverse,p->Fp);
    check(!fq_is_zero(fq_mat_entry(x,0,0),p->Fp)); fq_inv(inverse,fq_mat_entry(x,0,0),p->Fp);
    for(int j=0;j<=p->d;j++) {
        fmpz_powm(plaintext,req.ct_inputs[j*p->m].c[0],secret,key->pk.n2);
        fmpz_sub_ui(plaintext,plaintext,1); check(fmpz_divisible(plaintext,key->pk.n)); fmpz_divexact(plaintext,plaintext,key->pk.n);
        fq_set_fmpz(v,plaintext,p->Fp); fq_mul(v,v,inverse,p->Fp);
        if(!j) { check(fq_is_one(v,p->Fp)); fq_sub(y,expected,fq_mat_entry(F,0,0),p->Fp); } else scaled_eval(y,F,x,v,p);
        pv_random_field(responses[0].values[j],p->Fp); fq_sub(responses[1].values[j],y,responses[0].values[j],p->Fp);
    }
    double verify=0,decode=0;
    for(int r=-1;r<25;r++) { gamma_ov_hss_verify_out vout; start=tick(); int ok=gamma_ov_hss_verify(&vout,req.vk_x,responses,p); double elapsed=tick()-start; check(ok); if(r>=0)verify+=elapsed;
        gamma_ov_hss_reconstruct_out actual; start=tick(); gamma_ov_hss_reconstruct(&actual,vout.aux_R_prime,fq_mat_entry(F,0,0),p); elapsed=tick()-start;
        check(fq_equal(fq_mat_entry(actual.result,0,0),expected,p->Fp)); if(r>=0)decode+=elapsed;
        gamma_ov_hss_reconstruct_out_clear(&actual,p); gamma_ov_hss_verify_out_clear(&vout,p); }
    if(trial>=0)emit("Gamma_OV_HSS",p,trial,count,nr,native,pg,verify/25,decode/25,keygen,key->pk.modulus_bits,key->pk.digit_count,"reference_eval_additive_split");
    for(int b=0;b<2;b++)gamma_ov_hss_compute_out_clear(&responses[b],p);
    fq_clear(y,p->Fp); fq_clear(v,p->Fp); fq_clear(inverse,p->Fp); fmpz_clear(plaintext); gamma_ov_hss_probgen_out_clear(&req,p);
}
int main(int argc,char **argv) {
    if(argc!=5) { fprintf(stderr,"Usage: %s DEGREE MAX_M TRIALS WITH_HSS(0|1)\n",argv[0]); return 2; }
    int d=atoi(argv[1]),max_m=atoi(argv[2]),trials=atoi(argv[3]),with_hss=atoi(argv[4]);
    check(d>=1&&d<=10&&max_m>=1&&max_m<=15&&trials>=1&&trials<=100&&(with_hss==0||with_hss==1));
    gamma_ov_hss_keygen_out key; double keygen=0; fmpz_t secret,delta; fmpz_init(secret); fmpz_init(delta);
    if(with_hss) { pubpar maxpar; check(pv_params_init(&maxpar,max_m,d,1)); double start=tick();
        check(gamma_ov_hss_keygen(&key,&maxpar)); keygen=tick()-start;
        for(unsigned i=key.pk.digit_count;i>0;i--) { fmpz_mul_2exp(secret,secret,key.pk.base_bits); fmpz_sub(delta,key.ek[1].digits[i-1],key.ek[0].digits[i-1]); fmpz_add(secret,secret,delta); }
        pv_params_clear(&maxpar); }
    fprintf(stderr,"FLINT=%s GMP=%s sodium=%s degree=%d max_m=%d trials=%d keygen_ms=%.3f\n",FLINT_VERSION,gmp_version,sodium_version_string(),d,max_m,trials,keygen);
    puts("scheme,m,d,t,coefficients,trial,native_inner_runs,verify_decode_inner_runs,native_ms,probgen_ms,verify_ms,decode_ms,client_ms,client_ratio,keygen_ms,paillier_bits,paillier_digits,response_source,correct"); fflush(stdout);
    for(int m=1;m<=max_m;m++) {
        pubpar p; check(pv_params_init(&p,m,d,1)); slong count=mpoly_num_coefficients(m,d); check(count>0);
        fq_mat_t F,x; fq_mat_init(F,count,1,p.Fp); fq_mat_init(x,m,1,p.Fp);
        for(slong i=0;i<count;i++)fixture(fq_mat_entry(F,i,0),"coefficient",m,d,i,0,&p);
        for(int trial=-1;trial<trials;trial++) {
            for(int i=0;i<m;i++)fixture(fq_mat_entry(x,i,0),"input",m,d,i,trial,&p);
            double start=tick(); fq_t expected; Eval(expected,F,x,&p); double single=tick()-start;
            int nr=single>0?(int)(20.0/single):10000; if(nr<1)nr=1; if(nr>10000)nr=10000;
            start=tick(); for(int r=0;r<nr;r++) { fq_t temp; Eval(temp,F,x,&p); fq_clear(temp,p.Fp); }
            double native=(tick()-start)/nr;
            fe(&p,x,expected,trial,count,nr,native); ov(&p,F,x,expected,trial,count,nr,native);
            if(with_hss)hss(&p,F,x,expected,&key,secret,trial,count,nr,native,keygen);
            fq_clear(expected,p.Fp);
        }
        fprintf(stderr,"complete d=%d m=%d\n",d,m); fflush(stderr);
        fq_mat_clear(F,p.Fp); fq_mat_clear(x,p.Fp); pv_params_clear(&p);
    }
    if(with_hss)gamma_ov_hss_keygen_out_clear(&key); fmpz_clear(secret); fmpz_clear(delta); return 0;
}
