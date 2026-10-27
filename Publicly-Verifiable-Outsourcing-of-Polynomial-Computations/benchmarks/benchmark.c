#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <limits.h>
#include <time.h>
#include "gamma_fe_sh.h"
#include "gamma_ov_sh.h"
#include "gamma_ov_hss.h"

static double now_ms(void) {
    struct timespec ts;
    if(clock_gettime(CLOCK_MONOTONIC,&ts)) abort();
    return ts.tv_sec*1000.0+ts.tv_nsec/1000000.0;
}
/* Public deterministic workload fixtures; never used for masks or keys. */
static void fixture(fq_t out, const char *kind, int m, int d, slong index, int trial,
                    pubpar *p) {
    char label[128];
    int n=snprintf(label,sizeof(label),"pv-benchmark-v1:%s:%d:%d:%ld:%d",kind,m,d,(long)index,trial);
    unsigned char digest[64]; crypto_hash_sha512(digest,(unsigned char *)label,(unsigned long long)n);
    mpz_t z; mpz_init(z); mpz_import(z,64,-1,1,0,0,digest);
    fmpz_t value; fmpz_init(value); fmpz_set_mpz(value,z); fq_set_fmpz(out,value,p->Fp);
    fmpz_clear(value); mpz_clear(z);
}
static int number(const char *s) {
    char *end; errno=0; long n=strtol(s,&end,10);
    if(errno || *end || n<1 || n>INT_MAX) { fprintf(stderr,"Invalid positive integer: %s\n",s); exit(2); }
    return (int)n;
}
/* Canonical ASCII spellings of the manuscript's scheme names. */
static const char *scheme_names[] = { NULL, "Gamma_FE_Sh", "Gamma_OV_Sh", "Gamma_OV_HSS" };
static int parse_scheme(const char *name) {
    for (int i=1;i<=3;i++) {
        char legacy[2] = { (char)('0'+i), '\0' };
        if (!strcmp(name,scheme_names[i]) || !strcmp(name,legacy)) return i;
    }
    fprintf(stderr,"Unknown scheme: %s (use Gamma_FE_Sh, Gamma_OV_Sh, or Gamma_OV_HSS)\n",name);
    exit(2);
}
static void ensure(int valid) {
    if(!valid) { fputs("Correctness check failed; timing row discarded.\n",stderr); exit(1); }
}
int main(int argc,char **argv) {
    if(argc!=7) {
        fprintf(stderr,"Usage: %s SCHEME(Gamma_FE_Sh|Gamma_OV_Sh|Gamma_OV_HSS) M D T TRIALS dense\n",argv[0]); return 2;
    }
    int scheme=parse_scheme(argv[1]),m=number(argv[2]),d=number(argv[3]),t=number(argv[4]),trials=number(argv[5]);
    if(scheme>3 || (scheme==3 && t!=1) || strcmp(argv[6],"dense")) return 2;
    pubpar p; ensure(pv_params_init(&p,m,d,t));
    slong count=mpoly_num_coefficients(m,d); ensure(count>0);
    fq_mat_t F,x; fq_mat_init(F,count,1,p.Fp); fq_mat_init(x,m,1,p.Fp);
    for(slong i=0;i<count;i++) fixture(fq_mat_entry(F,i,0),"coefficient",m,d,i,0,&p);
    gamma_ov_hss_keygen_out key; double begin=now_ms(),keygen_ms=0;
    if(scheme==3) { ensure(gamma_ov_hss_keygen(&key,&p)); keygen_ms=now_ms()-begin; }
    fprintf(stderr,"FLINT=%s GMP=%s sodium=%s scheme=%s m=%d d=%d t=%d\n",FLINT_VERSION,gmp_version,sodium_version_string(),scheme_names[scheme],m,d,t);
    puts("workload,scheme,m,d,t,servers,coefficients,trial,native_inner_runs,keygen_ms,native_ms,probgen_ms,server_total_ms,verify_ms,decode_ms,client_ms,client_ratio,input_bytes_per_server,response_bytes_per_server,paillier_bits,correct");
    fflush(stdout);
    for(int trial=0;trial<trials;trial++) {
        for(int i=0;i<m;i++) fixture(fq_mat_entry(x,i,0),"input",m,d,i,trial,&p);
        fq_t expected; Eval(expected,F,x,&p);
        begin=now_ms();
        for(int i=0;i<100;i++) { fq_t temp; Eval(temp,F,x,&p); fq_clear(temp,p.Fp); }
        double native=(now_ms()-begin)/100.0,pg,server,verify,decode;
        size_t input_bytes=0,response_bytes=0;
        if(scheme==1) {
            gamma_fe_sh_probgen_out req; begin=now_ms(); ensure(gamma_fe_sh_probgen(&req,x,&p)); pg=now_ms()-begin;
            ext_field_elem_t *r=pv_alloc(p.k,sizeof(*r));
            begin=now_ms();
            for(int i=0;i<p.k;i++) { ext_field_elem_init(&r[i]); gamma_fe_sh_compute(&r[i],F,&req.shares[i],&req.masks[i],&req.vk.field,&p); }
            server=now_ms()-begin;
            gamma_fe_sh_verify_out v; begin=now_ms(); int ok=gamma_fe_sh_verify(&v,&req.vk,r,&p); verify=now_ms()-begin; ensure(ok);
            fq_t actual; fq_init(actual,p.Fp); begin=now_ms(); ok=gamma_fe_sh_decode(actual,&req,&v,&p); decode=now_ms()-begin;
            ensure(ok && fq_equal(actual,expected,p.Fp));
            for(int i=0;i<p.k;i++) ext_field_elem_clear(&r[i]); free(r);
            fq_clear(actual,p.Fp); gamma_fe_sh_verify_out_clear(&v,&p); gamma_fe_sh_probgen_out_clear(&req,&p);
            input_bytes=64*((size_t)m+1); response_bytes=64;
        } else if(scheme==2) {
            gamma_ov_sh_probgen_out req; begin=now_ms(); ensure(gamma_ov_sh_probgen(&req,x,&p)); pg=now_ms()-begin;
            gamma_ov_sh_compute_out *r=pv_alloc(p.k,sizeof(*r)); begin=now_ms();
            for(int i=0;i<p.k;i++) gamma_ov_sh_compute(&r[i],i,F,req.sigma[i],&p); server=now_ms()-begin;
            gamma_ov_sh_verify_out v; begin=now_ms(); int ok=gamma_ov_sh_verify(&v,req.vk_x,r,&p); verify=now_ms()-begin; ensure(ok);
            gamma_ov_sh_reconstruct_out decoded; begin=now_ms(); gamma_ov_sh_reconstruct(&decoded,v.aux_R_prime,fq_mat_entry(F,0,0),&p); decode=now_ms()-begin;
            ensure(fq_equal(fq_mat_entry(decoded.result,0,0),expected,p.Fp));
            for(int i=0;i<p.k;i++) gamma_ov_sh_compute_out_clear(&r[i],&p); free(r);
            gamma_ov_sh_verify_out_clear(&v,&p); gamma_ov_sh_reconstruct_out_clear(&decoded,&p); gamma_ov_sh_probgen_out_clear(&req,&p);
            input_bytes=32*(size_t)m*(d+1); response_bytes=32*((size_t)d+1);
        } else {
            gamma_ov_hss_probgen_out req; begin=now_ms(); ensure(gamma_ov_hss_probgen(&req,&key.pk,x,&p)); pg=now_ms()-begin;
            gamma_ov_hss_compute_out r[2]; begin=now_ms();
            for(int b=0;b<2;b++) ensure(gamma_ov_hss_compute(&r[b],b,&key.ek[b],&key.pk,F,req.ct_inputs,&p)); server=now_ms()-begin;
            gamma_ov_hss_verify_out v; begin=now_ms(); int ok=gamma_ov_hss_verify(&v,req.vk_x,r,&p); verify=now_ms()-begin; ensure(ok);
            gamma_ov_hss_reconstruct_out decoded; begin=now_ms(); gamma_ov_hss_reconstruct(&decoded,v.aux_R_prime,fq_mat_entry(F,0,0),&p); decode=now_ms()-begin;
            ensure(fq_equal(fq_mat_entry(decoded.result,0,0),expected,p.Fp));
            for(int b=0;b<2;b++) gamma_ov_hss_compute_out_clear(&r[b],&p);
            gamma_ov_hss_verify_out_clear(&v,&p); gamma_ov_hss_reconstruct_out_clear(&decoded,&p); gamma_ov_hss_probgen_out_clear(&req,&p);
            input_bytes=hss_input_bytes(&key.pk)*(size_t)m*(d+1); response_bytes=32*((size_t)d+1);
        }
        double client=pg+verify+decode;
        printf("%s,%s,%d,%d,%d,%d,%ld,%d,100,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%zu,%zu,%u,1\n",
               argv[6],scheme_names[scheme],m,d,t,scheme==3?2:p.k,(long)count,trial,keygen_ms,native,pg,server,verify,decode,client,client/native,input_bytes,response_bytes,scheme==3?key.pk.modulus_bits:0);
        fflush(stdout); fq_clear(expected,p.Fp);
    }
    if(scheme==3) gamma_ov_hss_keygen_out_clear(&key);
    fq_mat_clear(F,p.Fp); fq_mat_clear(x,p.Fp); pv_params_clear(&p); return 0;
}
