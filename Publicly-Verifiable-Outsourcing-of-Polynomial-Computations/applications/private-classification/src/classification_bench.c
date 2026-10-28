/* Standalone workload adapters; protocol setup, encoding, verification and
 * decoding are linked unchanged from the main implementation. */
#define _POSIX_C_SOURCE 200809L
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>
#include "gamma_fe_sh.h"
#include "gamma_ov_sh.h"

typedef __int128 i128;
typedef struct {
    int kind,m,d,n,nq;
    int *a,*b,*label;
    int64_t *coef,*V,*X,offset;
    double *dc,*dv;
    i128 bias;
    fmpz_t fbias,fzero;
} Model;
static volatile i128 integer_sink;
static volatile double floating_sink;
static volatile unsigned char wire_sink;
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec*1000.0+t.tv_nsec/1e6;}
static i128 parse128(const char *s){int neg=(*s=='-');if(neg)s++;i128 v=0;while(*s)v=v*10+(*s++-'0');return neg?-v:v;}
static void str128(char *s,i128 v){char b[64];int n=0,neg=v<0;if(neg)v=-v;do{b[n++]='0'+v%10;v/=10;}while(v);int p=0;if(neg)s[p++]='-';while(n)s[p++]=b[--n];s[p]=0;}
static void set128(fmpz_t z,i128 v){char b[64];str128(b,v);assert(!fmpz_set_str(z,b,10));}
static Model read_model(const char *path){
    Model M={0};FILE *f=fopen(path,"r");assert(f);
    assert(fscanf(f,"%d%d%d%d%d",&M.kind,&M.m,&M.d,&M.n,&M.nq)==5);
    M.coef=calloc(M.n,sizeof(int64_t));M.dc=calloc(M.n,sizeof(double));
    M.X=calloc((size_t)M.nq*M.m,sizeof(int64_t));M.label=calloc(M.nq,sizeof(int));
    fmpz_init(M.fbias);fmpz_init(M.fzero);
    if(!M.kind){
        M.a=calloc(M.n,sizeof(int));M.b=calloc(M.n,sizeof(int));
        for(int j=0;j<M.n;j++){assert(fscanf(f,"%"SCNd64"%d%d",&M.coef[j],&M.a[j],&M.b[j])==3);M.dc[j]=M.coef[j];}
        fmpz_set_si(M.fzero,M.coef[0]);
    }else{
        char s[128];assert(fscanf(f,"%"SCNd64"%127s",&M.offset,s)==2);
        M.bias=parse128(s);set128(M.fbias,M.bias);
        M.V=calloc((size_t)M.n*M.m,sizeof(int64_t));M.dv=calloc((size_t)M.n*M.m,sizeof(double));
        i128 f0=M.bias;
        for(int j=0;j<M.n;j++){
            assert(fscanf(f,"%"SCNd64,&M.coef[j])==1);M.dc[j]=M.coef[j];
            for(int i=0;i<M.m;i++){size_t r=(size_t)j*M.m+i;assert(fscanf(f,"%"SCNd64,&M.V[r])==1);M.dv[r]=M.V[r];}
            i128 t=1;for(int d=0;d<M.d;d++)t*=M.offset;f0+=t*M.coef[j];
        }
        set128(M.fzero,f0);
    }
    for(int j=0;j<M.nq;j++){assert(fscanf(f,"%d",&M.label[j])==1);for(int i=0;i<M.m;i++)assert(fscanf(f,"%"SCNd64,&M.X[(size_t)j*M.m+i])==1);}
    fclose(f);return M;
}

/* Exact native integer baseline; prepare.py checks a <2^126 absolute bound. */
__attribute__((noinline)) static i128 native_integer(const Model *M,const int64_t *x){
    i128 sum=M->kind?M->bias:0;
    if(!M->kind){
        for(int j=0;j<M->n;j++)sum+=(i128)M->coef[j]*(M->a[j]<0?1:x[M->a[j]])*(M->b[j]<0?1:x[M->b[j]]);
    }else{
        for(int j=0;j<M->n;j++){
            int64_t dot=M->offset;const int64_t *v=M->V+(size_t)j*M->m;
            for(int i=0;i<M->m;i++)dot+=v[i]*x[i];
            i128 p=1;for(int d=0;d<M->d;d++)p*=dot;sum+=M->coef[j]*p;
        }
    }
    return sum;
}
__attribute__((noinline)) static double native_float(const Model *M,const double *x){
    double sum=M->kind?(double)M->bias:0;
    if(!M->kind){
        for(int j=0;j<M->n;j++)sum+=M->dc[j]*(M->a[j]<0?1:x[M->a[j]])*(M->b[j]<0?1:x[M->b[j]]);
    }else{
        for(int j=0;j<M->n;j++){
            double dot=M->offset;const double *v=M->dv+(size_t)j*M->m;
            for(int i=0;i<M->m;i++)dot+=v[i]*x[i];
            double p=1;for(int d=0;d<M->d;d++)p*=dot;sum+=M->dc[j]*p;
        }
    }
    return sum;
}

static void eval_base(fmpz_t out,const Model *M,fmpz_t *x,const fmpz_t q){
    fmpz_t dot,p,term;fmpz_init(dot);fmpz_init(p);fmpz_init(term);fmpz_zero(out);
    if(!M->kind){
        for(int j=0;j<M->n;j++){
            if(!M->coef[j])continue;
            fmpz_set_si(term,M->coef[j]);
            if(M->a[j]>=0){fmpz_mul(term,term,x[M->a[j]]);fmpz_mod(term,term,q);}
            if(M->b[j]>=0){fmpz_mul(term,term,x[M->b[j]]);fmpz_mod(term,term,q);}
            fmpz_add(out,out,term);
        }
    }else{
        fmpz_set(out,M->fbias);
        for(int j=0;j<M->n;j++){
            if(!M->coef[j])continue;
            fmpz_set_si(dot,M->offset);
            for(int i=0;i<M->m;i++)fmpz_addmul_si(dot,x[i],M->V[(size_t)j*M->m+i]);
            fmpz_mod(dot,dot,q);fmpz_powm_ui(p,dot,M->d,q);fmpz_addmul_si(out,p,M->coef[j]);
        }
    }
    fmpz_mod(out,out,q);fmpz_clear(dot);fmpz_clear(p);fmpz_clear(term);
}
static void eval_ext(ext_field_elem_t *out,const Model *M,ext_field_mat_t *x,ext_field_ctx_t *ctx){
    ext_field_elem_t dot,p,term;ext_field_elem_init(&dot);ext_field_elem_init(&p);ext_field_elem_init(&term);ext_field_elem_zero(out);
    if(!M->kind){
        for(int j=0;j<M->n;j++){
            if(!M->coef[j])continue;
            fmpz_set_si(term.a,M->coef[j]);fmpz_zero(term.b);
            if(M->a[j]>=0)ext_field_elem_mul(&term,&term,ext_field_mat_entry(x,M->a[j],0),ctx);
            if(M->b[j]>=0)ext_field_elem_mul(&term,&term,ext_field_mat_entry(x,M->b[j],0),ctx);
            fmpz_add(out->a,out->a,term.a);fmpz_add(out->b,out->b,term.b);
        }
    }else{
        fmpz_set(out->a,M->fbias);
        for(int j=0;j<M->n;j++){
            if(!M->coef[j])continue;
            fmpz_set_si(dot.a,M->offset);fmpz_zero(dot.b);
            for(int i=0;i<M->m;i++){
                const ext_field_elem_t *e=ext_field_mat_entry(x,i,0);slong c=M->V[(size_t)j*M->m+i];
                fmpz_addmul_si(dot.a,e->a,c);fmpz_addmul_si(dot.b,e->b,c);
            }
            fmpz_mod(dot.a,dot.a,ctx->p);fmpz_mod(dot.b,dot.b,ctx->p);ext_field_elem_one(&p);
            for(int d=0;d<M->d;d++)ext_field_elem_mul(&p,&p,&dot,ctx);
            fmpz_addmul_si(out->a,p.a,M->coef[j]);fmpz_addmul_si(out->b,p.b,M->coef[j]);
        }
    }
    fmpz_mod(out->a,out->a,ctx->p);fmpz_mod(out->b,out->b,ctx->p);
    ext_field_elem_clear(&dot);ext_field_elem_clear(&p);ext_field_elem_clear(&term);
}
static void wire_write(unsigned char *out,const fmpz_t z){
    mpz_t v;mpz_init(v);fmpz_get_mpz(v,z);size_t n=0;memset(out,0,32);mpz_export(out,&n,-1,1,0,0,v);assert(n<=32);mpz_clear(v);
}
static int check(fq_t actual,i128 want,pubpar *p){
    fmpz_t expected,got;fmpz_init(expected);fmpz_init(got);set128(expected,want);fmpz_mod(expected,expected,p->p);fq2fmpz(got,actual,p->Fp);
    int ok=fmpz_equal(expected,got);fmpz_clear(expected);fmpz_clear(got);return ok;
}

int main(int argc,char **argv){
    if(argc!=6){fprintf(stderr,"usage: %s fixture scheme[1|2] t trials predictions.tsv\n",argv[0]);return 2;}
    Model M=read_model(argv[1]);int scheme=atoi(argv[2]),t=atoi(argv[3]),trials=atoi(argv[4]);assert(scheme==1||scheme==2);
    pubpar p;assert(pv_params_init(&p,M.m,M.d,t));
    FILE *pred=fopen(argv[5],"w");assert(pred);fprintf(pred,"query\tlabel\tinteger_score\n");
    for(int i=0;i<M.nq;i++){char b[64];str128(b,native_integer(&M,M.X+(size_t)i*M.m));fprintf(pred,"%d\t%d\t%s\n",i,M.label[i],b);}fclose(pred);
    fprintf(stderr,"FLINT=%s GMP=%s sodium=%s; model=%s; complete native predictions=%d\n",FLINT_VERSION,gmp_version,sodium_version_string(),argv[1],M.nq);
    puts("scheme,t,m,d,servers,trial,native_int_ms,native_float_ms,native_field_ms,input_convert_ms,probgen_ms,serialize_ms,server_sum_ms,server_max_ms,verify_ms,decode_ms,cleanup_ms,client_ms,request_bytes,response_bytes,vk_bytes,correct,tamper_tests,tamper_rejections");fflush(stdout);
    fq_mat_t x;fq_mat_init(x,M.m,1,p.Fp);fq_t constant;fq_init(constant,p.Fp);fq_set_fmpz(constant,M.fzero,p.Fp);
    fmpz_t tmp,oracle;fmpz_init(tmp);fmpz_init(oracle);
    fmpz_t *xb=calloc(M.m,sizeof(fmpz_t));for(int i=0;i<M.m;i++)fmpz_init(xb[i]);
    double *xd=calloc(M.m,sizeof(double));
    size_t request_bytes=scheme==1?(size_t)p.k*(M.m+1)*64:(size_t)p.k*(M.d+1)*M.m*32;
    size_t response_bytes=scheme==1?(size_t)p.k*64:(size_t)p.k*(M.d+1)*32;
    size_t vk_bytes=scheme==1?(size_t)p.k*64:(size_t)(M.d+1)*32;
    unsigned char *wire=malloc(request_bytes);assert(wire);
    for(int trial=-1;trial<trials;trial++){
        int qi=trial<0?M.nq-1:trial%M.nq;int64_t *xi=M.X+(size_t)qi*M.m;
        for(int i=0;i<M.m;i++)xd[i]=xi[i];
        i128 want=native_integer(&M,xi);double begin=now();
        for(int r=0;r<1000;r++){__asm__ volatile("" ::: "memory");integer_sink=native_integer(&M,xi);}
        double native_int=(now()-begin)/1000;begin=now();
        for(int r=0;r<1000;r++){__asm__ volatile("" ::: "memory");floating_sink=native_float(&M,xd);}
        double native_fp=(now()-begin)/1000;begin=now();
        for(int i=0;i<M.m;i++){fmpz_set_si(xb[i],xi[i]);fq_set_fmpz(fq_mat_entry(x,i,0),xb[i],p.Fp);}
        double convert=now()-begin;begin=now();
        for(int r=0;r<5;r++)eval_base(oracle,&M,xb,p.p);
        double native_field=(now()-begin)/5;set128(tmp,want);fmpz_mod(tmp,tmp,p.p);assert(fmpz_equal(oracle,tmp));
        double pg,ser,servers=0,slowest=0,verify,decode,cleanup;int ntamper=0,rejected=0;
        if(scheme==1){
            gamma_fe_sh_probgen_out req;begin=now();assert(gamma_fe_sh_probgen(&req,x,&p));pg=now()-begin;
            begin=now();size_t pos=0;
            for(int i=0;i<p.k;i++){
                for(int j=0;j<M.m;j++){ext_field_elem_t *v=ext_field_mat_entry(&req.shares[i],j,0);wire_write(wire+pos,v->a);wire_write(wire+pos+32,v->b);pos+=64;}
                wire_write(wire+pos,req.masks[i].a);wire_write(wire+pos+32,req.masks[i].b);pos+=64;
            }
            assert(pos==request_bytes);wire_sink=wire[pos-1];ser=now()-begin;
            ext_field_elem_t *res=calloc(p.k,sizeof(*res));
            for(int i=0;i<p.k;i++){
                ext_field_elem_init(&res[i]);begin=now();eval_ext(&res[i],&M,&req.shares[i],&req.vk.field);ext_field_elem_add(&res[i],&res[i],&req.masks[i],&req.vk.field);
                double tm=now()-begin;servers+=tm;if(tm>slowest)slowest=tm;
            }
            gamma_fe_sh_verify_out v;begin=now();int ok=gamma_fe_sh_verify(&v,&req.vk,res,&p);verify=now()-begin;assert(ok);
            fq_t out;fq_init(out,p.Fp);begin=now();ok=gamma_fe_sh_decode(out,&req,&v,&p);decode=now()-begin;assert(ok&&check(out,want,&p));
            if(trial>=0&&trial<3){
                for(int i=0;i<p.k;i++){
                    fmpz_set(tmp,res[i].a);fmpz_add_ui(res[i].a,res[i].a,1);fmpz_mod(res[i].a,res[i].a,p.p);
                    gamma_fe_sh_verify_out bad;int accept=gamma_fe_sh_verify(&bad,&req.vk,res,&p);ntamper++;rejected+=!accept;assert(!accept);gamma_fe_sh_verify_out_clear(&bad,&p);fmpz_set(res[i].a,tmp);
                }
            }
            begin=now();fq_clear(out,p.Fp);gamma_fe_sh_verify_out_clear(&v,&p);gamma_fe_sh_probgen_out_clear(&req,&p);cleanup=now()-begin;
            for(int i=0;i<p.k;i++)ext_field_elem_clear(&res[i]);free(res);
        }else{
            gamma_ov_sh_probgen_out req;begin=now();assert(gamma_ov_sh_probgen(&req,x,&p));pg=now()-begin;
            begin=now();size_t pos=0;
            for(int i=0;i<p.k;i++)for(int j=0;j<=M.d;j++)for(int h=0;h<M.m;h++){
                fq2fmpz(tmp,fq_mat_entry(req.sigma[i][j],h,0),p.Fp);wire_write(wire+pos,tmp);pos+=32;
            }
            assert(pos==request_bytes);wire_sink=wire[pos-1];ser=now()-begin;
            gamma_ov_sh_compute_out *res=calloc(p.k,sizeof(*res));
            for(int i=0;i<p.k;i++){
                begin=now();res[i].pi=calloc(M.d+1,sizeof(fq_mat_t));
                for(int j=0;j<=M.d;j++){
                    fq_mat_init(res[i].pi[j],1,1,p.Fp);
                    for(int h=0;h<M.m;h++)fq2fmpz(xb[h],fq_mat_entry(req.sigma[i][j],h,0),p.Fp);
                    eval_base(oracle,&M,xb,p.p);fmpz_sub(oracle,oracle,M.fzero);fmpz_mod(oracle,oracle,p.p);fq_set_fmpz(fq_mat_entry(res[i].pi[j],0,0),oracle,p.Fp);
                }
                double tm=now()-begin;servers+=tm;if(tm>slowest)slowest=tm;
            }
            gamma_ov_sh_verify_out v;begin=now();int ok=gamma_ov_sh_verify(&v,req.vk_x,res,&p);verify=now()-begin;assert(ok);
            gamma_ov_sh_reconstruct_out out;begin=now();gamma_ov_sh_reconstruct(&out,v.aux_R_prime,constant,&p);decode=now()-begin;assert(check(fq_mat_entry(out.result,0,0),want,&p));
            if(trial>=0&&trial<3){
                for(int i=0;i<p.k;i++)for(int j=0;j<=M.d;j++){
                    fq_t saved,one;fq_init(saved,p.Fp);fq_init(one,p.Fp);fq_one(one,p.Fp);fq_set(saved,fq_mat_entry(res[i].pi[j],0,0),p.Fp);fq_add(fq_mat_entry(res[i].pi[j],0,0),saved,one,p.Fp);
                    gamma_ov_sh_verify_out bad;int accept=gamma_ov_sh_verify(&bad,req.vk_x,res,&p);ntamper++;rejected+=!accept;assert(!accept);gamma_ov_sh_verify_out_clear(&bad,&p);fq_set(fq_mat_entry(res[i].pi[j],0,0),saved,p.Fp);fq_clear(saved,p.Fp);fq_clear(one,p.Fp);
                }
            }
            begin=now();gamma_ov_sh_reconstruct_out_clear(&out,&p);gamma_ov_sh_verify_out_clear(&v,&p);gamma_ov_sh_probgen_out_clear(&req,&p);cleanup=now()-begin;
            for(int i=0;i<p.k;i++)gamma_ov_sh_compute_out_clear(&res[i],&p);free(res);
        }
        if(trial>=0){printf("%s,%d,%d,%d,%d,%d,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%zu,%zu,%zu,1,%d,%d\n",scheme==1?"Gamma_FE_Sh":"Gamma_OV_Sh",t,M.m,M.d,p.k,trial,native_int,native_fp,native_field,convert,pg,ser,servers,slowest,verify,decode,cleanup,convert+pg+ser+verify+decode+cleanup,request_bytes,response_bytes,vk_bytes,ntamper,rejected);fflush(stdout);}
    }
    fq_clear(constant,p.Fp);fq_mat_clear(x,p.Fp);fmpz_clear(tmp);fmpz_clear(oracle);for(int i=0;i<M.m;i++)fmpz_clear(xb[i]);free(xb);free(xd);free(wire);pv_params_clear(&p);
    fmpz_clear(M.fbias);fmpz_clear(M.fzero);free(M.a);free(M.b);free(M.coef);free(M.V);free(M.X);free(M.label);free(M.dc);free(M.dv);return 0;
}
