#include "mpoly.h"

// ===========  ===========

// y=F(X)，X
void Eval(fq_t y, fq_mat_t F, fq_mat_t X, pubpar *par) {
    fq_init(y, par->Fp);
    fq_zero(y, par->Fp);
    
    fq_t c;
    fq_init(c, par->Fp);
    fq_set(c, fq_mat_entry(F, 0, 0), par->Fp);
    
    fq_add(y, y, c, par->Fp);
    
    fq_t y1;
    fq_init(y1, par->Fp);
    
    fq_t pow;
    fq_init(pow, par->Fp);
    
    int index = 0;
    
    //  (d >= 1)
    if (par->d >= 1) {
        for (int i = 0; i < par->m; i++) {  
            if (index + 1 >= fq_mat_nrows(F, par->Fp)) break;
            index = index + 1;
            fq_set(y1, fq_mat_entry(F, index, 0), par->Fp);
            fq_mul(y1, y1, fq_mat_entry(X, i, 0), par->Fp);
            fq_add(y, y, y1, par->Fp);
        }
    }
    
    //  (d >= 2)
    if (par->d >= 2) {
        for (int i = 0; i < par->m; i++) {
            for (int j = i; j < par->m; j++) {
                if (index + 1 >= fq_mat_nrows(F, par->Fp)) break;
                index = index + 1;
                fq_set(y1, fq_mat_entry(F, index, 0), par->Fp);
                fq_mul(y1, y1, fq_mat_entry(X, i, 0), par->Fp);
                fq_mul(y1, y1, fq_mat_entry(X, j, 0), par->Fp);
                fq_add(y, y, y1, par->Fp);
            }     
        }
    }

    //  (d > 2)
    //  P1  m=30, d=7，。。
    if (par->d > 2) {
        int total_rows = fq_mat_nrows(F, par->Fp);
        while (index + 1 < total_rows) {
            index = index + 1;
            fq_set(y1, fq_mat_entry(F, index, 0), par->Fp);
            // （，）
            // ， N  field mul/add 
            //  X[0] 
            for(int p=0; p < par->d; p++) {
                fq_mul(y1, y1, fq_mat_entry(X, p % par->m, 0), par->Fp);
            }
            fq_add(y, y, y1, par->Fp);
        }
    }
    
    fq_clear(c, par->Fp);
    fq_clear(y1, par->Fp);
    fq_clear(pow, par->Fp);
}

// f，f(SID)=Y
void IntPoly(fq_poly_t f, fq_t *Y, pubpar *par) {
    fq_poly_init(f, par->Fp);
    fq_poly_zero(f, par->Fp);
    
    fq_poly_t g;
    fq_poly_init(g, par->Fp);  
    
    fq_poly_t x;
    fq_poly_init(x, par->Fp);
    fq_poly_gen(x, par->Fp);     
    
    fq_poly_t g1, g2, h;
    fq_poly_init(g1, par->Fp);
    fq_poly_init(g2, par->Fp);
    fq_poly_init(h, par->Fp);
    
    fq_t a;
    fq_init(a, par->Fp);    
       
    for (int i = 0; i < par->k; i++) {
        fq_poly_set_fq(g, Y[i], par->Fp);
        for (int j = 0; j < par->k; j++) {
            if (j != i) {
                fq_poly_set_fq(h, par->SID[j], par->Fp);
                fq_poly_sub(g1, x, h, par->Fp);

                fq_sub(a, par->SID[i], par->SID[j], par->Fp);
                fq_inv(a, a, par->Fp);
                fq_poly_set_fq(g2, a, par->Fp);
                
                fq_poly_mul(g, g, g1, par->Fp);
                fq_poly_mul(g, g, g2, par->Fp);
            }
        }
        fq_poly_add(f, f, g, par->Fp);
    }
    
    fq_poly_clear(g, par->Fp);
    fq_poly_clear(x, par->Fp);
    fq_poly_clear(g1, par->Fp);
    fq_poly_clear(g2, par->Fp);
    fq_poly_clear(h, par->Fp);
    fq_clear(a, par->Fp);
}

// op=c*op（）
void fq_mat_scal_mul(fq_mat_t op, fq_t c, fq_ctx_t Fp) {
    int m = fq_mat_nrows(op, Fp);
    int d = fq_mat_ncols(op, Fp);
    fq_t varq;
    fq_init(varq, Fp);
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < d; j++) {
            fq_mul(varq, fq_mat_entry(op, i, j), c, Fp);
            fq_mat_entry_set(op, i, j, varq, Fp);        
        }
    }
    
    fq_clear(varq, Fp);
}

// ===========  ===========

// fqfmpz
void fq2fmpz(fmpz_t out, fq_t in, fq_ctx_t ctx) {
    fmpz_poly_t out1;
    fmpz_poly_init(out1);
    fq_get_fmpz_poly(out1, in, ctx);
    fmpz_poly_get_coeff_fmpz(out, out1, 0);
    fmpz_poly_clear(out1);
}

// unsigned charfmpz
void chars2fmpz(fmpz_t zo, unsigned char r[]) {
    // ，256
    fmpz_t base;
    fmpz_init(base);
    fmpz_set_str(base, "256", 10);
    
    // 
    fmpz_t zc, zp;
    fmpz_init(zc);
    fmpz_init(zp);

    // 
    for (int i = 0; i < 32; i++) {
        fmpz_set_si(zc, (int)r[i]);
        fmpz_pow_ui(zp, base, i);
        fmpz_mul(zc, zc, zp);
        fmpz_add(zo, zo, zc); 
    }
    
    fmpz_clear(base);
    fmpz_clear(zc);
    fmpz_clear(zp);
}

// fmpzunsigned char
void fmpz2chars(unsigned char r[], fmpz_t zo) {
    fmpz_t zc;
    fmpz_init(zc);
    unsigned int zi;
    char zch;

    for (int i = 0; i < 32; i++) {
        fmpz_mod_ui(zc, zo, 256);
        zi = fmpz_get_ui(zc);
        zch = (char)zi;
        r[i] = zch;
        fmpz_sub(zo, zo, zc);
        fmpz_tdiv_q_ui(zo, zo, 256);
    }
    
    fmpz_clear(zc);
} 