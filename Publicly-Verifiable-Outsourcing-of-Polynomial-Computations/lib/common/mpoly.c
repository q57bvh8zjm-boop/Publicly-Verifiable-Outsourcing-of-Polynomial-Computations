#include "mpoly.h"
#include <limits.h>

// ===========  ===========

slong mpoly_num_coefficients(int m, int d) {
    if (m < 0 || d < 0)
        return -1;

    /* C(m + d, d), computed by the exact binomial recurrence. */
    unsigned long long value = 1;
    for (int i = 1; i <= d; i++) {
        unsigned long long factor = (unsigned long long)m + (unsigned long long)i;
        if (value > ULLONG_MAX / factor)
            return -1;
        value = value * factor / (unsigned long long)i;
        if (value > (unsigned long long)LONG_MAX)
            return -1;
    }
    return (slong)value;
}

typedef struct {
    fq_struct *y;
    fq_mat_struct *F;
    fq_mat_struct *X;
    pubpar *par;
    slong next_coeff;
    slong coeff_count;
} eval_state_t;

static void eval_degree(eval_state_t *state, int start_var, int remaining,
                        const fq_t monomial) {
    if (state->next_coeff >= state->coeff_count)
        return;

    if (remaining == 0) {
        fq_t term;
        fq_init(term, state->par->Fp);
        fq_mul(term,
               fq_mat_entry(state->F, state->next_coeff, 0),
               monomial,
               state->par->Fp);
        fq_add(state->y, state->y, term, state->par->Fp);
        fq_clear(term, state->par->Fp);
        state->next_coeff++;
        return;
    }

    for (int i = start_var; i < state->par->m; i++) {
        fq_t next_monomial;
        fq_init(next_monomial, state->par->Fp);
        fq_mul(next_monomial, monomial,
               fq_mat_entry(state->X, i, 0), state->par->Fp);
        eval_degree(state, i, remaining - 1, next_monomial);
        fq_clear(next_monomial, state->par->Fp);
        if (state->next_coeff >= state->coeff_count)
            return;
    }
}

/*
 * Evaluate a dense polynomial in graded order:
 *   1; x_0,...,x_{m-1}; x_0^2,x_0x_1,...; ...
 * A short coefficient vector is accepted and interpreted as a prefix of that
 * order. Extra rows beyond total degree d are ignored.
 */
void Eval(fq_t y, fq_mat_t F, fq_mat_t X, pubpar *par) {
    fq_init(y, par->Fp);
    fq_zero(y, par->Fp);

    slong rows = fq_mat_nrows(F, par->Fp);
    if (rows == 0 || par->d < 0)
        return;

    fq_set(y, fq_mat_entry(F, 0, 0), par->Fp);

    fq_t one;
    fq_init(one, par->Fp);
    fq_one(one, par->Fp);

    eval_state_t state = { y, F, X, par, 1, rows };
    for (int degree = 1; degree <= par->d && state.next_coeff < rows; degree++)
        eval_degree(&state, 0, degree, one);

    fq_clear(one, par->Fp);
}

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
