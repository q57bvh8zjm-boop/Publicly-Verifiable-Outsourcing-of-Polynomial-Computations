/**
 * Protocol 2 - F_{q^2}
 * 
 * - KeyGen: 
 * - ProbGen:   
 * - Compute: 
 * - Verify: 
 * 
 * - F_{q^2}
 * - F_q
 * - 
 */

#include "protocol2.h"

// ===========  ===========

/**
 * 
 * 
 * @param y 
 * @param F 
 * @param X 
 * @param par 
 */
typedef struct {
    fq_struct *result;
    fq_mat_struct *F;
    fq_mat_struct *X;
    pubpar_p2 *par;
    slong next_coeff;
    slong coeff_count;
} protocol2_eval_state_t;

static void protocol2_eval_degree(protocol2_eval_state_t *state,
                                  int start_var, int remaining,
                                  const fq_t monomial) {
    if (state->next_coeff >= state->coeff_count)
        return;

    if (remaining == 0) {
        fq_t term;
        fq_init(term, state->par->Fq2);
        fq_mul(term,
               fq_mat_entry(state->F, state->next_coeff, 0),
               monomial,
               state->par->Fq2);
        fq_add(state->result, state->result, term, state->par->Fq2);
        fq_clear(term, state->par->Fq2);
        state->next_coeff++;
        return;
    }

    for (int i = start_var; i < state->par->m; i++) {
        fq_t next_monomial;
        fq_init(next_monomial, state->par->Fq2);
        fq_mul(next_monomial, monomial,
               fq_mat_entry(state->X, i, 0), state->par->Fq2);
        protocol2_eval_degree(state, i, remaining - 1, next_monomial);
        fq_clear(next_monomial, state->par->Fq2);
        if (state->next_coeff >= state->coeff_count)
            return;
    }
}

void protocol2_eval(fq_t y, fq_mat_t F, fq_mat_t X, pubpar_p2 *par) {
    fq_init(y, par->Fq2);
    fq_zero(y, par->Fq2);

    slong rows = fq_mat_nrows(F, par->Fq2);
    if (rows == 0 || par->d < 0)
        return;

    fq_set(y, fq_mat_entry(F, 0, 0), par->Fq2);
    fq_t one;
    fq_init(one, par->Fq2);
    fq_one(one, par->Fq2);

    protocol2_eval_state_t state = { y, F, X, par, 1, rows };
    for (int degree = 1; degree <= par->d && state.next_coeff < rows; degree++)
        protocol2_eval_degree(&state, 0, degree, one);

    fq_clear(one, par->Fq2);
}

/**
 * 
 * ff(SID[i]) = Y[i]
 * 
 * @param f 
 * @param Y 
 * @param par 
 */
void protocol2_interpolate(fq_poly_t f, fq_t *Y, pubpar_p2 *par) {
    fq_poly_init(f, par->Fq2);
    fq_poly_zero(f, par->Fq2);
    
    fq_poly_t g;
    fq_poly_init(g, par->Fq2);
    fq_poly_t x;
    fq_poly_init(x, par->Fq2);
    fq_poly_gen(x, par->Fq2);
    
    fq_poly_t g1, g2, h;
    fq_poly_init(g1, par->Fq2);
    fq_poly_init(g2, par->Fq2);
    fq_poly_init(h, par->Fq2);
    
    fq_t a;
    fq_init(a, par->Fq2);
    
    for (int i = 0; i < par->k; i++) {
        fq_poly_set_fq(g, Y[i], par->Fq2);
        for (int j = 0; j < par->k; j++) {
            if (j != i) {
                fq_poly_set_fq(h, par->SID[j], par->Fq2);
                fq_poly_sub(g1, x, h, par->Fq2);
                fq_sub(a, par->SID[i], par->SID[j], par->Fq2);
                fq_inv(a, a, par->Fq2);
                fq_poly_set_fq(g2, a, par->Fq2);
                fq_poly_mul(g, g, g1, par->Fq2);
                fq_poly_mul(g, g, g2, par->Fq2);
            }
        }
        fq_poly_add(f, f, g, par->Fq2);
    }
    
    // 
    fq_poly_clear(g, par->Fq2);
    fq_poly_clear(x, par->Fq2);
    fq_poly_clear(g1, par->Fq2);
    fq_poly_clear(g2, par->Fq2);
    fq_poly_clear(h, par->Fq2);
    fq_clear(a, par->Fq2);
}

// ===========  ===========

/**
 * 
 * Protocol 2
 * 
 * @param F 
 */
void protocol2_keygen(fq_mat_t F) {
    //  - 
}

// ===========  ===========

/**
 * 
 * 
 * @param alpha 
 * @param c 
 * @param X 
 * @param par 
 */
void protocol2_probgen(fq_t *alpha, fq_mat_t *c, fq_mat_t X, pubpar_p2 *par) {
    flint_rand_t state;
    flint_randinit(state);

    fq_init(*alpha, par->Fq2);
    do {
        fq_randtest(*alpha, state, par->Fq2);
    } while (fq_is_zero(*alpha, par->Fq2)); // [k]

    fq_mat_t *r = malloc(sizeof(fq_mat_t) * par->t);
    for (int i = 0; i < par->t; i++) {
        fq_mat_init(r[i], par->m, 1, par->Fq2);
        fq_mat_randtest(r[i], state, par->Fq2);
    }

    fq_t us, alphas;
    fq_init(us, par->Fq2);
    fq_init(alphas, par->Fq2);

    for (int i = 0; i < par->k; i++) {
        fq_mat_init(c[i], par->m, 1, par->Fq2);
        
        for (int j = 0; j < par->m; j++) {
            fq_set(fq_mat_entry(c[i], j, 0), fq_mat_entry(X, j, 0), par->Fq2);
        }
        
        for (int s = 0; s < par->t; s++) {
            fq_pow_ui(us, par->SID[i], s + 1, par->Fq2);
            fq_pow_ui(alphas, *alpha, s + 1, par->Fq2);
            fq_sub(us, us, alphas, par->Fq2);
            
            for (int j = 0; j < par->m; j++) {
                fq_t tmp;
                fq_init(tmp, par->Fq2);
                fq_mul(tmp, fq_mat_entry(r[s], j, 0), us, par->Fq2);
                fq_add(fq_mat_entry(c[i], j, 0), fq_mat_entry(c[i], j, 0), tmp, par->Fq2);
                fq_clear(tmp, par->Fq2);
            }
        }
    }

    // 
    fq_clear(us, par->Fq2);
    fq_clear(alphas, par->Fq2);
    for (int i = 0; i < par->t; i++) {
        fq_mat_clear(r[i], par->Fq2);
    }
    free(r);
    flint_randclear(state);
}

// ===========  ===========

/**
 * 
 * 
 * @param yi 
 * @param F 
 * @param ci 
 * @param par 
 */
void protocol2_compute(fq_t yi, fq_mat_t F, fq_mat_t ci, pubpar_p2* par) {
    protocol2_eval(yi, F, ci, par);
}

// ===========  ===========

/**
 * 
 * 
 * @param alpha 
 * @param y 
 * @param par 
 */
int protocol2_verify(fq_t alpha, fq_t *y, pubpar_p2 *par) {
    fq_poly_t phi;
    fq_poly_init(phi, par->Fq2);
    protocol2_interpolate(phi, y, par);

    fq_t val;
    fq_init(val, par->Fq2);
    fq_poly_evaluate_fq(val, phi, alpha, par->Fq2);

    // valF_q
    fq_t val_frobenius;
    fq_init(val_frobenius, par->Fq2);
    fq_frobenius(val_frobenius, val, 1, par->Fq2);
    int is_in_Fq = fq_equal(val_frobenius, val, par->Fq2);

    // 
    fq_poly_clear(phi, par->Fq2);
    fq_clear(val, par->Fq2);
    fq_clear(val_frobenius, par->Fq2);
    
    return is_in_Fq;
} 
