/**
 * Protocol 2 - F_{q^2}
 * 
 * F_{q^2}，：
 * - KeyGen: 
 * - ProbGen:   
 * - Compute: 
 * - Verify: 
 * 
 * ：
 * - F_{q^2}
 * - F_q
 * - 
 */

#include "protocol2.h"

// ===========  ===========

/**
 * 
 * y = F(X)，F
 * 
 * @param y 
 * @param F 
 * @param X 
 * @param par 
 */
void protocol2_eval(fq_t y, fq_mat_t F, fq_mat_t X, pubpar_p2 *par) {
    fq_init(y, par->Fq2);
    fq_zero(y, par->Fq2);
    
    fq_t c;
    fq_init(c, par->Fq2);
    
    // 
    fq_set(c, fq_mat_entry(F, 0, 0), par->Fq2);
    fq_add(y, y, c, par->Fq2);
    
    fq_t y1;
    fq_init(y1, par->Fq2);
    
    if (par->d == 2) {
        int index = 0;
        
        // ：Σᵢ aᵢxᵢ
        for (int i = 0; i < par->m; i++) {
            index = index + 1;
            fq_set(y1, fq_mat_entry(F, index, 0), par->Fq2);
            fq_mul(y1, y1, fq_mat_entry(X, i, 0), par->Fq2);
            fq_add(y, y, y1, par->Fq2);
        }
        
        // ：Σᵢ≤ⱼ aᵢⱼxᵢxⱼ
        for (int i = 0; i < par->m; i++)
        for (int j = i; j < par->m; j++) {
            index = index + 1;
            fq_set(y1, fq_mat_entry(F, index, 0), par->Fq2);
            fq_mul(y1, y1, fq_mat_entry(X, i, 0), par->Fq2);
            fq_mul(y1, y1, fq_mat_entry(X, j, 0), par->Fq2);
            fq_add(y, y, y1, par->Fq2);
        }
    }
    
    fq_clear(c, par->Fq2);
    fq_clear(y1, par->Fq2);
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
    
    // ：f(x) = Σᵢ Yᵢ · Πⱼ≠ᵢ (x - SIDⱼ)/(SIDᵢ - SIDⱼ)
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
 * F_{q^2}α，
 * 
 * @param alpha 
 * @param c 
 * @param X 
 * @param par 
 */
void protocol2_probgen(fq_t *alpha, fq_mat_t *c, fq_mat_t X, pubpar_p2 *par) {
    flint_rand_t state;
    flint_randinit(state);

    //  α ∈ F_{q^2}^* \ [k]
    fq_init(*alpha, par->Fq2);
    do {
        fq_randtest(*alpha, state, par->Fq2);
    } while (fq_is_zero(*alpha, par->Fq2)); // [k]

    // tr₁,...,rₜ
    fq_mat_t *r = malloc(sizeof(fq_mat_t) * par->t);
    for (int i = 0; i < par->t; i++) {
        fq_mat_init(r[i], par->m, 1, par->Fq2);
        fq_mat_randtest(r[i], state, par->Fq2);
    }

    // c(u) = x + Σₛ₌₁ᵗ rₛ(uˢ - αˢ)
    fq_t us, alphas;
    fq_init(us, par->Fq2);
    fq_init(alphas, par->Fq2);

    for (int i = 0; i < par->k; i++) {
        fq_mat_init(c[i], par->m, 1, par->Fq2);
        
        // c[i] = x （）
        for (int j = 0; j < par->m; j++) {
            fq_set(fq_mat_entry(c[i], j, 0), fq_mat_entry(X, j, 0), par->Fq2);
        }
        
        // ：Σₛ₌₁ᵗ rₛ(SID[i]ˢ - αˢ)
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
 * F(cᵢ)
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
 * φ(α)F_q
 * 
 * @param alpha 
 * @param y 
 * @param par 
 * @return （1，0）
 */
int protocol2_verify(fq_t alpha, fq_t *y, pubpar_p2 *par) {
    fq_poly_t phi;
    fq_poly_init(phi, par->Fq2);
    protocol2_interpolate(phi, y, par);

    fq_t val;
    fq_init(val, par->Fq2);
    fq_poly_evaluate_fq(val, phi, alpha, par->Fq2);

    // valF_q
    // ：Frobenius(val) == val
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