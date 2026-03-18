#include "scheme1.h"
#include "../lib/common/extended_field.h"




static ext_field_ctx_t g_scheme1_ext_ctx;
static ext_field_elem_t *g_scheme1_aux_R_m = NULL;
static ext_field_elem_t *g_scheme1_aux_R_n = NULL;
static ext_field_elem_t g_scheme1_alpha;
static int g_scheme1_dt = 0;
static int g_scheme1_aux_allocated_size = 0; 

static int g_scheme1_initialized = 0;



#define MAX_SCHEME1_DT 1000  

// Helper function for polynomial evaluation
static void scheme1_ext_field_poly_evaluate(ext_field_elem_t *result, ext_field_poly_t *poly, ext_field_elem_t *x, ext_field_ctx_t *ctx) {
    ext_field_poly_evaluate(result, poly, x, ctx);
}



void scheme1_ext_field_interpolate(ext_field_elem_t *phi_coeffs, ext_field_elem_t *y_values, pubpar *par) {
    // Initialize coefficients to zero
    for (int i = 0; i < par->k; i++) {
        ext_field_elem_init(&phi_coeffs[i]);
        ext_field_elem_zero(&phi_coeffs[i]);
    }
    
    // Lagrange Interpolation
    // f(u) = Σ_{j=0}^{k-1} y_j * L_j(u)
    // L_j(u) = Π_{m!=j} (u - x_m) / (x_j - x_m)
    // where x_j = j + 1
    
    ext_field_poly_t f_poly;
    ext_field_poly_init(&f_poly, par->k - 1); // Degree k-1
    // Initialize f_poly to 0
    for(int i=0; i<=f_poly.degree; i++) ext_field_elem_zero(&f_poly.coeffs[i]);

    ext_field_poly_t L_j;
    ext_field_poly_t term_poly;
    ext_field_poly_init(&term_poly, 1); // Degree 1: u - x_m
    
    ext_field_elem_t denom, temp, x_j_elem, x_m_elem, diff;
    ext_field_elem_init(&denom);
    ext_field_elem_init(&temp);
    ext_field_elem_init(&x_j_elem);
    ext_field_elem_init(&x_m_elem);
    ext_field_elem_init(&diff);
    
    fmpz_t zero_fmpz, val_fmpz;
    fmpz_init(zero_fmpz);
    fmpz_zero(zero_fmpz);
    fmpz_init(val_fmpz);
    
    for (int j = 0; j < par->k; j++) {
        // Construct L_j(u)
        ext_field_poly_init(&L_j, 0);
        ext_field_elem_one(&L_j.coeffs[0]); // L_j = 1
        
        ext_field_elem_one(&denom); // Denominator accumulator
        
        // x_j = j + 1
        fmpz_set_ui(val_fmpz, j + 1);
        ext_field_elem_set(&x_j_elem, val_fmpz, zero_fmpz, &g_scheme1_ext_ctx);
        
        for (int m = 0; m < par->k; m++) {
            if (m == j) continue;
            
            // x_m = m + 1
            fmpz_set_ui(val_fmpz, m + 1);
            ext_field_elem_set(&x_m_elem, val_fmpz, zero_fmpz, &g_scheme1_ext_ctx);
            
            // term_poly = u - x_m
            // coeff[1] = 1, coeff[0] = -x_m
            ext_field_elem_one(&term_poly.coeffs[1]);
            ext_field_elem_zero(&term_poly.coeffs[0]);
            ext_field_elem_sub(&term_poly.coeffs[0], &term_poly.coeffs[0], &x_m_elem, &g_scheme1_ext_ctx);
            
            // L_j = L_j * term_poly
            ext_field_poly_t prod;
            ext_field_poly_mul(&prod, &L_j, &term_poly, &g_scheme1_ext_ctx);
            ext_field_poly_clear(&L_j);
            L_j = prod;
            
            // denom = denom * (x_j - x_m)
            ext_field_elem_sub(&diff, &x_j_elem, &x_m_elem, &g_scheme1_ext_ctx);
            ext_field_elem_mul(&denom, &denom, &diff, &g_scheme1_ext_ctx);
        }
        
        // Scale L_j by y_j / denom
        ext_field_elem_inv(&temp, &denom, &g_scheme1_ext_ctx); // temp = 1/denom
        ext_field_elem_mul(&temp, &temp, &y_values[j], &g_scheme1_ext_ctx); // temp = y_j / denom
        
        // Add temp * L_j to f_poly
        for (int i = 0; i <= L_j.degree; i++) {
            ext_field_elem_mul(&diff, &L_j.coeffs[i], &temp, &g_scheme1_ext_ctx); // Reuse diff as temp var
            ext_field_elem_add(&f_poly.coeffs[i], &f_poly.coeffs[i], &diff, &g_scheme1_ext_ctx);
        }
        
        ext_field_poly_clear(&L_j);
    }
    
    // Copy coefficients to result
    for (int i = 0; i < par->k; i++) {
        if (i <= f_poly.degree) {
            ext_field_elem_copy(&phi_coeffs[i], &f_poly.coeffs[i]);
        } else {
            ext_field_elem_zero(&phi_coeffs[i]);
        }
    }
    
    ext_field_poly_clear(&f_poly);
    ext_field_poly_clear(&term_poly);
    ext_field_elem_clear(&denom);
    ext_field_elem_clear(&temp);
    ext_field_elem_clear(&x_j_elem);
    ext_field_elem_clear(&x_m_elem);
    ext_field_elem_clear(&diff);
    fmpz_clear(zero_fmpz);
    fmpz_clear(val_fmpz);
}




ext_field_ctx_t* scheme1_get_ext_ctx(void) {
    return &g_scheme1_ext_ctx;
}



void scheme1_init_ext_ctx(pubpar *par) {
    
    if (!g_scheme1_initialized) {
        
        fmpz_t c;
        fmpz_init(c);
        fmpz_set_str(c, "2973054622723446757380186526872783240626948608321353055178391322035449304762", 10);
        
        
        ext_field_ctx_init(&g_scheme1_ext_ctx, par->p, c);
        
        
        // 
        ext_field_elem_init(&g_scheme1_alpha);
        
        
        fmpz_clear(c);
        g_scheme1_initialized = 1;
        
    } else {
    }
}


void scheme1_init_field_extension(scheme1_field_ext_param *fep, pubpar *par) {
    scheme1_init_ext_ctx(par);
}


void scheme1_clear_field_extension(scheme1_field_ext_param *fep, pubpar *par) {
    ext_field_elem_clear(&g_scheme1_alpha);
    ext_field_ctx_clear(&g_scheme1_ext_ctx);
}


void scheme1_clear_global_aux_vars(pubpar *par) {
    if (g_scheme1_aux_R_m != NULL) {
        for (int i = 0; i < g_scheme1_aux_allocated_size; i++) {
            ext_field_elem_clear(&g_scheme1_aux_R_m[i]);
            ext_field_elem_clear(&g_scheme1_aux_R_n[i]);
        }
        free(g_scheme1_aux_R_m);
        free(g_scheme1_aux_R_n);
        g_scheme1_aux_R_m = NULL;
        g_scheme1_aux_R_n = NULL;
        g_scheme1_aux_allocated_size = 0;  
    }
}


void scheme1_eval_extended(ext_field_elem_t *result, fq_mat_t F, ext_field_mat_t *X_ext, ext_field_ctx_t *ctx, pubpar *par) {
    ext_field_elem_t temp, term;
    ext_field_elem_init(&temp);
    ext_field_elem_init(&term);

    
    fq_t c_val;
    fq_init(c_val, par->Fp);
    fq_set(c_val, fq_mat_entry(F, 0, 0), par->Fp);
    

    fmpz_t c_fmpz;
    fmpz_init(c_fmpz);
    fq2fmpz(c_fmpz, c_val, par->Fp);
    fmpz_t zero;
    fmpz_init(zero);
    fmpz_zero(zero);
    ext_field_elem_set(result, c_fmpz, zero, ctx);
    
    int index = 0;
    
    if (par->d >= 1) {
        for (int i = 0; i < par->m; i++) {
            index++;
            fq_set(c_val, fq_mat_entry(F, index, 0), par->Fp);
            fq2fmpz(c_fmpz, c_val, par->Fp);
            

            ext_field_elem_t coeff;
            ext_field_elem_init(&coeff);
            ext_field_elem_set(&coeff, c_fmpz, zero, ctx);
            
            ext_field_elem_mul(&term, &coeff, ext_field_mat_entry(X_ext, i, 0), ctx);
            ext_field_elem_add(result, result, &term, ctx);
            
            ext_field_elem_clear(&coeff);
        }
    }
    

    if (par->d >= 2) {
        for (int i = 0; i < par->m; i++) {
            for (int j = i; j < par->m; j++) {
                index++;
                fq_set(c_val, fq_mat_entry(F, index, 0), par->Fp);
                fq2fmpz(c_fmpz, c_val, par->Fp);
                
                ext_field_elem_t coeff;
                ext_field_elem_init(&coeff);
                ext_field_elem_set(&coeff, c_fmpz, zero, ctx);
                
                // term = c_{ij} * x_i * x_j
                ext_field_elem_mul(&term, ext_field_mat_entry(X_ext, i, 0), ext_field_mat_entry(X_ext, j, 0), ctx);
                ext_field_elem_mul(&term, &term, &coeff, ctx);
                ext_field_elem_add(result, result, &term, ctx);
                
                ext_field_elem_clear(&coeff);
            }
        }
    }
    

    ext_field_elem_clear(&temp);
    ext_field_elem_clear(&term);
    fq_clear(c_val, par->Fp);
    fmpz_clear(c_fmpz);
    fmpz_clear(zero);
}


void scheme1_keygen(fq_mat_t F, pubpar *par) {
    

    
    g_scheme1_dt = par->d * par->t;
    

    
    if (g_scheme1_dt > MAX_SCHEME1_DT) {
        g_scheme1_dt = 0;
        return;
    }
}



void scheme1_probgen(unsigned char *ga, ext_field_mat_t *c, ext_field_elem_t *b, fq_mat_t X, pubpar *par) {
    
    flint_rand_t state;
    flint_randinit(state);
    
    
    g_scheme1_dt = par->d * par->t;
    
    

    if (g_scheme1_dt > MAX_SCHEME1_DT) {
        flint_randclear(state);
        return;
    }
    
    // Initialize extension field context
    scheme1_init_field_extension(NULL, par);
    

    ext_field_elem_t zero_elem;
    ext_field_elem_init(&zero_elem);
    ext_field_elem_zero(&zero_elem);
    int alpha_iter = 0;
    do {
        alpha_iter++;
        ext_field_elem_random(&g_scheme1_alpha, &g_scheme1_ext_ctx);
        
        if (ext_field_elem_equal(&g_scheme1_alpha, &zero_elem)) {
            continue;
        }
        

        
        if (fmpz_is_zero(g_scheme1_alpha.b)) {
            if (fmpz_cmp_ui(g_scheme1_alpha.a, par->k) <= 0 && fmpz_cmp_ui(g_scheme1_alpha.a, 1) >= 0) {
                continue;
            }
        }
        break;
    } while (1);
    
    

    ext_field_mat_t *U = malloc(sizeof(ext_field_mat_t) * par->t);
    if (!U) {
        return;
    }
    for (int s = 0; s < par->t; s++) {
        ext_field_mat_init(&U[s], par->m, 1);
        for (int j = 0; j < par->m; j++) {
            ext_field_elem_t rand_elem;
            ext_field_elem_init(&rand_elem);
            ext_field_elem_random(&rand_elem, &g_scheme1_ext_ctx);
            ext_field_mat_entry_set(&U[s], j, 0, &rand_elem);
            ext_field_elem_clear(&rand_elem);
        }
    }
    
    // 4. 
    // Note: The original ProbGen had a different step 4 (selecting h_s)
    // and computed c(i) and a(i) later in step 7.
    // This new step 4 and 5 seem to be a simplified or alternative ProbGen logic.
    // Assuming the user intends to replace the original logic with this.
    
    // The original `c` parameter is `ext_field_mat_t *c`, which is an array of matrices.
    // The new logic uses `ext_field_mat_t *c = malloc(...)` which redefines `c`.
    // This is problematic as `c` is an output parameter.
    // I will adapt this to use the existing `c` output parameter.
    
    // Original step 4: Select h_s
    fq_t *h = malloc(sizeof(fq_t) * par->k); // index 1 to k-1
    for (int s = 1; s < par->k; s++) {
        fq_init(h[s], par->Fp);
        fq_rand(h[s], state, par->Fp);
    }
    
    // 5. α^i = m_i + n_i ω
    scheme1_clear_global_aux_vars(par);
    
    // We need powers up to degree of f(u), which is max(dt, k-1).
    // Since k >= dt+1, degree is k-1.
    // So we allocate k elements (indices 0 to k-1).
    int max_deg = par->k;
    
    g_scheme1_aux_R_m = malloc(sizeof(ext_field_elem_t) * max_deg);
    g_scheme1_aux_R_n = malloc(sizeof(ext_field_elem_t) * max_deg);
    g_scheme1_aux_allocated_size = max_deg;
    
    ext_field_elem_t alpha_pow;
    ext_field_elem_init(&alpha_pow);
    ext_field_elem_one(&alpha_pow); // α^0 = 1
    
    for (int i = 0; i < max_deg; i++) {
        ext_field_elem_init(&g_scheme1_aux_R_m[i]);
        ext_field_elem_init(&g_scheme1_aux_R_n[i]);
        
        // m_i = alpha_pow.a, n_i = alpha_pow.b
        fmpz_t zero;
        fmpz_init(zero);
        fmpz_zero(zero);
        
        ext_field_elem_set(&g_scheme1_aux_R_m[i], alpha_pow.a, zero, &g_scheme1_ext_ctx);
        ext_field_elem_set(&g_scheme1_aux_R_n[i], alpha_pow.b, zero, &g_scheme1_ext_ctx);
        
        fmpz_clear(zero);
        
        // next power
        ext_field_elem_mul(&alpha_pow, &alpha_pow, &g_scheme1_alpha, &g_scheme1_ext_ctx);
    }
    ext_field_elem_clear(&alpha_pow);
    

    
    ext_field_elem_t temp1, temp2, temp3, u_val, u_pow_s, alpha_pow_s;
    ext_field_elem_init(&temp1);
    ext_field_elem_init(&temp2);
    ext_field_elem_init(&temp3);
    ext_field_elem_init(&u_val);
    ext_field_elem_init(&u_pow_s);
    ext_field_elem_init(&alpha_pow_s);
    
    fmpz_t zero_fmpz;
    fmpz_init(zero_fmpz);
    fmpz_zero(zero_fmpz);
    
    for (int i = 0; i < par->k; i++) {
        // i+1 as element in F_{q^2}
        fmpz_t i_val;
        fmpz_init(i_val);
        fmpz_set_ui(i_val, i + 1); // SID[i] usually i+1
        ext_field_elem_set(&u_val, i_val, zero_fmpz, &g_scheme1_ext_ctx);
        fmpz_clear(i_val);
        

        // ext_field_mat_init(&c[i], par->m, 1);
        
        for (int j = 0; j < par->m; j++) {
            // Start with x_j
            fq_t x_j;
            fq_init(x_j, par->Fp);
            fq_set(x_j, fq_mat_entry(X, j, 0), par->Fp);
            fmpz_t x_j_fmpz;
            fmpz_init(x_j_fmpz);
            fq2fmpz(x_j_fmpz, x_j, par->Fp);
            
            ext_field_elem_t sum_val;
            ext_field_elem_init(&sum_val);
            ext_field_elem_set(&sum_val, x_j_fmpz, zero_fmpz, &g_scheme1_ext_ctx);
            
            for (int s = 1; s <= par->t; s++) {
                // u^s
                ext_field_elem_one(&u_pow_s);
                // Manual power for now
                for(int k=0; k<s; k++) ext_field_elem_mul(&u_pow_s, &u_pow_s, &u_val, &g_scheme1_ext_ctx);
                
                // α^s
                ext_field_elem_one(&alpha_pow_s);
                for(int k=0; k<s; k++) ext_field_elem_mul(&alpha_pow_s, &alpha_pow_s, &g_scheme1_alpha, &g_scheme1_ext_ctx);
                
                // u^s - α^s
                ext_field_elem_sub(&temp1, &u_pow_s, &alpha_pow_s, &g_scheme1_ext_ctx);
                
                // U^{(s)}_j * (u^s - α^s)
                ext_field_elem_mul(&temp2, ext_field_mat_entry(&U[s-1], j, 0), &temp1, &g_scheme1_ext_ctx);
                
                // Add to sum
                ext_field_elem_add(&sum_val, &sum_val, &temp2, &g_scheme1_ext_ctx);
            }
            
            ext_field_mat_entry_set(&c[i], j, 0, &sum_val);
            
            ext_field_elem_clear(&sum_val);
            fq_clear(x_j, par->Fp);
            fmpz_clear(x_j_fmpz);
        }
        

        // a(u) = Σ_{s=1}^{k-1} h_s (u - α)^s
        ext_field_elem_init(&b[i]);
        ext_field_elem_zero(&b[i]);
        
        // u - α
        ext_field_elem_sub(&temp1, &u_val, &g_scheme1_alpha, &g_scheme1_ext_ctx);
        
        ext_field_elem_one(&temp2); // Stores (u - α)^s
        
        for (int s = 1; s < par->k; s++) {
            // Update (u - α)^s
            ext_field_elem_mul(&temp2, &temp2, &temp1, &g_scheme1_ext_ctx);
            
            // h_s * (u - α)^s
            fmpz_t h_s_fmpz;
            fmpz_init(h_s_fmpz);
            fq2fmpz(h_s_fmpz, h[s], par->Fp);
            
            ext_field_elem_t h_elem;
            ext_field_elem_init(&h_elem);
            ext_field_elem_set(&h_elem, h_s_fmpz, zero_fmpz, &g_scheme1_ext_ctx);
            
            ext_field_elem_mul(&temp3, &h_elem, &temp2, &g_scheme1_ext_ctx);
            ext_field_elem_add(&b[i], &b[i], &temp3, &g_scheme1_ext_ctx);
            
            fmpz_clear(h_s_fmpz);
            ext_field_elem_clear(&h_elem);
        }
    }
    
    // Cleanup
    ext_field_elem_clear(&temp1);
    ext_field_elem_clear(&temp2);
    ext_field_elem_clear(&temp3);
    ext_field_elem_clear(&u_val);
    ext_field_elem_clear(&u_pow_s);
    ext_field_elem_clear(&alpha_pow_s);
    fmpz_clear(zero_fmpz);
    
    for (int s = 0; s < par->t; s++) {
        ext_field_mat_clear(&U[s]);
    }
    free(U);
    
    for (int s = 1; s < par->k; s++) {
        fq_clear(h[s], par->Fp);
    }
    free(h);
    
    flint_randclear(state);
}


void scheme1_compute(ext_field_elem_t *yi, ext_field_elem_t *zi, fq_mat_t F, ext_field_mat_t *ci, ext_field_elem_t *bi, pubpar *par) {
    // yi = F(ci)
    scheme1_eval_extended(yi, F, ci, &g_scheme1_ext_ctx, par);
    
    // yi = yi + bi (where bi = a(i))
    ext_field_elem_add(yi, yi, bi, &g_scheme1_ext_ctx);
    
    // zi = bi (optional, for consistency or debugging)
    ext_field_elem_copy(zi, bi);
}


int scheme1_verify(unsigned char *ga, ext_field_elem_t *y, ext_field_elem_t *z, pubpar *par) {
    // 1. 
    ext_field_elem_t *phi_coeffs = malloc(sizeof(ext_field_elem_t) * par->k);
    scheme1_ext_field_interpolate(phi_coeffs, y, par);
    
    
    // 2. 
    
    fmpz_t D, term, temp;
    fmpz_init(D);
    fmpz_init(term);
    fmpz_init(temp);
    fmpz_zero(D);
    
    // b_0
    fmpz_add(D, D, phi_coeffs[0].b);
    
    for (int i = 1; i < par->k; i++) {
        // a_i * n_i
        fmpz_mul(term, phi_coeffs[i].a, g_scheme1_aux_R_n[i].a);
        fmpz_add(D, D, term);
        
        // b_i * m_i
        fmpz_mul(term, phi_coeffs[i].b, g_scheme1_aux_R_m[i].a);
        fmpz_add(D, D, term);
    }
    
    fmpz_mod(D, D, g_scheme1_ext_ctx.p);
    
    int result = fmpz_is_zero(D);
    
    fmpz_clear(D);
    fmpz_clear(term);
    fmpz_clear(temp);
    
    for (int i = 0; i < par->k; i++) ext_field_elem_clear(&phi_coeffs[i]);
    free(phi_coeffs);
    
    return result;
}



void scheme1_reconstruct(fq_t result, unsigned char *ga, ext_field_elem_t *y, ext_field_elem_t *z, pubpar *par) {

    ext_field_elem_t *phi_coeffs = malloc(sizeof(ext_field_elem_t) * par->k);
    scheme1_ext_field_interpolate(phi_coeffs, y, par);
    

    
    
    fmpz_t res_val, term;
    fmpz_init(res_val);
    fmpz_init(term);
    fmpz_zero(res_val);
    
    // a_0
    fmpz_add(res_val, res_val, phi_coeffs[0].a);
    
    for (int i = 1; i < par->k; i++) {
        
        // a_i * m_i
        fmpz_mul(term, phi_coeffs[i].a, g_scheme1_aux_R_m[i].a);
        fmpz_add(res_val, res_val, term);
        
        // - c * b_i * n_i
        fmpz_mul(term, phi_coeffs[i].b, g_scheme1_aux_R_n[i].a);
        fmpz_mul(term, term, g_scheme1_ext_ctx.c);
        fmpz_sub(res_val, res_val, term);
    }
    
    fmpz_mod(res_val, res_val, g_scheme1_ext_ctx.p);
    
    fq_set_fmpz(result, res_val, par->Fp);
    
    fmpz_clear(res_val);
    fmpz_clear(term);
    
    for (int i = 0; i < par->k; i++) ext_field_elem_clear(&phi_coeffs[i]);
    free(phi_coeffs);
}
