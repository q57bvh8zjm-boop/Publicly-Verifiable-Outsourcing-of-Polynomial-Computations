#include "scheme2.h"

// ========== Protocol 4  ==========

/**
 * Protocol 4 KeyGen - MSVC
 * scheme2.md Figure KeyGen2KeyGen
 */
void scheme2_keygen(fq_mat_t F, pubpar *par) {
    
    
    
    // - 
    // - F
    // - 
}

/**
 * Protocol 4 ProbGen
 * scheme2.md Figure ProbGen2ProbGen
 */
void scheme2_probgen(scheme2_probgen_out *out, fq_mat_t x, pubpar *par) {
    int d = par->d;
    int t = par->t;
    int k = par->k;
    int m = par->m;
    
    out->vk_x = (unsigned char**)malloc((d+1) * sizeof(unsigned char*));
    if (!out->vk_x) {
        // 
        return;
    }
    
    for(int j = 0; j <= d; j++) {
        out->vk_x[j] = (unsigned char*)malloc(32 * sizeof(unsigned char)); // libsodium group element size
        if (!out->vk_x[j]) {
            // 
            for(int i = 0; i < j; i++) {
                free(out->vk_x[i]);
            }
            free(out->vk_x);
            return;
        }
    }
    
    out->sigma = (fq_mat_t**)malloc(k * sizeof(fq_mat_t*));
    if (!out->sigma) {
        //  vk_x
        for(int j = 0; j <= d; j++) {
            free(out->vk_x[j]);
        }
        free(out->vk_x);
        return;
    }
    
    for(int i = 0; i < k; i++) {
        out->sigma[i] = (fq_mat_t*)malloc((d+1) * sizeof(fq_mat_t));
        if (!out->sigma[i]) {
            // 
            for(int j = 0; j < i; j++) {
                for(int l = 0; l <= d; l++) {
                    fq_mat_clear(out->sigma[j][l], par->Fp);
                }
                free(out->sigma[j]);
            }
            free(out->sigma);
            for(int j = 0; j <= d; j++) {
                free(out->vk_x[j]);
            }
            free(out->vk_x);
            return;
        }
        
        for(int j = 0; j <= d; j++) {
            fq_mat_init(out->sigma[i][j], m, 1, par->Fp);
        }
    }
    
    out->v_values = (fq_t*)malloc(d * sizeof(fq_t));
    if (!out->v_values) {
        // 
        for(int i = 0; i < k; i++) {
            for(int j = 0; j <= d; j++) {
                fq_mat_clear(out->sigma[i][j], par->Fp);
            }
            free(out->sigma[i]);
        }
        free(out->sigma);
        for(int j = 0; j <= d; j++) {
            free(out->vk_x[j]);
        }
        free(out->vk_x);
        return;
    }
    
    for(int j = 0; j < d; j++) {
        fq_init(out->v_values[j], par->Fp);
    }
    
    for(int j = 0; j < d; j++) {
        uint32_t rand_val = randombytes_random();
        fq_set_ui(out->v_values[j], (rand_val % 1000) + 2, par->Fp); // 01
    }
    
    fq_t one;
    fq_init(one, par->Fp);
    fq_one(one, par->Fp);
    
    for(int i = 0; i < d; i++) {
        // 1
        if(fq_equal(out->v_values[i], one, par->Fp)) {
            uint32_t rand_val = randombytes_random();
            fq_set_ui(out->v_values[i], (rand_val % 1000) + 100, par->Fp);
        }
        
        for(int j = i + 1; j < d; j++) {
            if(fq_equal(out->v_values[i], out->v_values[j], par->Fp)) {
                // v_j
                uint32_t rand_val = randombytes_random();
                fq_set_ui(out->v_values[j], (rand_val % 1000) + 100, par->Fp);
                j = i; // 
            }
        }
    }
    
    fq_clear(one, par->Fp);
    
    fq_mat_t **U = (fq_mat_t**)malloc((d+1) * sizeof(fq_mat_t*));
    for(int j = 0; j <= d; j++) {
        U[j] = (fq_mat_t*)malloc(t * sizeof(fq_mat_t));
        for(int i = 0; i < t; i++) {
            fq_mat_init(U[j][i], m, 1, par->Fp);
            fq_t temp_val;
            fq_init(temp_val, par->Fp);
            for(int row = 0; row < m; row++) {
                uint32_t rand_val = randombytes_random();
                fq_set_ui(temp_val, (rand_val % 500) + 1, par->Fp);
                fq_mat_entry_set(U[j][i], row, 0, temp_val, par->Fp);
            }
            fq_clear(temp_val, par->Fp);
        }
    }
    
    // 2. d+1:
    
    fq_t ui_power, temp;
    fq_init(ui_power, par->Fp);
    fq_init(temp, par->Fp);
    
    for(int i = 1; i <= k; i++) {
        fq_mat_set(out->sigma[i-1][0], x, par->Fp);
        
        for(int poly_deg = 1; poly_deg <= t; poly_deg++) {
            // i^poly_deg
            fq_set_si(ui_power, i, par->Fp);
            for(int exp = 1; exp < poly_deg; exp++) {
                fq_set_si(temp, i, par->Fp);
                fq_mul(ui_power, ui_power, temp, par->Fp);
            }
            
            fq_mat_t temp_mat;
            fq_mat_init(temp_mat, m, 1, par->Fp);
            fq_mat_set(temp_mat, U[0][poly_deg-1], par->Fp);
            fq_mat_scal_mul(temp_mat, ui_power, par->Fp);
            fq_mat_add(out->sigma[i-1][0], out->sigma[i-1][0], temp_mat, par->Fp);
            fq_mat_clear(temp_mat, par->Fp);
        }
        
        for(int j = 1; j <= d; j++) {
            // v_j * x
            fq_mat_set(out->sigma[i-1][j], x, par->Fp);
            fq_mat_scal_mul(out->sigma[i-1][j], out->v_values[j-1], par->Fp);
            
            for(int poly_deg = 1; poly_deg <= t; poly_deg++) {
                // i^poly_deg
                fq_set_si(ui_power, i, par->Fp);
                for(int exp = 1; exp < poly_deg; exp++) {
                    fq_set_si(temp, i, par->Fp);
                    fq_mul(ui_power, ui_power, temp, par->Fp);
                }
                
                fq_mat_t temp_mat;
                fq_mat_init(temp_mat, m, 1, par->Fp);
                fq_mat_set(temp_mat, U[j][poly_deg-1], par->Fp);
                fq_mat_scal_mul(temp_mat, ui_power, par->Fp);
                fq_mat_add(out->sigma[i-1][j], out->sigma[i-1][j], temp_mat, par->Fp);
                fq_mat_clear(temp_mat, par->Fp);
            }
        }
    }
    
    fq_t *lambda_prime = (fq_t*)malloc((d+1) * sizeof(fq_t));
    for(int j = 0; j <= d; j++) {
        fq_init(lambda_prime[j], par->Fp);
    }
    scheme2_compute_lambda_prime(lambda_prime, out->v_values, d, par);
    
    for(int j = 0; j <= d; j++) {
        fmpz_t lambda_prime_fmpz;
        fmpz_init(lambda_prime_fmpz);
        fq2fmpz(lambda_prime_fmpz, lambda_prime[j], par->Fp);
        
        unsigned char lambda_prime_chars[32];
        fmpz2chars(lambda_prime_chars, lambda_prime_fmpz);
        
        int ms = crypto_scalarmult_ristretto255_base(out->vk_x[j], lambda_prime_chars);
        
        if (ms != 0) {
            memset(out->vk_x[j], 0, 32);
        }
        
        fmpz_clear(lambda_prime_fmpz);
    }
    
    // 
    for(int j = 0; j <= d; j++) {
        for(int i = 0; i < t; i++) {
            fq_mat_clear(U[j][i], par->Fp);
        }
        free(U[j]);
        fq_clear(lambda_prime[j], par->Fp);
    }
    free(U);
    free(lambda_prime);
    fq_clear(ui_power, par->Fp);
    fq_clear(temp, par->Fp);
}

/**
 * Protocol 4 Compute
 * scheme2.md Figure Compute2Compute
 */
void scheme2_compute(scheme2_compute_out *out, int server_id, fq_mat_t F, 
                       fq_mat_t *sigma_i, pubpar *par) {
    int d = par->d;
    
    // 
    out->pi = (fq_mat_t*)malloc((d+1) * sizeof(fq_mat_t));
    for(int j = 0; j <= d; j++) {
        fq_mat_init(out->pi[j], 1, 1, par->Fp); // 1x1
    }
    
    // scheme2.md Figure Compute2:
    
    for(int j = 0; j <= d; j++) {
        fq_t result;
        fq_init(result, par->Fp);
        
        Eval(result, F, sigma_i[j], par);
        
        // 
        fq_mat_entry_set(out->pi[j], 0, 0, result, par->Fp);
        
        fq_clear(result, par->Fp);
    }
    
}

/**
 * Protocol 4 Verify
 * scheme2.md Figure Verify2Verify
 */
int scheme2_verify(scheme2_verify_out *out, unsigned char **vk_x, 
                     scheme2_compute_out *pi_results,                      fq_t *v_values, 
                     fq_mat_t F, pubpar *par) {
    int d = par->d;
    int k = par->k;
    
    // 
    if (!out || !vk_x || !pi_results || !v_values || !par) {
        if (out) {
            out->valid = 0;
        }
        return 0;
    }
    
    if (d <= 0 || k <= 0) {
        out->valid = 0;
        return 0;
    }
    
    // 
    fq_mat_init(out->aux_R_prime, 1, 1, par->Fp);
    out->valid = 0;
    
    // scheme2.md Figure Verify2:
    // 3. f(0), f_{v_1}(0), ..., f_{v_d}(0)
    
    fq_poly_t *f_polys = (fq_poly_t*)malloc((d+1) * sizeof(fq_poly_t));
    if (!f_polys) {
        out->valid = 0;
        return out->valid;
    }
    
    for(int j = 0; j <= d; j++) {
        fq_poly_init(f_polys[j], par->Fp);
    }
    
    for(int j = 0; j <= d; j++) {
        fq_t *y_values = (fq_t*)malloc(k * sizeof(fq_t));
        for(int i = 0; i < k; i++) {
            fq_init(y_values[i], par->Fp);
            fq_set(y_values[i], fq_mat_entry(pi_results[i].pi[j], 0, 0), par->Fp);
        }
        
        // pi_resultsSID
        IntPoly(f_polys[j], y_values, par);
        
        // 
        for(int i = 0; i < k; i++) {
            fq_clear(y_values[i], par->Fp);
        }
        free(y_values);
    }
    
    // f(0), f_{v_1}(0), ..., f_{v_d}(0)
    fq_t *f_at_zero = (fq_t*)malloc((d+1) * sizeof(fq_t));
    fq_t zero;
    fq_init(zero, par->Fp);
    fq_zero(zero, par->Fp);
    
    for(int j = 0; j <= d; j++) {
        fq_init(f_at_zero[j], par->Fp);
        fq_poly_evaluate_fq(f_at_zero[j], f_polys[j], zero, par->Fp); // 0
    }
    
    
    // 1. FF_0
    fq_t F_0;
    fq_init(F_0, par->Fp);
    fq_set(F_0, fq_mat_entry(F, 0, 0), par->Fp);
    
    fq_t *exponents = (fq_t*)malloc((d+1) * sizeof(fq_t));
    for(int i = 0; i <= d; i++) {
        fq_init(exponents[i], par->Fp);
        fq_sub(exponents[i], f_at_zero[i], F_0, par->Fp);
    }
    
    // libsodium
    unsigned char result[32];
    unsigned char temp[32];
    
    memset(result, 0, 32);
    
    for(int i = 0; i <= d; i++) {
        // 
        fmpz_t exp_fmpz;
        fmpz_init(exp_fmpz);
        fq2fmpz(exp_fmpz, exponents[i], par->Fp);
        
        unsigned char exp_chars[32];
        fmpz2chars(exp_chars, exp_fmpz);
        
        int ms = crypto_scalarmult_ristretto255(temp, exp_chars, vk_x[i]);
        
        if (ms == 0) {
            crypto_core_ristretto255_add(result, result, temp);
        }
        
        fmpz_clear(exp_fmpz);
    }
    
    // 4. 
    unsigned char identity[32];
    memset(identity, 0, 32);
    
    if (memcmp(result, identity, 32) == 0) {
        out->valid = 1;
    } else {
        out->valid = 0;
    }
    
    // 
    fq_clear(F_0, par->Fp);
    for(int i = 0; i <= d; i++) {
        fq_clear(exponents[i], par->Fp);
    }
    free(exponents);
    
    if(out->valid) {
        // f(0)
        fq_mat_entry_set(out->aux_R_prime, 0, 0, f_at_zero[0], par->Fp);
    }
    
    // 
    for(int j = 0; j <= d; j++) {
        fq_poly_clear(f_polys[j], par->Fp);
        fq_clear(f_at_zero[j], par->Fp);
    }
    free(f_polys);
    free(f_at_zero);
    fq_clear(zero, par->Fp);
    
    return out->valid;
}

/**
 * Protocol 4 Reconstruct
 * scheme2.md Figure Reconstruct2Reconstruct
 * : F(x)
 */
void scheme2_reconstruct(scheme2_reconstruct_out *out, fq_mat_t aux_R_prime, pubpar *par) {
    // 1x1
    fq_mat_init(out->result, 1, 1, par->Fp);
    
    // scheme2.md Figure Reconstruct2:
    // 2. F(x) = aux'_R
    
    
    // aux'_Rresult
    fq_mat_set(out->result, aux_R_prime, par->Fp);
    
}

// ========== Protocol 4  ==========

/**
 * tilde_lambda_j = (Product_{i=0,i!=j}^d v_i) * (Product_{i,k in {0,1,...,d}, i<k} (v_k - v_i)) / (Product_{i=0,i!=j}^d (v_j - v_i))
 */
void scheme2_compute_lambda_prime(fq_t *lambda_prime, fq_t *v_values, int d, pubpar *par) {
    fq_t *v = (fq_t*)malloc((d+1) * sizeof(fq_t));
    for(int i = 0; i <= d; i++) {
        fq_init(v[i], par->Fp);
        if (i == 0) {
            fq_one(v[i], par->Fp); // v_0 = 1
        } else {
            fq_set(v[i], v_values[i-1], par->Fp); // v_1, v_2, ..., v_d
        }
    }
    
    for(int j = 0; j <= d; j++) {
        fq_t numerator, denominator, temp;
        fq_init(numerator, par->Fp);
        fq_init(denominator, par->Fp);
        fq_init(temp, par->Fp);
        
        // 1
        fq_one(numerator, par->Fp);
        fq_one(denominator, par->Fp);
        
        for(int i = 0; i <= d; i++) {
            if (i != j) {
                fq_mul(numerator, numerator, v[i], par->Fp);
            }
        }
        
        fq_t second_part;
        fq_init(second_part, par->Fp);
        fq_one(second_part, par->Fp);
        
        for(int i = 0; i <= d; i++) {
            for(int k = i + 1; k <= d; k++) {
                fq_sub(temp, v[k], v[i], par->Fp);
                fq_mul(second_part, second_part, temp, par->Fp);
            }
        }
        
        // 
        fq_mul(numerator, numerator, second_part, par->Fp);
        
        for(int i = 0; i <= d; i++) {
            if (i != j) {
                fq_sub(temp, v[j], v[i], par->Fp);
                fq_mul(denominator, denominator, temp, par->Fp);
            }
        }
        
        if (!fq_is_zero(denominator, par->Fp)) {
            fq_inv(temp, denominator, par->Fp);
            fq_mul(lambda_prime[j], numerator, temp, par->Fp);
        } else {
            fq_one(lambda_prime[j], par->Fp);
        }
        
        // 
        fq_clear(second_part, par->Fp);
        fq_clear(numerator, par->Fp);
        fq_clear(denominator, par->Fp);
        fq_clear(temp, par->Fp);
    }
    
    // v
    for(int i = 0; i <= d; i++) {
        fq_clear(v[i], par->Fp);
    }
    free(v);
}

/**
 * 
 */
void scheme2_interpolate_polynomials(fq_poly_t *f_polys, scheme2_compute_out *pi_results, 
                                       fq_t *v_values, pubpar *par) {
    int d = par->d;
    int k = par->k;
    
    for(int j = 0; j <= d; j++) {
        fq_t *y_values = (fq_t*)malloc(k * sizeof(fq_t));
        for(int i = 0; i < k; i++) {
            fq_init(y_values[i], par->Fp);
            fq_set(y_values[i], fq_mat_entry(pi_results[i].pi[j], 0, 0), par->Fp);
        }
        
        // 
        IntPoly(f_polys[j], y_values, par);
        
        // 
        for(int i = 0; i < k; i++) {
            fq_clear(y_values[i], par->Fp);
        }
        free(y_values);
    }
}

// ========== Protocol 4  ==========

/**
 * probgen
 */
void scheme2_probgen_out_clear(scheme2_probgen_out *out, pubpar *par) {
    int d = par->d;
    int k = par->k;
    
    // vk_x
    if(out->vk_x) {
        for(int j = 0; j <= d; j++) {
            if(out->vk_x[j]) {
                free(out->vk_x[j]);
            }
        }
        free(out->vk_x);
    }
    
    // sigma
    if(out->sigma) {
        for(int i = 0; i < k; i++) {
            if(out->sigma[i]) {
                for(int j = 0; j <= d; j++) {
                    fq_mat_clear(out->sigma[i][j], par->Fp);
                }
                free(out->sigma[i]);
            }
        }
        free(out->sigma);
    }
    
    // v_values
    if(out->v_values) {
        for(int j = 0; j < d; j++) {
            fq_clear(out->v_values[j], par->Fp);
        }
        free(out->v_values);
    }
}

/**
 * compute
 */
void scheme2_compute_out_clear(scheme2_compute_out *out, pubpar *par) {
    int d = par->d;
    
    if(out->pi) {
        for(int j = 0; j <= d; j++) {
            fq_mat_clear(out->pi[j], par->Fp);
        }
        free(out->pi);
    }
}

/**
 * verify
 */
void scheme2_verify_out_clear(scheme2_verify_out *out, pubpar *par) {
    fq_mat_clear(out->aux_R_prime, par->Fp);
}

/**
 * reconstruct
 */
void scheme2_reconstruct_out_clear(scheme2_reconstruct_out *out, pubpar *par) {
    fq_mat_clear(out->result, par->Fp);
} 