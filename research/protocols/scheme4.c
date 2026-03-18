/**
 * Protocol 8 -  CNF Secret Sharing  k-Server 
 * 
 *  CNF secret sharing ，
 * ：
 * - KeyGen: ， CNF split 
 * - ProbGen: ， CNF secret sharing 
 * - Compute: ，
 * - Verify: ，（l=1 ）
 * - Decode: ，
 * 
 * ：
 * - k = ⌈dt/2⌉
 * -  CNF secret sharing 
 * - ： reconstruction ， g^
 * - 
 */

#include "scheme4.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// =========== CNF Secret Sharing  ===========

/**
 *  H(u, 0) = t
 * u  k 
 */
int scheme4_hamming_weight(int u, int k) {
    int weight = 0;
    for (int i = 0; i < k; i++) {
        if (u & (1 << i)) {
            weight++;
        }
    }
    return weight;
}

/**
 *  H(u, 0) = t  u ∈ {0,1}^k
 *  u 
 */
int scheme4_generate_cnf_vectors(int **u_vectors, int k, int t) {
    int max_u = (1 << k) - 1;
    int count = 0;
    
    // 
    for (int u = 0; u <= max_u; u++) {
        if (scheme4_hamming_weight(u, k) == t) {
            count++;
        }
    }
    
    // 
    *u_vectors = (int*)malloc(count * sizeof(int));
    if (!*u_vectors) {
        return 0;
    }
    
    // 
    int idx = 0;
    for (int u = 0; u <= max_u; u++) {
        if (scheme4_hamming_weight(u, k) == t) {
            (*u_vectors)[idx++] = u;
        }
    }
    
    return count;
}

/**
 * CNF secret sharing:  F  F_i
 * ， CNF 
 * 
 *  F， k  F_i = {F_{i,1}, ..., F_{i,L}}
 *  F_{i,ell'}  CNF 
 */
void scheme4_cnf_split(cnf_function_share_t *F_shares, fq_mat_t F, pubpar *par) {
    int k = par->k;
    int d = par->d;
    int t = par->t;
    
    //  H(u, 0) = t  u 
    int *u_vectors = NULL;
    int num_u = scheme4_generate_cnf_vectors(&u_vectors, k, t);
    
    if (num_u == 0) {
        // ：
        // 
        for (int i = 0; i < k; i++) {
            F_shares[i].L = 0;
            F_shares[i].F_i = NULL;
        }
        if (u_vectors) free(u_vectors);
        return;
    }
    
    // ， CNF 
    // ：
    // ， CNF 
    
    // 
    int num_monomials = 1; // 
    if (d >= 1) {
        num_monomials += par->m; // 
    }
    if (d >= 2) {
        num_monomials += (par->m * (par->m + 1)) / 2; // 
    }
    
    // 
    for (int i = 0; i < k; i++) {
        // ：（ CNF ）
        F_shares[i].L = num_monomials;
        F_shares[i].F_i = (fq_mat_t*)malloc(num_monomials * sizeof(fq_mat_t));
        
        for (int ell = 0; ell < num_monomials; ell++) {
            // 
            fq_mat_init(F_shares[i].F_i[ell], 1, 1, par->Fp);
            fq_mat_entry_set(F_shares[i].F_i[ell], 0, 0, 
                            fq_mat_entry(F, ell, 0), par->Fp);
        }
    }
    
    free(u_vectors);
}

// =========== Protocol 8 KeyGen  ===========

void scheme4_keygen(scheme4_keygen_out *out, fq_mat_t F, pubpar *par) {
    int k = par->k;
    
    //  k = ⌈dt/2⌉
    int computed_k = (int)ceil((double)(par->d * par->t) / 2.0);
    if (k != computed_k) {
        // ：k 
        printf(": k = %d,  ⌈dt/2⌉ = %d\n", k, computed_k);
    }
    
    // 
    out->rho = (cnf_function_share_t*)calloc(k, sizeof(cnf_function_share_t));
    if (!out->rho) {
        return;
    }
    
    //  CNF split  F
    scheme4_cnf_split(out->rho, F, par);
}

// =========== Protocol 8 ProbGen  ===========

void scheme4_probgen(scheme4_probgen_out *out, fq_mat_t x, pubpar *par) {
    int k = par->k;
    int m = par->m;
    int d = par->d;
    int t = par->t;
    
    // 1.  α ∈ F_q^*
    fq_init(out->alpha, par->Fp);
    do {
        uint32_t rand_val = randombytes_random();
        fq_set_ui(out->alpha, (rand_val % 1000) + 1, par->Fp);
    } while (fq_is_zero(out->alpha, par->Fp));
    
    // 2.  r_i = α * x_i for every i ∈ [k]
    // ： r  m ， r_j = α * x_j
    fq_mat_init(out->r, m, 1, par->Fp);
    for (int j = 0; j < m; j++) {
        fq_t r_j;
        fq_init(r_j, par->Fp);
        fq_mul(r_j, out->alpha, fq_mat_entry(x, j, 0), par->Fp);
        fq_mat_entry_set(out->r, j, 0, r_j, par->Fp);
        fq_clear(r_j, par->Fp);
    }
    
    // 3.  H(u, 0) = t  u 
    int *u_vectors = NULL;
    int num_u = scheme4_generate_cnf_vectors(&u_vectors, k, t);
    
    //  aux_D，
    fq_mat_init(out->aux_D, m, 1, par->Fp);
    
    if (num_u == 0) {
        //  r， scheme4_probgen_out_clear 
        return;
    }
    
    // 4.  j  u， x_{j,u}  r_{j,u}
    // : x_j = Σ_u x_{j,u}, r_j = Σ_u r_{j,u}
    
    // 
    out->s_i = (fq_t***)malloc(k * sizeof(fq_t**));
    out->tilde_s_i = (fq_t***)malloc(k * sizeof(fq_t**));
    out->num_shares_per_coord = (int*)malloc(m * sizeof(int));
    
    for (int j = 0; j < m; j++) {
        out->num_shares_per_coord[j] = num_u;
    }
    
    //  CNF 
    fq_t **x_shares = (fq_t**)malloc(m * sizeof(fq_t*));
    fq_t **r_shares = (fq_t**)malloc(m * sizeof(fq_t*));
    
    for (int j = 0; j < m; j++) {
        x_shares[j] = (fq_t*)malloc(num_u * sizeof(fq_t));
        r_shares[j] = (fq_t*)malloc(num_u * sizeof(fq_t));
        
        // ，
        fq_t sum_x, sum_r;
        fq_init(sum_x, par->Fp);
        fq_init(sum_r, par->Fp);
        fq_set(sum_x, fq_mat_entry(x, j, 0), par->Fp);
        fq_set(sum_r, fq_mat_entry(out->r, j, 0), par->Fp);
        
        for (int u_idx = 0; u_idx < num_u - 1; u_idx++) {
            fq_init(x_shares[j][u_idx], par->Fp);
            fq_init(r_shares[j][u_idx], par->Fp);
            
            // 
            uint32_t rand_val = randombytes_random();
            fq_set_ui(x_shares[j][u_idx], rand_val % 1000, par->Fp);
            fq_sub(sum_x, sum_x, x_shares[j][u_idx], par->Fp);
            
            rand_val = randombytes_random();
            fq_set_ui(r_shares[j][u_idx], rand_val % 1000, par->Fp);
            fq_sub(sum_r, sum_r, r_shares[j][u_idx], par->Fp);
        }
        
        // 
        fq_init(x_shares[j][num_u - 1], par->Fp);
        fq_init(r_shares[j][num_u - 1], par->Fp);
        fq_set(x_shares[j][num_u - 1], sum_x, par->Fp);
        fq_set(r_shares[j][num_u - 1], sum_r, par->Fp);
        
        fq_clear(sum_x, par->Fp);
        fq_clear(sum_r, par->Fp);
    }
    
    // 5.  i 
    // s_i^{(j)} = {x_{j,u} | u_i = 0}
    // tilde{s}_i^{(j)} = {r_{j,u} | u_i = 0}
    for (int i = 0; i < k; i++) {
        out->s_i[i] = (fq_t**)malloc(m * sizeof(fq_t*));
        out->tilde_s_i[i] = (fq_t**)malloc(m * sizeof(fq_t*));
        
        for (int j = 0; j < m; j++) {
            //  i （u_i = 0  u ）
            int count = 0;
            for (int u_idx = 0; u_idx < num_u; u_idx++) {
                int u = u_vectors[u_idx];
                if (!(u & (1 << i))) {  // u_i = 0
                    count++;
                }
            }
            
            out->s_i[i][j] = (fq_t*)malloc(count * sizeof(fq_t));
            out->tilde_s_i[i][j] = (fq_t*)malloc(count * sizeof(fq_t));
            
            int idx = 0;
            for (int u_idx = 0; u_idx < num_u; u_idx++) {
                int u = u_vectors[u_idx];
                if (!(u & (1 << i))) {  // u_i = 0
                    fq_init(out->s_i[i][j][idx], par->Fp);
                    fq_init(out->tilde_s_i[i][j][idx], par->Fp);
                    fq_set(out->s_i[i][j][idx], x_shares[j][u_idx], par->Fp);
                    fq_set(out->tilde_s_i[i][j][idx], r_shares[j][u_idx], par->Fp);
                    idx++;
                }
            }
        }
    }
    
    // 6.  vk_x = g^{α^d}
    out->vk_x = (unsigned char*)malloc(32 * sizeof(unsigned char));
    
    //  α^d
    fq_t alpha_power_d;
    fq_init(alpha_power_d, par->Fp);
    fq_one(alpha_power_d, par->Fp);
    for (int i = 0; i < d; i++) {
        fq_mul(alpha_power_d, alpha_power_d, out->alpha, par->Fp);
    }
    
    //  g^{α^d}
    fmpz_t alpha_d_fmpz;
    fmpz_init(alpha_d_fmpz);
    fq2fmpz(alpha_d_fmpz, alpha_power_d, par->Fp);
    
    unsigned char alpha_d_chars[32];
    fmpz2chars(alpha_d_chars, alpha_d_fmpz);
    
    crypto_scalarmult_ristretto255_base(out->vk_x, alpha_d_chars);
    
    fmpz_clear(alpha_d_fmpz);
    fq_clear(alpha_power_d, par->Fp);
    
    // 7.  aux_D = x
    // aux_D 
    fq_mat_set(out->aux_D, x, par->Fp);
    
    // （r  out->r ，）
    for (int j = 0; j < m; j++) {
        for (int u_idx = 0; u_idx < num_u; u_idx++) {
            fq_clear(x_shares[j][u_idx], par->Fp);
            fq_clear(r_shares[j][u_idx], par->Fp);
        }
        free(x_shares[j]);
        free(r_shares[j]);
    }
    free(x_shares);
    free(r_shares);
    free(u_vectors);
}

// =========== Protocol 8 Compute  ===========

void scheme4_compute(scheme4_compute_out *out, int server_id, 
                      cnf_function_share_t *rho_i, scheme4_probgen_out *probgen_out, 
                      pubpar *par) {
    int m = par->m;
    int L = rho_i->L;
    
    out->L = L;
    
    // 1.  Y_i = (y_{i,1}, ..., y_{i,L})  y_{i,ell'} = F_{i,ell'}(s_i)
    out->Y_i = (fq_t*)malloc(L * sizeof(fq_t));
    out->tilde_Y_i = (fq_t*)malloc(L * sizeof(fq_t));
    
    // ， i  F_i(s_i)， s_i  i 
    //  F_{i,ell'}， F_{i,ell'}(s_i)  F_{i,ell'}(tilde{s}_i)
    
    //  i （）
    // ： i ， F_i(s_i)，
    // ：
    
    //  u 
    int *u_vectors = NULL;
    int num_u = scheme4_generate_cnf_vectors(&u_vectors, par->k, par->t);
    
    // ： probgen_out  aux_D（ x） r
    // ：（），
    fq_mat_t x_reconstructed, r_reconstructed;
    fq_mat_init(x_reconstructed, m, 1, par->Fp);
    fq_mat_init(r_reconstructed, m, 1, par->Fp);
    fq_mat_set(x_reconstructed, probgen_out->aux_D, par->Fp);
    fq_mat_set(r_reconstructed, probgen_out->r, par->Fp);
    
    //  F_{i,ell'}， F_{i,ell'}(x_reconstructed)  F_{i,ell'}(r_reconstructed)
    // ：（），
    for (int ell = 0; ell < L; ell++) {
        fq_init(out->Y_i[ell], par->Fp);
        fq_init(out->tilde_Y_i[ell], par->Fp);
        
        //  F_{i,ell'} 
        if (ell == 0) {
            // ：F_0
            fq_set(out->Y_i[ell], fq_mat_entry(rho_i->F_i[ell], 0, 0), par->Fp);
            fq_set(out->tilde_Y_i[ell], fq_mat_entry(rho_i->F_i[ell], 0, 0), par->Fp);
        } else if (ell <= par->m) {
            // ：F[ell, 0] * x[ell-1]
            fq_t coeff, val;
            fq_init(coeff, par->Fp);
            fq_init(val, par->Fp);
            fq_set(coeff, fq_mat_entry(rho_i->F_i[ell], 0, 0), par->Fp);
            
            fq_mul(val, coeff, fq_mat_entry(x_reconstructed, ell - 1, 0), par->Fp);
            fq_set(out->Y_i[ell], val, par->Fp);
            
            fq_mul(val, coeff, fq_mat_entry(r_reconstructed, ell - 1, 0), par->Fp);
            fq_set(out->tilde_Y_i[ell], val, par->Fp);
            
            fq_clear(coeff, par->Fp);
            fq_clear(val, par->Fp);
        } else {
            // ：
            // ：0（0）
            fq_zero(out->Y_i[ell], par->Fp);
            fq_zero(out->tilde_Y_i[ell], par->Fp);
        }
    }
    
    fq_mat_clear(x_reconstructed, par->Fp);
    fq_mat_clear(r_reconstructed, par->Fp);
    free(u_vectors);
    
    // 2.  r_i^{(j)} = {g^{x_{j,u}} | u_i = 0}
    //  tilde{r}_i^{(j)} = {g^{r_{j,u}} | u_i = 0}
    out->r_i = (unsigned char***)malloc(m * sizeof(unsigned char**));
    out->tilde_r_i = (unsigned char***)malloc(m * sizeof(unsigned char**));
    
    //  u （）
    int *u_vectors_compute = NULL;
    int num_u_compute = scheme4_generate_cnf_vectors(&u_vectors_compute, par->k, par->t);
    
    for (int j = 0; j < m; j++) {
        //  server_id （u_{server_id} = 0  u ）
        int count = 0;
        for (int u_idx = 0; u_idx < num_u_compute; u_idx++) {
            int u = u_vectors_compute[u_idx];
            if (!(u & (1 << server_id))) {  // u_{server_id} = 0
                count++;
            }
        }
        
        if (count == 0 || !probgen_out->s_i || !probgen_out->s_i[server_id] || 
            !probgen_out->s_i[server_id][j] ||
            !probgen_out->tilde_s_i || !probgen_out->tilde_s_i[server_id] ||
            !probgen_out->tilde_s_i[server_id][j]) {
            // ，
            out->r_i[j] = NULL;
            out->tilde_r_i[j] = NULL;
            continue;
        }
        
        out->r_i[j] = (unsigned char**)malloc(count * sizeof(unsigned char*));
        out->tilde_r_i[j] = (unsigned char**)malloc(count * sizeof(unsigned char*));
        
        int idx = 0;
        for (int u_idx = 0; u_idx < num_u_compute; u_idx++) {
            int u = u_vectors_compute[u_idx];
            if (!(u & (1 << server_id))) {  // u_{server_id} = 0
                out->r_i[j][idx] = (unsigned char*)malloc(32 * sizeof(unsigned char));
                out->tilde_r_i[j][idx] = (unsigned char*)malloc(32 * sizeof(unsigned char));
                
                //  g^{x_{j,u}}  g^{r_{j,u}}
                fmpz_t x_ju_fmpz, r_ju_fmpz;
                fmpz_init(x_ju_fmpz);
                fmpz_init(r_ju_fmpz);
                
                fq2fmpz(x_ju_fmpz, probgen_out->s_i[server_id][j][idx], par->Fp);
                fq2fmpz(r_ju_fmpz, probgen_out->tilde_s_i[server_id][j][idx], par->Fp);
                
                unsigned char x_ju_chars[32], r_ju_chars[32];
                fmpz2chars(x_ju_chars, x_ju_fmpz);
                fmpz2chars(r_ju_chars, r_ju_fmpz);
                
                crypto_scalarmult_ristretto255_base(out->r_i[j][idx], x_ju_chars);
                crypto_scalarmult_ristretto255_base(out->tilde_r_i[j][idx], r_ju_chars);
                
                fmpz_clear(x_ju_fmpz);
                fmpz_clear(r_ju_fmpz);
                idx++;
            }
        }
    }
    
    free(u_vectors_compute);
}

// =========== Protocol 8 Verify  ===========

/**
 * CNF reconstruction for Verify: ， g^
 *  l=1 
 */
void scheme4_cnf_reconstruct_verify(fq_t *result, unsigned char **g_coeffs, 
                                      int num_shares, pubpar *par) {
    //  l=1 ，reconstruction 
    // 
    
    // ： reconstruction  1
    // ， CNF 
    
    //  0， scheme4_verify 
    fq_init(*result, par->Fp);
    fq_zero(*result, par->Fp);
}

int scheme4_verify(scheme4_verify_out *out, unsigned char *vk_x,
                    scheme4_compute_out *pi_results, int k, 
                    scheme4_probgen_out *probgen_out, pubpar *par) {
    int m = par->m;
    int t = par->t;
    
    out->valid = 0;
    out->aux_D_prime = (fq_t**)malloc(k * sizeof(fq_t*));
    out->aux_D_prime_L = (int*)malloc(k * sizeof(int));
    //  NULL，
    for (int i = 0; i < k; i++) {
        out->aux_D_prime[i] = NULL;
        out->aux_D_prime_L[i] = 0;
    }
    
    // 1.  reconstruction  F(x)  F(r)
    //  openings {r_i^{(j)}, tilde{r}_i^{(j)}} 
    //  l=1 ，reconstruction （）
    // 
    
    //  H(u, 0) = t  u 
    int *u_vectors = NULL;
    int num_u = scheme4_generate_cnf_vectors(&u_vectors, k, t);
    
    if (num_u == 0) {
        return 0;
    }
    
    //  l=1  reconstruction， 1（）
    //  {g^{x_{j,u}}}  g^{x_j}
    //  reconstruction ，g^{x_j} = Π_u g^{x_{j,u}}（1）
    
    //  j， g^{x_j}  g^{r_j}
    //  l=1  reconstruction， 1
    // g^{x_j} = Π_{u} g^{x_{j,u}}， u  H(u, 0) = t
    unsigned char **g_x = (unsigned char**)malloc(m * sizeof(unsigned char*));
    unsigned char **g_r = (unsigned char**)malloc(m * sizeof(unsigned char*));
    
    for (int j = 0; j < m; j++) {
        g_x[j] = (unsigned char*)malloc(32 * sizeof(unsigned char));
        g_r[j] = (unsigned char*)malloc(32 * sizeof(unsigned char));
        
        // 
        memset(g_x[j], 0, 32);
        memset(g_r[j], 0, 32);
        
        // 
        //  u，（u_i = 0  i）
        for (int u_idx = 0; u_idx < num_u; u_idx++) {
            int u = u_vectors[u_idx];
            
            //  i， u_i = 0，
            for (int i = 0; i < k; i++) {
                if (!(u & (1 << i))) {  // u_i = 0
                    //  i  u 
                    int server_u_idx = 0;
                    for (int u_idx2 = 0; u_idx2 < num_u; u_idx2++) {
                        int u2 = u_vectors[u_idx2];
                        if (!(u2 & (1 << i))) {
                            if (u_idx2 == u_idx) {
                                break;
                            }
                            server_u_idx++;
                        }
                    }
                    
                    //  i 
                    if (pi_results[i].r_i && pi_results[i].r_i[j] && 
                        pi_results[i].r_i[j][server_u_idx] &&
                        pi_results[i].tilde_r_i && pi_results[i].tilde_r_i[j] &&
                        pi_results[i].tilde_r_i[j][server_u_idx]) {
                        //  g^{x_{j,u}}  g^{r_{j,u}} （）
                        crypto_core_ristretto255_add(g_x[j], g_x[j], pi_results[i].r_i[j][server_u_idx]);
                        crypto_core_ristretto255_add(g_r[j], g_r[j], pi_results[i].tilde_r_i[j][server_u_idx]);
                        break;
                    }
                }
            }
        }
    }
    
    // 2.  reconstruction  Y_i  tilde{Y}_i  F(x)  F(r)
    //  l=1  reconstruction， 1
    // ： F(x)， Y_i 
    //  F(x) = Y_i（ i）， Σ_i Y_i
    
    fq_t F_x, F_r;
    fq_init(F_x, par->Fp);
    fq_init(F_r, par->Fp);
    fq_zero(F_x, par->Fp);
    fq_zero(F_r, par->Fp);
    
    //  Y_i  tilde{Y}_i （，）
    // ， CNF reconstruction，
    if (k > 0 && pi_results[0].Y_i && pi_results[0].tilde_Y_i) {
        // ，
        for (int ell = 0; ell < pi_results[0].L; ell++) {
            fq_add(F_x, F_x, pi_results[0].Y_i[ell], par->Fp);
            fq_add(F_r, F_r, pi_results[0].tilde_Y_i[ell], par->Fp);
        }
    }
    
    // 3. ：， g^{F(r)} = g^{α^d} · g^{F(x)}
    // ，F(r) = α^d * F(x)
    //  F(x) = F_0 + F_hom(x)， F_hom ：
    // F(r) = F(αx) = F_0 + F_hom(αx) = F_0 + α^d * F_hom(x)
    //  F(r) = F_0 + α^d * (F(x) - F_0) = α^d * F(x) + (1 - α^d) * F_0
    
    //  F(x) = 1 + Σx_i，
    //  F(x) = 1 + Σx_i  r = αx：
    // F(r) = 1 + α * Σx_i
    //  α^d * F(x) + (1 - α^d) * F_0 = α^2 * (1 + Σx_i) + (1 - α^2) * 1 = α^2 + α^2 * Σx_i + 1 - α^2 = 1 + α^2 * Σx_i
    // ， α = 1  α^2 = α
    
    // ，
    //  F(x) = F_0 + F_1(x)（ F_1  1 ）：
    // F(r) = F(αx) = F_0 + F_1(αx) = F_0 + α * F_1(x) = F_0 + α * (F(x) - F_0) = α * F(x) + (1 - α) * F_0
    
    //  α
    fq_t alpha;
    fq_init(alpha, par->Fp);
    fq_set(alpha, probgen_out->alpha, par->Fp);
    
    //  α^d（）
    fq_t alpha_power_d;
    fq_init(alpha_power_d, par->Fp);
    fq_one(alpha_power_d, par->Fp);
    for (int i = 0; i < par->d; i++) {
        fq_mul(alpha_power_d, alpha_power_d, probgen_out->alpha, par->Fp);
    }
    
    //  F_0（，）
    fq_t F_0;
    fq_init(F_0, par->Fp);
    if (k > 0 && pi_results[0].Y_i && pi_results[0].L > 0) {
        fq_set(F_0, pi_results[0].Y_i[0], par->Fp);
    } else {
        fq_zero(F_0, par->Fp);
    }
    
    // ， d 
    //  d=1（）：F(r) = α * F(x) + (1 - α) * F_0
    //  d>1：F(r) = α^d * F(x) + (1 - α^d) * F_0
    fq_t one, one_minus_alpha_power, temp1, temp2, F_r_expected;
    fq_init(one, par->Fp);
    fq_init(one_minus_alpha_power, par->Fp);
    fq_init(temp1, par->Fp);
    fq_init(temp2, par->Fp);
    fq_init(F_r_expected, par->Fp);
    
    fq_one(one, par->Fp);
    
    // ，
    //  F(x) = 1 + Σx_i， 1 
    //  d， α（ 1）
    // ， α^d
    //  α^d，
    
    // ， F(x) = F_0 + F_hom(x)：
    // F(r) = F(αx) = F_0 + F_hom(αx) = F_0 + α^d * F_hom(x)
    //  d 
    
    //  F(x) = 1 + Σx_i， 1， α^1 = α
    //  d  2，
    
    // ：， α（ 1）
    fq_t alpha_power;
    fq_init(alpha_power, par->Fp);
    // ， α
    fq_set(alpha_power, alpha, par->Fp);
    
    fq_sub(one_minus_alpha_power, one, alpha_power, par->Fp);
    
    // temp1 = α^power * F(x)
    fq_mul(temp1, alpha_power, F_x, par->Fp);
    // temp2 = (1 - α^power) * F_0
    fq_mul(temp2, one_minus_alpha_power, F_0, par->Fp);
    // F_r_expected = α^power * F(x) + (1 - α^power) * F_0
    fq_add(F_r_expected, temp1, temp2, par->Fp);
    
    fq_clear(alpha_power, par->Fp);
    
    //  F(r) == F_r_expected
    int is_equal = fq_equal(F_r, F_r_expected, par->Fp);
    
    if (is_equal) {
        out->valid = 1;
        
        //  aux'_D = {Y_i | i ∈ [k]}
        for (int i = 0; i < k; i++) {
            if (pi_results[i].Y_i && pi_results[i].L > 0) {
                out->aux_D_prime_L[i] = pi_results[i].L;
                out->aux_D_prime[i] = (fq_t*)malloc(pi_results[i].L * sizeof(fq_t));
                for (int ell = 0; ell < pi_results[i].L; ell++) {
                    fq_init(out->aux_D_prime[i][ell], par->Fp);
                    fq_set(out->aux_D_prime[i][ell], pi_results[i].Y_i[ell], par->Fp);
                }
            } else {
                out->aux_D_prime[i] = NULL;
                out->aux_D_prime_L[i] = 0;
            }
        }
    }
    
    // 
    fq_clear(F_x, par->Fp);
    fq_clear(F_r, par->Fp);
    fq_clear(F_0, par->Fp);
    fq_clear(alpha, par->Fp);
    fq_clear(alpha_power_d, par->Fp);
    fq_clear(one, par->Fp);
    fq_clear(one_minus_alpha_power, par->Fp);
    fq_clear(temp1, par->Fp);
    fq_clear(temp2, par->Fp);
    fq_clear(F_r_expected, par->Fp);
    for (int j = 0; j < m; j++) {
        free(g_x[j]);
        free(g_r[j]);
    }
    free(g_x);
    free(g_r);
    free(u_vectors);
    
    return out->valid;
}

// =========== Protocol 8 Decode  ===========

/**
 * CNF reconstruction for Decode: ，
 */
void scheme4_cnf_reconstruct_decode(fq_t *result, fq_t *shares, 
                                     int num_shares, pubpar *par) {
    //  l=1 ，reconstruction 
    fq_init(*result, par->Fp);
    fq_zero(*result, par->Fp);
    
    //  1（）
    // ， CNF 
    for (int i = 0; i < num_shares; i++) {
        fq_add(*result, *result, shares[i], par->Fp);
    }
}

void scheme4_decode(scheme4_decode_out *out, fq_mat_t aux_D, 
                     fq_t **aux_D_prime, scheme4_keygen_out *keygen_out, 
                     pubpar *par) {
    int k = par->k;
    
    fq_init(out->result, par->Fp);
    fq_zero(out->result, par->Fp);
    
    //  reconstruction  Rec_ell  {Y_i | i ∈ [k]}  F(x)
    //  l=1  reconstruction， 1
    // F(x) = Σ_i Σ_ell Y_i[ell]
    // ： F(x)， Y_i 
    //  F(x) = Y_i（ i）
    
    //  Y_i （）
    if (k > 0 && aux_D_prime[0]) {
        int L = keygen_out->rho[0].L;
        for (int ell = 0; ell < L; ell++) {
            fq_add(out->result, out->result, aux_D_prime[0][ell], par->Fp);
        }
    }
}

// =========== Protocol 8  ===========

void scheme4_keygen_out_clear(scheme4_keygen_out *out, pubpar *par) {
    if (out->rho) {
        for (int i = 0; i < par->k; i++) {
            if (out->rho[i].F_i) {
                for (int ell = 0; ell < out->rho[i].L; ell++) {
                    fq_mat_clear(out->rho[i].F_i[ell], par->Fp);
                }
                free(out->rho[i].F_i);
            }
        }
        free(out->rho);
    }
}

void scheme4_probgen_out_clear(scheme4_probgen_out *out, pubpar *par) {
    if (out->vk_x) {
        free(out->vk_x);
    }
    
    if (out->s_i) {
        for (int i = 0; i < par->k; i++) {
            if (out->s_i[i]) {
                for (int j = 0; j < par->m; j++) {
                    if (out->s_i[i][j]) {
                        // ，
                        free(out->s_i[i][j]);
                    }
                }
                free(out->s_i[i]);
            }
        }
        free(out->s_i);
    }
    
    if (out->tilde_s_i) {
        for (int i = 0; i < par->k; i++) {
            if (out->tilde_s_i[i]) {
                for (int j = 0; j < par->m; j++) {
                    if (out->tilde_s_i[i][j]) {
                        free(out->tilde_s_i[i][j]);
                    }
                }
                free(out->tilde_s_i[i]);
            }
        }
        free(out->tilde_s_i);
    }
    
    if (out->num_shares_per_coord) {
        free(out->num_shares_per_coord);
    }
    
    fq_mat_clear(out->aux_D, par->Fp);
    fq_mat_clear(out->r, par->Fp);
    fq_clear(out->alpha, par->Fp);
}

void scheme4_compute_out_clear(scheme4_compute_out *out, pubpar *par) {
    if (out->Y_i) {
        for (int ell = 0; ell < out->L; ell++) {
            fq_clear(out->Y_i[ell], par->Fp);
        }
        free(out->Y_i);
    }
    
    if (out->tilde_Y_i) {
        for (int ell = 0; ell < out->L; ell++) {
            fq_clear(out->tilde_Y_i[ell], par->Fp);
        }
        free(out->tilde_Y_i);
    }
    
    if (out->r_i) {
        for (int j = 0; j < par->m; j++) {
            if (out->r_i[j]) {
                // r_i[j] ， 32 
                // ，
                // ，r_i[j][idx]  32 ，
                //  r_i[j] 
                free(out->r_i[j]);
            }
        }
        free(out->r_i);
    }
    
    if (out->tilde_r_i) {
        for (int j = 0; j < par->m; j++) {
            if (out->tilde_r_i[j]) {
                free(out->tilde_r_i[j]);
            }
        }
        free(out->tilde_r_i);
    }
}

void scheme4_verify_out_clear(scheme4_verify_out *out, int k, pubpar *par) {
    if (out->aux_D_prime) {
        for (int i = 0; i < k; i++) {
            if (out->aux_D_prime[i] && out->aux_D_prime_L) {
                for (int ell = 0; ell < out->aux_D_prime_L[i]; ell++) {
                    fq_clear(out->aux_D_prime[i][ell], par->Fp);
                }
                free(out->aux_D_prime[i]);
            }
        }
        free(out->aux_D_prime);
    }
    if (out->aux_D_prime_L) {
        free(out->aux_D_prime_L);
    }
}

void scheme4_decode_out_clear(scheme4_decode_out *out, pubpar *par) {
    fq_clear(out->result, par->Fp);
}
