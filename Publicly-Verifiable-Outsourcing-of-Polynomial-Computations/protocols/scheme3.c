#include "scheme3.h"

// ========== Paillier ==========

// ========== Paillier/HSS  (Based on GroupHSS/HSSElg.cpp) ==========

// Helper to get random state
void init_rand_state(flint_rand_t state) {
    flint_randinit(state);
}

void clear_rand_state(flint_rand_t state) {
    flint_randclear(state);
}

/**
 * Paillier (Elgamal_Gen)
 */
void paillier_keygen(paillier_pk_t *pk, fmpz_t sk, int security_bits) {
    fmpz_t p, q;
    fmpz_init(p);
    fmpz_init(q);
    
    // GroupHSS use fixed large primes. We use the same for consistency/simplicity in this port.
    // p = ...
     fmpz_set_str(p, "2410312426921032588580116606028314112912093247945688951359675039065257391591803200669085024107346049663448766280888004787862416978794958324969612987890774651455213339381625224770782077917681499676845543137387820057597345857904599109461387122099507964997815641342300677629473355281617428411794163967785870370368969109221591943054232011562758450080579587850900993714892283476646631181515063804873375182260506246992837898705971012525843324401232986857004760339321639", 10);
    fmpz_set_str(q, "2410312426921032588580116606028314112912093247945688951359675039065257391591803200669085024107346049663448766280888004787862416978794958324969612987890774651455213339381625224770782077917681499676845543137387820057597345857904599109461387122099507964997815641342300677629473355281617428411794163967785870370368969109221591943054232011562758450080579587850900993714892283476646631181515063804873375182260506246992837898705971012525843324401232986857004760339319223", 10);
    
    // N = p * q
    fmpz_init(pk->n);
    fmpz_mul(pk->n, p, q);
    
    // N^2
    fmpz_init(pk->n2);
    fmpz_mul(pk->n2, pk->n, pk->n);
    
    // g
    fmpz_init(pk->g);
    
    flint_rand_t state;
    flint_randinit(state);
    
    fmpz_t r_rand;
    fmpz_init(r_rand);
    fmpz_randm(r_rand, state, pk->n2); // Use FLINT random
    fmpz_set(pk->g, r_rand);
    
    // f = n + 1
    fmpz_init(pk->f);
    fmpz_add_ui(pk->f, pk->n, 1);
    
    // sk = random bits
    fmpz_randm(sk, state, pk->n); // Random in N
    
    // h = g^sk mod n^2
    fmpz_init(pk->h);
    fmpz_powm(pk->h, pk->g, sk, pk->n2);
    
    fmpz_init(pk->mu); // Unused in this scheme potentially
    
    flint_randclear(state);
    fmpz_clear(p);
    fmpz_clear(q);
    fmpz_clear(r_rand);
}

/**
 * Elgamal Encryption (Elgamal_Enc)
 * ct = [g^r, h^r * (1 + N*x)] mod N^2
 */
void elgamal_encrypt(fmpz_t c1, fmpz_t c2, paillier_pk_t *pk, fmpz_t x) {
    fmpz_t r, temp;
    fmpz_init(r);
    fmpz_init(temp);
    
    flint_rand_t state;
    flint_randinit(state);
    
    // Random r in [0, N-1]
    fmpz_randm(r, state, pk->n);
    
    // c1 = g^r mod N^2
    fmpz_powm(c1, pk->g, r, pk->n2);
    
    // term2 = h^r mod N^2
    fmpz_powm(c2, pk->h, r, pk->n2);
    
    // term3 = 1 + N*x mod N^2
    fmpz_mul(temp, pk->n, x);
    fmpz_add_ui(temp, temp, 1);
    fmpz_mod(temp, temp, pk->n2);
    
    // c2 = h^r * (1+N*x) mod N^2
    fmpz_mul(c2, c2, temp);
    fmpz_mod(c2, c2, pk->n2);
    
    flint_randclear(state);
    fmpz_clear(r);
    fmpz_clear(temp);
}

/**
 * Elgamal Secret Key Encryption (Elgamal_skEnc)
 * ct = [g^r * (1 - N*x), h^r] mod N^2
 * Note: GroupHSS implementation uses this for the second component of HSS Input
 */
void elgamal_sk_encrypt(fmpz_t c1, fmpz_t c2, paillier_pk_t *pk, fmpz_t x) {
    fmpz_t r, temp;
    fmpz_init(r);
    fmpz_init(temp);
    
    flint_rand_t state;
    flint_randinit(state);
    
    // Random r
    fmpz_randm(r, state, pk->n);
    
    // c1 part 1 = g^r
    fmpz_powm(c1, pk->g, r, pk->n2);
    
    // part 2 = 1 - N*x = 1 + N*(-x)
    fmpz_neg(temp, x);
    fmpz_mul(temp, temp, pk->n);
    fmpz_add_ui(temp, temp, 1);
    fmpz_mod(temp, temp, pk->n2); // Canonicalize
    
    // c1 = g^r * (1-Nx)
    fmpz_mul(c1, c1, temp);
    fmpz_mod(c1, c1, pk->n2);
    
    // c2 = h^r
    fmpz_powm(c2, pk->h, r, pk->n2);
    
    flint_randclear(state);
    fmpz_clear(r);
    fmpz_clear(temp);
}

// Wrapper for old interface if needed, or we implement correct HSS_Input directly
// We don't need paillier_encrypt/add/mul exposed likely, but keeping interface clean.

// Helper to init/clear 2D array in struct is annoying in C.
// We assume caller or init functions handle it.


// ========== HSS (Based on GroupHSS) ==========

/**
 * HSS (HSS_Gen)
 */
void hss_gen(paillier_pk_t *pk, fmpz_t sk, hss_ek_t *ek0, hss_ek_t *ek1, int security_bits) {
    // Generate Paillier keys (and implicit sk = d)
    paillier_keygen(pk, sk, security_bits);
    
    fmpz_init(ek0->ek);
    fmpz_init(ek1->ek);
    
    flint_rand_t state;
    flint_randinit(state);
    
    // Random bits for ek0 (random in N)
    fmpz_randm(ek0->ek, state, pk->n);
    
    // ek1 = ek0 + sk mod N? Or just add?
    // HSSElg line 63: add(ek1, ek0, s) -> ek1 = ek0 + s (ZZ addition)
    // NOTE: GroupHSS does not seem to mod by N in HSS_Gen.
    // BUT memory values are usually mod N.
    // If we want consistency, we should maybe Mod N?
    // Let's follow HSS_Gen: add(ek1, ek0, s).
    // Since sk and ek0 are in [0, N-1], sum is < 2N.
    // We will stick to strict implementation from HSSElg.cpp
    fmpz_add(ek1->ek, ek0->ek, sk);
    // fmpz_mod(ek1->ek, ek1->ek, pk->n); // Unleashed?
    
    flint_randclear(state);
}

/**
 * HSS (HSS_Input)
 * ct = [Elgamal_Enc(x), Elgamal_skEnc(x)]
 */
void hss_input(paillier_ct_t *ct, paillier_pk_t *pk, fq_t x, pubpar *par) {
    // fq_tfmpz_t
    fmpz_t x_fmpz;
    fmpz_init(x_fmpz);
    fq2fmpz(x_fmpz, x, par->Fp);
    
    // Initialize output
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
            fmpz_init(ct->c[i][j]);
        }
    }
    
    // ct[0] = Elgamal_Enc(x)
    elgamal_encrypt(ct->c[0][0], ct->c[0][1], pk, x_fmpz);
    
    // ct[1] = Elgamal_skEnc(x)
    elgamal_sk_encrypt(ct->c[1][0], ct->c[1][1], pk, x_fmpz);
    
    fmpz_clear(x_fmpz);
}

// ========== Paillier HSS  / Core Logic ==========

/**
 * HSS DDLog (HSS_DDLog)
 * Solves x from (1+Nx) = M mod N^2
 * But input g is whatever we have.
 * z = (g/N) * (g%N)^-1 mod N
 */
void hss_ddlog(fmpz_t z, paillier_pk_t *pk, fmpz_t g) {
    fmpz_t h1, h, temp1, rem;
    fmpz_init(h1);
    fmpz_init(h);
    fmpz_init(temp1);
    fmpz_init(rem);
    
    // DivRem: g = h1 * N + h  => h = g % N, h1 = g / N
    fmpz_tdiv_qr(h1, h, g, pk->n);
    
    // temp1 = h^-1 mod N
    fmpz_invmod(temp1, h, pk->n);
    
    // z = h1 * temp1 mod N
    fmpz_mul(z, h1, temp1);
    fmpz_mod(z, z, pk->n);
    
    fmpz_clear(h1);
    fmpz_clear(h);
    fmpz_clear(temp1);
    fmpz_clear(rem);
}

/**
 * PRF_ZZ
 * Returns random integer in [0, N-1] deterministically from key
 */
void prf_zz(fmpz_t res, int *prf_key, fmpz_t modulus) {
    // Generate deterministic randomness from prf_key
    unsigned char seed[32];
    crypto_hash_sha256(seed, (unsigned char*)prf_key, sizeof(int));
    
    // If N is 3072 bits, we need more bytes.
    // HSSElg uses NTL PRF logic.
    // We will simulate a stream cipher.
    
    // Simple implementation: Key expansion using ChaCha20 or just hashing counter
    // For large N, we need multiple blocks.
    int n_bytes = fmpz_sizeinbase(modulus, 2) / 8 + 64; // + safety
    unsigned char *buf = malloc(n_bytes);
    
    // Use randombytes_buf_deterministic
    randombytes_buf_deterministic(buf, n_bytes, seed);
    
    // Convert to fmpz
    fmpz_t val, base;
    fmpz_init(val);
    fmpz_init(base);
    fmpz_set_ui(base, 1);
    fmpz_zero(res);
    
    for(int i=0; i<n_bytes; i++) {
        fmpz_set_ui(val, buf[i]);
        fmpz_mul(val, val, base);
        fmpz_add(res, res, val);
        fmpz_mul_ui(base, base, 256);
    }
    
    fmpz_mod(res, res, modulus);
    
    fmpz_clear(val);
    fmpz_clear(base);
    free(buf);
    
    (*prf_key)++; // Increment key
}

/**
 * HSS_Mul implementation
 * Mz = HSS_Mul(b, pk, Ix, My, prf_key)
 * Mz[0] = DDLog( Ix[0][1]^My[0] * Ix[0][0]^-My[1] ) + PRF
 * Mz[1] = DDLog( Ix[1][1]^My[0] * Ix[1][0]^-My[1] ) + PRF
 */
void hss_mul_paillier(hss_mv_t *result, int server_id, 
                      paillier_pk_t *pk, paillier_ct_t *Ix,
                      hss_mv_t *My, int *prf_key) {
    fmpz_t temp1, temp2, term;
    fmpz_t prf_val;
    fmpz_init(temp1);
    fmpz_init(temp2);
    fmpz_init(term);
    fmpz_init(prf_val);
    
    fmpz_init(result->mv[0]);
    fmpz_init(result->mv[1]);
    
    // Loop for 2 components of result (Share has 2 parts)
    for(int k=0; k<2; k++) {
        // temp1 = Ix[k][1] ^ My[0] mod N^2
        fmpz_powm(temp1, Ix->c[k][1], My->mv[0], pk->n2);
        
        // temp2 = Ix[k][0] ^ (-My[1]) mod N^2
        fmpz_t neg_my1;
        fmpz_init(neg_my1);
        fmpz_neg(neg_my1, My->mv[1]); 
        // Exponent can be negative. fmpz_powm handles positive exponents usually? 
        // FLINT powm takes integer exponent.
        // If negative, we compute inverse.
        // Actually fmpz_powm expects positive?
        // Let's optimize: temp2 = (Ix[k][0] ^ My[1])^-1
        fmpz_powm(temp2, Ix->c[k][0], My->mv[1], pk->n2);
        fmpz_invmod(temp2, temp2, pk->n2);
        
        fmpz_clear(neg_my1); // logic simplified
        
        // term = temp1 * temp2
        fmpz_mul(term, temp1, temp2);
        fmpz_mod(term, term, pk->n2);
        
        // DDLog
        hss_ddlog(result->mv[k], pk, term);
        
        // Add PRF
        prf_zz(prf_val, prf_key, pk->n);
        fmpz_add(result->mv[k], result->mv[k], prf_val);
        // Result is in Z (integer). Usually shares are mod N?
        // GroupHSS uses "ZZ" and DDLog is mod N.
        // But HSS_Mul adds PRF (mod N).
        // Does it mod N afterwards?
        // HSSElg.cpp: Mz[0] = PRF(...) + Mz[0]; No explicit mod.
        // But DDLog result is < N. PRF is < N. Sum < 2N.
        // If we don't mod, the shares grow?
        // HSS_AddMemory adds them.
        // Eventually Decrypt subtracts them.
        // (M0 + M1) mod N?
        // Usually shares are in Z_N.
        // I will mod N to be safe and consistent with "AddMemory".
        fmpz_mod(result->mv[k], result->mv[k], pk->n);
    }
    
    fmpz_clear(temp1);
    fmpz_clear(temp2);
    fmpz_clear(term);
    fmpz_clear(prf_val);
}

void hss_add_memory(hss_mv_t *res, paillier_pk_t *pk, hss_mv_t *mx, hss_mv_t *my) {
    if(res != mx) { // Support in-place if res==mx
        fmpz_set(res->mv[0], mx->mv[0]);
        fmpz_set(res->mv[1], mx->mv[1]);
    }
    fmpz_add(res->mv[0], res->mv[0], my->mv[0]);
    fmpz_mod(res->mv[0], res->mv[0], pk->n);
    
    fmpz_add(res->mv[1], res->mv[1], my->mv[1]);
    fmpz_mod(res->mv[1], res->mv[1], pk->n);
}


/**
 * HSS
 * linear combination for Protocol 9
 */
void hss_evaluate(hss_mv_t *result, int server_id, hss_ek_t *ek_b,
                  paillier_pk_t *pk, fq_mat_t F, paillier_ct_t *ct_inputs, pubpar *par) {
    int m = par->m;
    int prf_key = 0; // Initialize PRF key
    
    // Initialize result to 0
    fmpz_init(result->mv[0]);
    fmpz_init(result->mv[1]);
    fmpz_zero(result->mv[0]);
    fmpz_zero(result->mv[1]);
    
    // Unity share (represents "1")
    // M = (b, ek_b)
    hss_mv_t unity;
    fmpz_init(unity.mv[0]);
    fmpz_init(unity.mv[1]);
    fmpz_set_si(unity.mv[0], server_id);
    fmpz_set(unity.mv[1], ek_b->ek);
    
    hss_mv_t term_share, tmp_share;
    // Init temp shares
    for(int k=0; k<2; k++) {
        fmpz_init(term_share.mv[k]);
        fmpz_init(tmp_share.mv[k]);
    }
    
    // Loop over input terms (linear)
    // F[0] = const, F[k] for k in 1..m
    
    fq_t coeff;
    fq_init(coeff, par->Fp);
    fmpz_t coeff_fmpz;
    fmpz_init(coeff_fmpz);
    
    for(int k=0; k<m; k++) {
        // Coeff for x_k is F[k+1]
        // If F has enough rows
         if(k+1 < fq_mat_nrows(F, par->Fp)) {
            fq_set(coeff, fq_mat_entry(F, k+1, 0), par->Fp);
        } else {
            fq_zero(coeff, par->Fp);
        }
        fq2fmpz(coeff_fmpz, coeff, par->Fp);
        
        if(fmpz_is_zero(coeff_fmpz)) continue;
        
        // Compute share of x_k: Mul(unity, ct[k])
        hss_mul_paillier(&term_share, server_id, pk, &ct_inputs[k], &unity, &prf_key);
        
        // Multiply by coeff (scalar mul)
        // term_share * coeff
        fmpz_mul(term_share.mv[0], term_share.mv[0], coeff_fmpz);
        fmpz_mod(term_share.mv[0], term_share.mv[0], pk->n);
        fmpz_mul(term_share.mv[1], term_share.mv[1], coeff_fmpz);
        fmpz_mod(term_share.mv[1], term_share.mv[1], pk->n);
        
        // Add to result
        hss_add_memory(result, pk, result, &term_share);
    }
    
    // Constant term F[0]
    fq_set(coeff, fq_mat_entry(F, 0, 0), par->Fp);
    fq2fmpz(coeff_fmpz, coeff, par->Fp);
    
    if(!fmpz_is_zero(coeff_fmpz)) {
        // Share of const C is just C * Unity_Share
        // unity = (b, ek)
        fmpz_mul(term_share.mv[0], unity.mv[0], coeff_fmpz);
        fmpz_mod(term_share.mv[0], term_share.mv[0], pk->n);
        fmpz_mul(term_share.mv[1], unity.mv[1], coeff_fmpz);
        fmpz_mod(term_share.mv[1], term_share.mv[1], pk->n);
        
        hss_add_memory(result, pk, result, &term_share);
    }
    
    fq_clear(coeff, par->Fp);
    fmpz_clear(coeff_fmpz);
    fmpz_clear(unity.mv[0]); fmpz_clear(unity.mv[1]);
    fmpz_clear(term_share.mv[0]); fmpz_clear(term_share.mv[1]);
    fmpz_clear(tmp_share.mv[0]); fmpz_clear(tmp_share.mv[1]);
}

/**
 * HSS
 * result = mv1 - mv0 mod n
 * Note: HSSElg says sub(z, Mx1[0], Mx0[0]);
 * Wait, HSSElg uses vector index 0 as result?
 * HSS_Dec(ZZ &z, const HSS_MV &Mx0, const HSS_MV &Mx1) { sub(z, Mx1[0], Mx0[0]); }
 * HSS_MV has 2 elements.
 * It seems element 0 carries the value?
 * Let's check HSS_Mul again.
 * Mz[0] = ... + PRF
 * Mz[1] = ... + PRF
 * Unity share: M[0]=b, M[1]=ek.
 * If b=0 (Server 0), M[0]=0.
 * If b=1 (Server 1), M[0]=1.
 * M[1] is ek.
 * The value seems to be accumulated in both components but extraction uses component 0?
 * Let's trust HSS_Dec implementation: sub(z, Mx1[0], Mx0[0]).
 */
void hss_decrypt(fq_t result, paillier_pk_t *pk, hss_mv_t *mv0, hss_mv_t *mv1, pubpar *par) {
    fmpz_t sum;
    fmpz_init(sum);
    
    // z = mv1[0] - mv0[0] mod N
    fmpz_sub(sum, mv1->mv[0], mv0->mv[0]);
    fmpz_mod(sum, sum, pk->n);
    
    // Convert to fq
    // Ensure it's modulo p if p < n
    // Paillier space is large (N), result x is small (Fp).
    // We assume result fits in p.
    fmpz_mod(sum, sum, par->p); // mod p just in case
    fq_set_fmpz(result, sum, par->Fp);
    
    fmpz_clear(sum);
}

// ========== Protocol 9  ==========

/**
 * Protocol 9 KeyGen
 */
void scheme3_keygen(scheme3_keygen_out *out, fq_mat_t F, pubpar *par) {
    // HSS
    int security_bits = 1024;
    fmpz_init(out->sk);
    hss_gen(&out->pk, out->sk, &out->ek[0], &out->ek[1], security_bits);
}

/**
 * Protocol 9 ProbGen
 */
void scheme3_probgen(scheme3_probgen_out *out, paillier_pk_t *pk, 
                       fq_mat_t x, pubpar *par) {
    int d = par->d;
    int m = par->m;
    
    // 
    out->vk_x = (unsigned char**)malloc((d+1) * sizeof(unsigned char*));
    for(int j = 0; j <= d; j++) {
        out->vk_x[j] = (unsigned char*)malloc(32 * sizeof(unsigned char));
    }
    
    out->ct_inputs = (paillier_ct_t*)malloc((d+1) * m * sizeof(paillier_ct_t));
    
    out->v_values = (fq_t*)malloc(d * sizeof(fq_t));
    for(int j = 0; j < d; j++) {
        fq_init(out->v_values[j], par->Fp);
    }
    
    for(int j = 0; j < d; j++) {
        uint32_t rand_val = randombytes_random();
        fq_set_ui(out->v_values[j], (rand_val % 1000) + 2, par->Fp);
    }
    
    // v1
    fq_t one;
    fq_init(one, par->Fp);
    fq_one(one, par->Fp);
    
    for(int i = 0; i < d; i++) {
        if(fq_equal(out->v_values[i], one, par->Fp)) {
            uint32_t rand_val = randombytes_random();
            fq_set_ui(out->v_values[i], (rand_val % 1000) + 100, par->Fp);
        }
        
        for(int j = i + 1; j < d; j++) {
            if(fq_equal(out->v_values[i], out->v_values[j], par->Fp)) {
                uint32_t rand_val = randombytes_random();
                fq_set_ui(out->v_values[j], (rand_val % 1000) + 100, par->Fp);
                j = i;
            }
        }
    }
    
    fq_clear(one, par->Fp);
    
    // ct_0 = HSS.Input(x)  (v_0 = 1)
    
    fq_t scaled_val;
    fq_init(scaled_val, par->Fp);
    fq_t x_elem;
    fq_init(x_elem, par->Fp);
    
    for(int j = 0; j <= d; j++) {
        for(int k = 0; k < m; k++) {
            fq_set(x_elem, fq_mat_entry(x, k, 0), par->Fp);
            
            if(j == 0) {
                // v_0 = 1
                hss_input(&out->ct_inputs[j*m + k], pk, x_elem, par);
            } else {
                // v_j * x_k
                fq_mul(scaled_val, out->v_values[j-1], x_elem, par->Fp);
                hss_input(&out->ct_inputs[j*m + k], pk, scaled_val, par);
            }
        }
    }
    
    fq_clear(scaled_val, par->Fp);
    fq_clear(x_elem, par->Fp);
    
    fq_t *lambda_prime = (fq_t*)malloc((d+1) * sizeof(fq_t));
    for(int j = 0; j <= d; j++) {
        fq_init(lambda_prime[j], par->Fp);
    }
    scheme3_compute_lambda_prime(lambda_prime, out->v_values, d, par);
    
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
        fq_clear(lambda_prime[j], par->Fp);
    }
    free(lambda_prime);
}

/**
 * Protocol 9 Compute
 * bFHSS
 */
void scheme3_compute(scheme3_compute_out *out, int server_id, 
                       hss_ek_t *ek_b, paillier_pk_t *pk, fq_mat_t F,
                       paillier_ct_t *ct_inputs, pubpar *par) {
    int d = par->d;
    
    // 
    out->mv = (hss_mv_t*)malloc((d+1) * sizeof(hss_mv_t));
    
    int m = par->m;
    for(int j = 0; j <= d; j++) {
        hss_evaluate(&out->mv[j], server_id, ek_b, pk, F, &ct_inputs[j*m], par);
    }
}

/**
 * Protocol 9 Verify
 * 
 */
int scheme3_verify(scheme3_verify_out *out, unsigned char **vk_x,
                     scheme3_compute_out *compute_results, fq_t *v_values,
                     fq_mat_t F, paillier_pk_t *pk, pubpar *par) {
    int d = par->d;
    
    // 
    fq_mat_init(out->aux_R_prime, 1, 1, par->Fp);
    out->valid = 0;
    
    // 1. F(x), F(v_1*x), ..., F(v_d*x)
    fq_t *f_values = (fq_t*)malloc((d+1) * sizeof(fq_t));
    for(int j = 0; j <= d; j++) {
        fq_init(f_values[j], par->Fp);
        // f_values[j] = HSS.Dec(mv0[j], mv1[j])
        hss_decrypt(f_values[j], pk, &compute_results[0].mv[j], 
                   &compute_results[1].mv[j], par);
    }
    
    // 2. 
    fq_t *lambda_prime = (fq_t*)malloc((d+1) * sizeof(fq_t));
    for(int j = 0; j <= d; j++) {
        fq_init(lambda_prime[j], par->Fp);
    }
    scheme3_compute_lambda_prime(lambda_prime, v_values, d, par);
    
    fq_t F_0;
    fq_init(F_0, par->Fp);
    fq_set(F_0, fq_mat_entry(F, 0, 0), par->Fp);
    
    fq_t *exponents = (fq_t*)malloc((d+1) * sizeof(fq_t));
    for(int i = 0; i <= d; i++) {
        fq_init(exponents[i], par->Fp);
        fq_sub(exponents[i], f_values[i], F_0, par->Fp);
    }
    
    // 
    unsigned char result[32];
    unsigned char temp[32];
    memset(result, 0, 32);
    
    for(int i = 0; i <= d; i++) {
        fmpz_t exp_fmpz;
        fmpz_init(exp_fmpz);
        fq2fmpz(exp_fmpz, exponents[i], par->Fp);
        
        unsigned char exp_chars[32];
        fmpz2chars(exp_chars, exp_fmpz);
        
        int ms = crypto_scalarmult_ristretto255(temp, exp_chars, vk_x[i]);
        if (ms == 0) {
            crypto_core_ristretto255_add(result, result, temp);
        } else {
            // printf("DEBUG Verify: scalarmult failed for i=%d\n", i);
        }
        
        fmpz_clear(exp_fmpz);
    }
    
    // printf("DEBUG Verify: Final Result = ");
    // for(int b=0; b<32; b++) printf("%02x", result[b]);
    // printf("\n");
    
    // 
    unsigned char identity[32];
    memset(identity, 0, 32);
    
    if (memcmp(result, identity, 32) == 0) {
        out->valid = 1;
        // aux'_R = f(0) = f_values[0]
        fq_mat_entry_set(out->aux_R_prime, 0, 0, f_values[0], par->Fp);
    } else {
        out->valid = 0;
    }
    
    // 
    fq_clear(F_0, par->Fp);
    for(int i = 0; i <= d; i++) {
        fq_clear(exponents[i], par->Fp);
        fq_clear(lambda_prime[i], par->Fp);
        fq_clear(f_values[i], par->Fp);
    }
    free(exponents);
    free(lambda_prime);
    free(f_values);
    
    return out->valid;
}

/**
 * Protocol 9 Reconstruct
 * 
 */
void scheme3_reconstruct(scheme3_reconstruct_out *out, fq_mat_t aux_R_prime, pubpar *par) {
    // 
    fq_mat_init(out->result, 1, 1, par->Fp);
    
    // F(x) = aux'_R
    fq_mat_set(out->result, aux_R_prime, par->Fp);
}

// ========== Protocol 9  ==========

/**
 */
void scheme3_compute_lambda_prime(fq_t *lambda_prime, fq_t *v_values, int d, pubpar *par) {
    fq_t *v = (fq_t*)malloc((d+1) * sizeof(fq_t));
    for(int i = 0; i <= d; i++) {
        fq_init(v[i], par->Fp);
        if (i == 0) {
            fq_one(v[i], par->Fp);
        } else {
            fq_set(v[i], v_values[i-1], par->Fp);
        }
    }
    
    for(int j = 0; j <= d; j++) {
        fq_t numerator, denominator, temp, second_part;
        fq_init(numerator, par->Fp);
        fq_init(denominator, par->Fp);
        fq_init(temp, par->Fp);
        fq_init(second_part, par->Fp);
        
        fq_one(numerator, par->Fp);
        fq_one(denominator, par->Fp);
        fq_one(second_part, par->Fp);
        
        for(int i = 0; i <= d; i++) {
            if (i != j) {
                fq_mul(numerator, numerator, v[i], par->Fp);
            }
        }
        
        for(int i = 0; i <= d; i++) {
            for(int k = i + 1; k <= d; k++) {
                fq_sub(temp, v[k], v[i], par->Fp);
                fq_mul(second_part, second_part, temp, par->Fp);
            }
        }
        
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
        
        fq_clear(numerator, par->Fp);
        fq_clear(denominator, par->Fp);
        fq_clear(temp, par->Fp);
        fq_clear(second_part, par->Fp);
    }
    
    for(int i = 0; i <= d; i++) {
        fq_clear(v[i], par->Fp);
    }
    free(v);
}

// ==========  ==========

void paillier_pk_clear(paillier_pk_t *pk) {
    fmpz_clear(pk->n);
    fmpz_clear(pk->n2);
    fmpz_clear(pk->g);
    fmpz_clear(pk->h);
    fmpz_clear(pk->f);
    fmpz_clear(pk->mu);
}

void paillier_ct_clear(paillier_ct_t *ct) {
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
            fmpz_clear(ct->c[i][j]);
        }
    }
}

void hss_ek_clear(hss_ek_t *ek) {
    fmpz_clear(ek->ek);
}

void hss_mv_clear(hss_mv_t *mv) {
    fmpz_clear(mv->mv[0]);
    fmpz_clear(mv->mv[1]);
}

void scheme3_keygen_out_clear(scheme3_keygen_out *out) {
    paillier_pk_clear(&out->pk);
    hss_ek_clear(&out->ek[0]);
    hss_ek_clear(&out->ek[1]);
    fmpz_clear(out->sk);
}

void scheme3_probgen_out_clear(scheme3_probgen_out *out, pubpar *par) {
    int d = par->d;
    
    if(out->vk_x) {
        for(int j = 0; j <= d; j++) {
            if(out->vk_x[j]) {
                free(out->vk_x[j]);
            }
        }
        free(out->vk_x);
    }
    
    if(out->ct_inputs) {
        int m = par->m;
        for(int j = 0; j < (d+1)*m; j++) {
            paillier_ct_clear(&out->ct_inputs[j]);
        }
        free(out->ct_inputs);
    }
    
    if(out->v_values) {
        for(int j = 0; j < d; j++) {
            fq_clear(out->v_values[j], par->Fp);
        }
        free(out->v_values);
    }
}

void scheme3_compute_out_clear(scheme3_compute_out *out, pubpar *par) {
    int d = par->d;
    
    if(out->mv) {
        for(int j = 0; j <= d; j++) {
            hss_mv_clear(&out->mv[j]);
        }
        free(out->mv);
    }
}

void scheme3_verify_out_clear(scheme3_verify_out *out, pubpar *par) {
    fq_mat_clear(out->aux_R_prime, par->Fp);
}

void scheme3_reconstruct_out_clear(scheme3_reconstruct_out *out, pubpar *par) {
    fq_mat_clear(out->result, par->Fp);
}
