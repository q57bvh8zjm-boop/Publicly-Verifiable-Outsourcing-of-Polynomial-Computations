#include "hss_paillier.h"
#include <limits.h>

static void safe_prime(fmpz_t prime, unsigned bits) {
    mpz_t half, candidate;
    mpz_init(half); mpz_init(candidate);
    fmpz_t random;
    fmpz_init(random);
    do {
        pv_random_bits(random, bits - 1);
        fmpz_get_mpz(half, random);
        mpz_setbit(half, bits - 2); mpz_setbit(half, 0);
        mpz_nextprime(half, half);
        mpz_mul_2exp(candidate, half, 1); mpz_add_ui(candidate, candidate, 1);
    } while (mpz_sizeinbase(candidate, 2) != bits ||
             mpz_probab_prime_p(candidate, 64) == 0);
    fmpz_set_mpz(prime, candidate);
    fmpz_clear(random); mpz_clear(half); mpz_clear(candidate);
}

static unsigned ceil_log2_u64(uint64_t n) {
    unsigned bits = 0;
    if (n) n--;
    while (n) { n >>= 1; bits++; }
    return bits;
}

static void encrypt(fmpz_t out, const fmpz_t message, const hss_public_key *pk) {
    fmpz_t r, gcd, lift;
    fmpz_init(r); fmpz_init(gcd); fmpz_init(lift);
    do {
        pv_random_below(r, pk->n);
        fmpz_gcd(gcd, r, pk->n);
    } while (!fmpz_is_one(gcd));
    fmpz_powm(out, r, pk->n, pk->n2);
    fmpz_mul(lift, message, pk->n); fmpz_add_ui(lift, lift, 1);
    fmpz_mul(out, out, lift); fmpz_mod(out, out, pk->n2);
    fmpz_clear(r); fmpz_clear(gcd); fmpz_clear(lift);
}

int hss_setup(hss_public_key *pk, hss_evaluation_key ek[2], const pubpar *par) {
    memset(pk, 0, sizeof(*pk)); memset(ek, 0, 2 * sizeof(*ek));
    slong count = mpoly_num_coefficients(par->m, par->d);
    if (par->m < 1 || par->d < 1 || count <= 0 ||
        (uint64_t)count > UINT64_MAX / (uint64_t)par->d / ((uint64_t)par->d + 1) / 128)
        return 0;
    fmpz_init(pk->n); fmpz_init(pk->n2); fmpz_init(pk->field); fmpz_init(pk->value_bound);
    fmpz_set(pk->field, par->p);
    pk->max_m = par->m; pk->max_d = par->d;
    /* The bound covers every input and every coefficient in F_q. */
    fmpz_pow_ui(pk->value_bound, par->p, (ulong)par->d + 1);
    fmpz_mul_ui(pk->value_bound, pk->value_bound, (ulong)count);
    uint64_t gates = (uint64_t)count * (uint64_t)par->d * ((uint64_t)par->d + 1);
    /* L <= 9 with base_bits = modulus_bits / 4. Reserve 16 components. */
    pk->conversion_budget = 16 * gates;
    unsigned range_bits = (unsigned)fmpz_bits(pk->value_bound);
    unsigned margin = PV_STAT_BITS + ceil_log2_u64(pk->conversion_budget) + 4;
    unsigned bits = 3072;
    while (bits - 2 < range_bits + bits / 4 + margin) bits += 256;
    pk->modulus_bits = bits; pk->base_bits = bits / 4;

    fmpz_t p, q, phi, pm1, qm1, secret, base, digit, remainder, share_bound;
    fmpz_init(p); fmpz_init(q); fmpz_init(phi); fmpz_init(pm1); fmpz_init(qm1);
    fmpz_init(secret); fmpz_init(base); fmpz_init(digit); fmpz_init(remainder);
    fmpz_init(share_bound);
    safe_prime(p, bits / 2);
    do { safe_prime(q, bits / 2); } while (fmpz_equal(p, q));
    fmpz_mul(pk->n, p, q); fmpz_mul(pk->n2, pk->n, pk->n);
    fmpz_sub_ui(pm1, p, 1); fmpz_sub_ui(qm1, q, 1); fmpz_mul(phi, pm1, qm1);
    /* secret = 0 mod phi(N), 1 mod N, so C^secret = (1+N)^plaintext. */
    if (!fmpz_invmod(secret, phi, pk->n)) abort();
    fmpz_mul(secret, secret, phi);
    pk->digit_count = (unsigned)((fmpz_bits(secret) + pk->base_bits - 1) / pk->base_bits);
    if (pk->digit_count + 1 > 16) abort();
    pk->encrypted_digits = pv_alloc(pk->digit_count, sizeof(fmpz_t));
    for (int server = 0; server < 2; server++) {
        ek[server].count = pk->digit_count;
        ek[server].digits = pv_alloc(pk->digit_count, sizeof(fmpz_t));
    }
    randombytes_buf(ek[0].prf_key, 32); memcpy(ek[1].prf_key, ek[0].prf_key, 32);
    fmpz_one(base); fmpz_mul_2exp(base, base, pk->base_bits);
    /* Statistical hiding of all digit shares, including either disclosed key. */
    fmpz_mul_2exp(share_bound, base, PV_STAT_BITS + ceil_log2_u64(pk->digit_count) + 2);
    fmpz_set(remainder, secret);
    for (unsigned i = 0; i < pk->digit_count; i++) {
        fmpz_init(pk->encrypted_digits[i]);
        fmpz_init(ek[0].digits[i]); fmpz_init(ek[1].digits[i]);
        fmpz_mod(digit, remainder, base);
        fmpz_fdiv_q_2exp(remainder, remainder, pk->base_bits);
        encrypt(pk->encrypted_digits[i], digit, pk);
        pv_random_below(ek[0].digits[i], share_bound);
        fmpz_add(ek[1].digits[i], ek[0].digits[i], digit);
    }
    fmpz_clear(p); fmpz_clear(q); fmpz_clear(phi); fmpz_clear(pm1); fmpz_clear(qm1);
    fmpz_clear(secret); fmpz_clear(base); fmpz_clear(digit); fmpz_clear(remainder);
    fmpz_clear(share_bound);
    return 1;
}

void hss_input(hss_ciphertext *out, const hss_public_key *pk, const fq_t input,
                const fq_ctx_t field) {
    out->count = pk->digit_count + 1;
    out->c = pv_alloc(out->count, sizeof(fmpz_t));
    for (unsigned i = 0; i < out->count; i++) fmpz_init(out->c[i]);
    fmpz_t value, zero, rerandomization;
    fmpz_init(value); fmpz_init(zero); fmpz_init(rerandomization);
    fq2fmpz(value, input, field);
    encrypt(out->c[0], value, pk);
    for (unsigned i = 0; i < pk->digit_count; i++) {
        encrypt(rerandomization, zero, pk);
        fmpz_powm(out->c[i+1], pk->encrypted_digits[i], value, pk->n2);
        fmpz_mul(out->c[i+1], out->c[i+1], rerandomization);
        fmpz_mod(out->c[i+1], out->c[i+1], pk->n2);
    }
    fmpz_clear(value); fmpz_clear(zero); fmpz_clear(rerandomization);
}

static void hash_integer(crypto_generichash_state *state, const fmpz_t value) {
    mpz_t z; mpz_init(z); fmpz_get_mpz(z, value);
    size_t n = (mpz_sizeinbase(z, 2) + 7) / 8, actual = 0;
    unsigned char *bytes = pv_alloc(n, 1), length[8];
    mpz_export(bytes, &actual, -1, 1, 0, 0, z);
    for (int i = 0; i < 8; i++) length[i] = (unsigned char)((uint64_t)actual >> (8*i));
    crypto_generichash_update(state, length, 8);
    crypto_generichash_update(state, bytes, actual);
    free(bytes); mpz_clear(z);
}

static void evaluation_context(unsigned char out[32], hss_ciphertext *inputs,
                                fq_mat_t F, const hss_public_key *pk, pubpar *par) {
    crypto_generichash_state hash;
    crypto_generichash_init(&hash, NULL, 0, 32);
    const unsigned char domain[] = "pv-hss-paillier-eval-v1";
    crypto_generichash_update(&hash, domain, sizeof(domain));
    hash_integer(&hash, pk->n); hash_integer(&hash, par->p);
    fmpz_t value; fmpz_init(value);
    fmpz_set_ui(value, (ulong)par->m); hash_integer(&hash, value);
    fmpz_set_ui(value, (ulong)par->d); hash_integer(&hash, value);
    fmpz_set_ui(value, (ulong)fq_mat_nrows(F, par->Fp)); hash_integer(&hash, value);
    for (slong i = 0; i < fq_mat_nrows(F, par->Fp); i++) {
        fq2fmpz(value, fq_mat_entry(F, i, 0), par->Fp); hash_integer(&hash, value);
    }
    for (int i = 0; i < par->m; i++)
        for (unsigned j = 0; j < inputs[i].count; j++) hash_integer(&hash, inputs[i].c[j]);
    fmpz_clear(value); crypto_generichash_final(&hash, out, 32);
}

/* Deterministic, keyed rejection sampling; the ciphertext/function digest and
 * gate/component identifiers prevent mask reuse across evaluations.
 */
static void mask_for(fmpz_t out, const hss_public_key *pk, const hss_evaluation_key *ek,
                     const unsigned char context[32], uint64_t gate, unsigned component) {
    unsigned char message[52], seed[32]; memcpy(message, context, 32);
    for (int i = 0; i < 8; i++) message[32+i] = (unsigned char)(gate >> (8*i));
    for (int i = 0; i < 4; i++) message[40+i] = (unsigned char)(component >> (8*i));
    unsigned bits = (unsigned)fmpz_bits(pk->n);
    size_t length = (bits + 7U) / 8U;
    unsigned char *bytes = pv_alloc(length, 1);
    mpz_t value; mpz_init(value);
    for (uint64_t retry = 0;; retry++) {
        for (int i = 0; i < 8; i++) message[44+i] = (unsigned char)(retry >> (8*i));
        crypto_generichash(seed, 32, message, sizeof(message), ek->prf_key, 32);
        randombytes_buf_deterministic(bytes, length, seed);
        if (bits % 8) bytes[length-1] &= (1U << (bits % 8)) - 1;
        mpz_import(value, length, -1, 1, 0, 0, bytes);
        fmpz_set_mpz(out, value);
        if (fmpz_cmp(out, pk->n) < 0) break;
    }
    sodium_memzero(seed, sizeof(seed)); sodium_memzero(bytes, length);
    free(bytes); mpz_clear(value);
}

static void multiply(fmpz_t *out, fmpz_t *memory, const hss_ciphertext *input,
                      const hss_public_key *pk, const hss_evaluation_key *ek,
                      const unsigned char context[32], uint64_t gate) {
    fmpz_t exponent, term, quotient, remainder, inverse, mask;
    fmpz_init(exponent); fmpz_init(term); fmpz_init(quotient);
    fmpz_init(remainder); fmpz_init(inverse); fmpz_init(mask);
    for (unsigned i = pk->digit_count; i > 0; i--) {
        fmpz_mul_2exp(exponent, exponent, pk->base_bits);
        fmpz_add(exponent, exponent, memory[i]);
    }
    for (unsigned j = 0; j <= pk->digit_count; j++) {
        fmpz_powm(term, input->c[j], exponent, pk->n2);
        fmpz_fdiv_qr(quotient, remainder, term, pk->n);
        if (!fmpz_invmod(inverse, remainder, pk->n)) abort();
        fmpz_mul(out[j], quotient, inverse);
        mask_for(mask, pk, ek, context, gate, j);
        fmpz_add(out[j], out[j], mask); fmpz_mod(out[j], out[j], pk->n);
    }
    fmpz_clear(exponent); fmpz_clear(term); fmpz_clear(quotient);
    fmpz_clear(remainder); fmpz_clear(inverse); fmpz_clear(mask);
}

typedef struct {
    const hss_public_key *pk;
    const hss_evaluation_key *ek;
    hss_ciphertext *inputs;
    fq_mat_struct *F;
    pubpar *par;
    int server, *variables;
    slong coefficient, count;
    uint64_t gate;
    unsigned char context[32];
    fmpz_t sum;
} eval_state;

static void evaluate_monomial(eval_state *s, int degree) {
    fmpz_t coefficient, term;
    fmpz_init(coefficient); fmpz_init(term);
    fq2fmpz(coefficient, fq_mat_entry(s->F, s->coefficient++, 0), s->par->Fp);
    if (fmpz_is_zero(coefficient)) { fmpz_clear(coefficient); fmpz_clear(term); return; }
    unsigned length = s->pk->digit_count + 1;
    fmpz_t *memory = pv_alloc(length, sizeof(fmpz_t));
    fmpz_t *next = pv_alloc(length, sizeof(fmpz_t));
    for (unsigned j = 0; j < length; j++) { fmpz_init(memory[j]); fmpz_init(next[j]); }
    fmpz_set_ui(memory[0], (ulong)s->server);
    for (unsigned j = 1; j < length; j++) fmpz_set(memory[j], s->ek->digits[j-1]);
    for (int i = 0; i < degree; i++) {
        multiply(next, memory, &s->inputs[s->variables[i]], s->pk, s->ek, s->context, s->gate++);
        fmpz_t *swap = memory; memory = next; next = swap;
    }
    fmpz_mul(term, coefficient, memory[0]);
    /* Preserve integer differences. Reducing each share modulo N here would
     * destroy the final conversion to F_q when a share wraps around N.
     */
    fmpz_add(s->sum, s->sum, term);
    for (unsigned j = 0; j < length; j++) { fmpz_clear(memory[j]); fmpz_clear(next[j]); }
    free(memory); free(next); fmpz_clear(coefficient); fmpz_clear(term);
}

static void visit_degree(eval_state *s, int start, int depth, int degree) {
    if (s->coefficient >= s->count) return;
    if (depth == degree) { evaluate_monomial(s, degree); return; }
    for (int variable = start; variable < s->par->m && s->coefficient < s->count; variable++) {
        s->variables[depth] = variable;
        visit_degree(s, variable, depth + 1, degree);
    }
}

int hss_eval(fq_t out, int server, const hss_public_key *pk,
             const hss_evaluation_key *ek, hss_ciphertext *inputs,
             fq_mat_t polynomial, pubpar *par) {
    slong count = mpoly_num_coefficients(par->m, par->d);
    if (server < 0 || server > 1 || par->m > pk->max_m || par->d > pk->max_d ||
        count <= 0 || fq_mat_nrows(polynomial, par->Fp) != count ||
        fq_mat_ncols(polynomial, par->Fp) != 1 || ek->count != pk->digit_count ||
        !fmpz_equal(par->p, pk->field)) return 0;
    fmpz_t gcd; fmpz_init(gcd);
    for (int i = 0; i < par->m; i++) {
        if (inputs[i].count != pk->digit_count + 1) { fmpz_clear(gcd); return 0; }
        for (unsigned j = 0; j < inputs[i].count; j++) {
            if (fmpz_sgn(inputs[i].c[j]) <= 0 || fmpz_cmp(inputs[i].c[j], pk->n2) >= 0) {
                fmpz_clear(gcd); return 0;
            }
            fmpz_gcd(gcd, inputs[i].c[j], pk->n);
            if (!fmpz_is_one(gcd)) { fmpz_clear(gcd); return 0; }
        }
    }
    fmpz_clear(gcd);
    eval_state state = {.pk=pk, .ek=ek, .inputs=inputs, .F=polynomial, .par=par,
                        .server=server, .coefficient=0, .count=count, .gate=0};
    state.variables = pv_alloc((size_t)par->d, sizeof(int));
    fmpz_init(state.sum);
    evaluation_context(state.context, inputs, polynomial, pk, par);
    for (int degree = 0; degree <= par->d; degree++) visit_degree(&state, 0, 0, degree);
    if (!server) fmpz_neg(state.sum, state.sum);
    fmpz_mod(state.sum, state.sum, par->p); fq_set_fmpz(out, state.sum, par->Fp);
    free(state.variables); fmpz_clear(state.sum);
    return 1;
}

size_t hss_input_bytes(const hss_public_key *pk) {
    return ((size_t)pk->digit_count + 1) * ((fmpz_bits(pk->n2) + 7) / 8);
}
void hss_ciphertext_clear(hss_ciphertext *ct) {
    for (unsigned i = 0; i < ct->count; i++) fmpz_clear(ct->c[i]);
    free(ct->c); memset(ct, 0, sizeof(*ct));
}
void hss_evaluation_key_clear(hss_evaluation_key *ek) {
    for (unsigned i = 0; i < ek->count; i++) fmpz_clear(ek->digits[i]);
    free(ek->digits); sodium_memzero(ek, sizeof(*ek));
}
void hss_public_key_clear(hss_public_key *pk) {
    for (unsigned i = 0; i < pk->digit_count; i++) fmpz_clear(pk->encrypted_digits[i]);
    free(pk->encrypted_digits);
    fmpz_clear(pk->n); fmpz_clear(pk->n2); fmpz_clear(pk->field); fmpz_clear(pk->value_bound);
    memset(pk, 0, sizeof(*pk));
}
