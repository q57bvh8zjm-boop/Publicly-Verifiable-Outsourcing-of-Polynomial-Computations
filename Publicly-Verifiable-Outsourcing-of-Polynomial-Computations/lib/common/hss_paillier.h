#ifndef PV_HSS_PAILLIER_H
#define PV_HSS_PAILLIER_H
#include "crypto_util.h"

/* OSY21, Construction 4.4: Paillier HSS with encrypted secret-key digits.
 * Values are nonnegative integer lifts; each server reduces its signed output
 * modulo the application field only at the final Output step.
 */
typedef struct {
    fmpz_t n, n2, field, value_bound;
    fmpz_t *encrypted_digits;
    unsigned modulus_bits, base_bits, digit_count;
    int max_m, max_d;
    uint64_t conversion_budget;
} hss_public_key;
typedef struct {
    fmpz_t *digits;
    unsigned count;
    unsigned char prf_key[32];
} hss_evaluation_key;
typedef struct { fmpz_t *c; unsigned count; } hss_ciphertext;

int hss_setup(hss_public_key *pk, hss_evaluation_key ek[2], const pubpar *par);
void hss_public_key_clear(hss_public_key *pk);
void hss_evaluation_key_clear(hss_evaluation_key *ek);
void hss_input(hss_ciphertext *out, const hss_public_key *pk, const fq_t input,
                const fq_ctx_t field);
void hss_ciphertext_clear(hss_ciphertext *ct);
int hss_eval(fq_t out, int server, const hss_public_key *pk,
             const hss_evaluation_key *ek, hss_ciphertext *inputs,
             fq_mat_t polynomial, pubpar *par);
size_t hss_input_bytes(const hss_public_key *pk);
#endif
