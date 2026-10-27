#define _POSIX_C_SOURCE 200809L
#include "vm/wire.h"

/* Each independent trial samples fresh full-field operands outside the timer.
 * Repetition amortizes timer overhead, without subtracting an estimated cost. */
int main(int argc, char **argv) {
    if (argc != 4) {
        fputs("Usage: microbenchmark TRIALS INNER_REPEATS INCLUDE_PAILLIER(0|1)\n", stderr);
        return 2;
    }
    int trials = vm_number(argv[1], 1, 100000);
    int inner = vm_number(argv[2], 1, 1000000);
    int paillier = vm_number(argv[3], 0, 1);
    pubpar p;
    vm_check(pv_params_init(&p, 1, 1, 1), "parameters");
    fq_t a, b, c;
    fq_init(a, p.Fp);
    fq_init(b, p.Fp);
    fq_init(c, p.Fp);
    ext_field_ctx_t field;
    ext_field_ctx_init_default(&field, p.p);
    ext_field_elem_t ea, eb, ec;
    ext_field_elem_init(&ea);
    ext_field_elem_init(&eb);
    ext_field_elem_init(&ec);
    hss_public_key pk;
    hss_evaluation_key ek[2];
    uint64_t keygen_ns = 0;
    if (paillier) {
        uint64_t begin = vm_ns();
        vm_check(hss_setup(&pk, ek, &p), "Paillier key generation");
        keygen_ns = vm_ns() - begin;
    }
    fprintf(stderr, "FLINT=%s GMP=%s sodium=%s; Paillier setup_ms=%.6f\n", FLINT_VERSION,
            gmp_version, sodium_version_string(), keygen_ns / 1e6);
    puts("operation,trial,inner_runs,ns_per_operation,modulus_bits,setup_ms");
    const char *names[] = {"Add_q", "Mul_q", "Add_q2", "Mul_q2", "Exp_G", "Paillier_encrypt"};
    fmpz_t message, ciphertext;
    fmpz_init(message);
    fmpz_init(ciphertext);
    unsigned char scalar[32], point[32], result[32];
    for (int op = 0; op < (paillier ? 6 : 5); op++) {
        for (int trial = -1; trial < trials; trial++) {
            pv_random_field(a, p.Fp);
            pv_random_field(b, p.Fp);
            ext_field_elem_random(&ea, &field);
            ext_field_elem_random(&eb, &field);
            do {
                crypto_core_ristretto255_scalar_random(scalar);
            } while (sodium_is_zero(scalar, sizeof(scalar)));
            crypto_scalarmult_ristretto255_base(point, scalar);
            do {
                crypto_core_ristretto255_scalar_random(scalar);
            } while (sodium_is_zero(scalar, sizeof(scalar)));
            fq2fmpz(message, a, p.Fp);
            int repetitions = op == 5 ? 1 : inner;
            uint64_t begin = vm_ns();
            int valid = 1;
            for (int i = 0; i < repetitions; i++) {
                switch (op) {
                case 0:
                    fq_add(c, a, b, p.Fp);
                    break;
                case 1:
                    fq_mul(c, a, b, p.Fp);
                    break;
                case 2:
                    ext_field_elem_add(&ec, &ea, &eb, &field);
                    break;
                case 3:
                    ext_field_elem_mul(&ec, &ea, &eb, &field);
                    break;
                case 4:
                    valid &= crypto_scalarmult_ristretto255(result, scalar, point) == 0;
                    break;
                case 5:
                    hss_paillier_encrypt(ciphertext, message, &pk);
                    break;
                }
            }
            double elapsed = (double)(vm_ns() - begin) / repetitions;
            vm_check(valid, "group multiplication");
            if (trial >= 0)
                printf("%s,%d,%d,%.9f,%u,%.6f\n", names[op], trial, repetitions, elapsed,
                       op == 5 ? pk.modulus_bits : 0, keygen_ns / 1e6);
        }
    }
    fmpz_clear(message);
    fmpz_clear(ciphertext);
    if (paillier) {
        hss_evaluation_key_clear(&ek[0]);
        hss_evaluation_key_clear(&ek[1]);
        hss_public_key_clear(&pk);
    }
    fq_clear(a, p.Fp);
    fq_clear(b, p.Fp);
    fq_clear(c, p.Fp);
    ext_field_elem_clear(&ea);
    ext_field_elem_clear(&eb);
    ext_field_elem_clear(&ec);
    ext_field_ctx_clear(&field);
    pv_params_clear(&p);
    return 0;
}
