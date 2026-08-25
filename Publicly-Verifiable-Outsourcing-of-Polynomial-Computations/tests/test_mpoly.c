#include <stdio.h>
#include <stdlib.h>

#include "../lib/common/mpoly.h"
#include "../lib/common/extended_field.h"
#include "../protocols/scheme1.h"
#include "../protocols/protocol2.h"

static void fail(const char *message) {
    fprintf(stderr, "FAIL: %s\n", message);
    exit(EXIT_FAILURE);
}

static void fill_dense_example(fq_mat_t F, fq_mat_t X, const fq_ctx_t ctx) {
    for (slong i = 0; i < 10; i++)
        fq_set_ui(fq_mat_entry(F, i, 0), (ulong)(i + 1), ctx);
    fq_set_ui(fq_mat_entry(X, 0, 0), 2, ctx);
    fq_set_ui(fq_mat_entry(X, 1, 0), 3, ctx);
}

static void test_coefficient_counts(void) {
    if (mpoly_num_coefficients(2, 3) != 10)
        fail("C(2 + 3, 3) must be 10");
    if (mpoly_num_coefficients(200, 2) != 20301)
        fail("Figure 2 coefficient count is incorrect");
    if (mpoly_num_coefficients(15, 10) != 3268760)
        fail("Figure 4 coefficient count is incorrect");
    if (mpoly_num_coefficients(0, 3) != 1)
        fail("a zero-variable polynomial must contain only its constant");
    if (mpoly_num_coefficients(-1, 2) != -1)
        fail("invalid dimensions must be rejected");
}

static void test_base_field_eval(void) {
    pubpar par = {0};
    par.m = 2;
    par.d = 3;
    fmpz_init(par.p);
    fmpz_set_ui(par.p, 101);
    fq_ctx_init(par.Fp, par.p, 1, "a");

    fq_mat_t F, X;
    fq_mat_init(F, 10, 1, par.Fp);
    fq_mat_init(X, 2, 1, par.Fp);
    fill_dense_example(F, X, par.Fp);

    fq_t result, expected;
    Eval(result, F, X, &par);
    fq_init(expected, par.Fp);
    fq_set_ui(expected, 92, par.Fp);
    if (!fq_equal(result, expected, par.Fp))
        fail("base-field degree-3 evaluation returned the wrong value");

    fq_clear(expected, par.Fp);
    fq_clear(result, par.Fp);
    fq_mat_clear(X, par.Fp);
    fq_mat_clear(F, par.Fp);
    fq_ctx_clear(par.Fp);
    fmpz_clear(par.p);
}

static void test_extension_field_eval(void) {
    pubpar par = {0};
    par.m = 2;
    par.d = 3;
    fmpz_init(par.p);
    fmpz_set_ui(par.p, 101);
    fq_ctx_init(par.Fp, par.p, 1, "a");

    fq_mat_t F, base_X;
    fq_mat_init(F, 10, 1, par.Fp);
    fq_mat_init(base_X, 2, 1, par.Fp);
    fill_dense_example(F, base_X, par.Fp);

    fmpz_t modulus_constant, value, zero;
    fmpz_init(modulus_constant);
    fmpz_init(value);
    fmpz_init(zero);
    fmpz_set_ui(modulus_constant, 1);
    fmpz_zero(zero);
    ext_field_ctx_t ext_ctx;
    ext_field_ctx_init(&ext_ctx, par.p, modulus_constant);

    ext_field_mat_t X;
    ext_field_mat_init(&X, 2, 1);
    for (int i = 0; i < 2; i++) {
        fmpz_set_ui(value, (ulong)(i + 2));
        ext_field_elem_set(ext_field_mat_entry(&X, i, 0), value, zero, &ext_ctx);
    }

    ext_field_elem_t result;
    ext_field_elem_init(&result);
    scheme1_eval_extended(&result, F, &X, &ext_ctx, &par);
    if (fmpz_cmp_ui(result.a, 92) != 0 || !fmpz_is_zero(result.b))
        fail("extension-field degree-3 evaluation returned the wrong value");

    ext_field_elem_clear(&result);
    ext_field_mat_clear(&X);
    ext_field_ctx_clear(&ext_ctx);
    fmpz_clear(zero);
    fmpz_clear(value);
    fmpz_clear(modulus_constant);
    fq_mat_clear(base_X, par.Fp);
    fq_mat_clear(F, par.Fp);
    fq_ctx_clear(par.Fp);
    fmpz_clear(par.p);
}

static void test_protocol2_eval(void) {
    pubpar_p2 par = {0};
    par.m = 2;
    par.d = 3;
    fmpz_init(par.p);
    fmpz_set_ui(par.p, 101);
    fq_ctx_init(par.Fq2, par.p, 2, "z");

    fq_mat_t F, X;
    fq_mat_init(F, 10, 1, par.Fq2);
    fq_mat_init(X, 2, 1, par.Fq2);
    fill_dense_example(F, X, par.Fq2);

    fq_t result, expected;
    protocol2_eval(result, F, X, &par);
    fq_init(expected, par.Fq2);
    fq_set_ui(expected, 92, par.Fq2);
    if (!fq_equal(result, expected, par.Fq2))
        fail("Pi_3 degree-3 evaluation returned the wrong value");

    fq_clear(expected, par.Fq2);
    fq_clear(result, par.Fq2);
    fq_mat_clear(X, par.Fq2);
    fq_mat_clear(F, par.Fq2);
    fq_ctx_clear(par.Fq2);
    fmpz_clear(par.p);
}

int main(void) {
    test_coefficient_counts();
    test_base_field_eval();
    test_extension_field_eval();
    test_protocol2_eval();
    puts("All polynomial evaluation tests passed.");
    return EXIT_SUCCESS;
}
