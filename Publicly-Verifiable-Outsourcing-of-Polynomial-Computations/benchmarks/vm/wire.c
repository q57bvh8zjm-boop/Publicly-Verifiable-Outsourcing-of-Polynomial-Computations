#define _POSIX_C_SOURCE 200809L
#include "wire.h"
#include <errno.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>
void vm_fail(const char *s) {
    fprintf(stderr, "experiment error: %s\n", s);
    exit(1);
}
void vm_check(int ok, const char *s) {
    if (!ok)
        vm_fail(s);
}
uint64_t vm_ns(void) {
    struct timespec t;
    vm_check(!clock_gettime(CLOCK_MONOTONIC, &t), "clock_gettime");
    return (uint64_t)t.tv_sec * 1000000000ULL + t.tv_nsec;
}
int vm_number(const char *s, int min, int max) {
    char *end;
    errno = 0;
    long n = strtol(s, &end, 10);
    vm_check(!errno && *s && !*end && n >= min && n <= max, "invalid integer argument");
    return (int)n;
}
const char *vm_scheme_name(int s) {
    static const char *names[] = {"", "Gamma_FE_Sh", "Gamma_OV_Sh", "Gamma_OV_HSS", "Pi3"};
    vm_check(s >= 1 && s <= 4, "invalid scheme");
    return names[s];
}
int vm_scheme(const char *s) {
    for (int i = 1; i <= 4; i++)
        if (!strcmp(s, vm_scheme_name(i)))
            return i;
    vm_fail("unknown scheme");
    return 0;
}
void wire_clear(wire_buf *b) {
    if (b->data) {
        sodium_memzero(b->data, b->len);
        free(b->data);
    }
    memset(b, 0, sizeof(*b));
}
void wire_put(wire_buf *b, const void *p, size_t n) {
    vm_check(n <= VM_FRAME_LIMIT && b->len <= VM_FRAME_LIMIT - n, "frame limit exceeded");
    size_t need = b->len + n;
    if (need > b->cap) {
        size_t cap = b->cap ? b->cap : 256;
        while (cap < need)
            cap = cap > VM_FRAME_LIMIT / 2 ? VM_FRAME_LIMIT : cap * 2;
        unsigned char *v = realloc(b->data, cap);
        vm_check(v != NULL, "allocation");
        b->data = v;
        b->cap = cap;
    }
    if (n)
        memcpy(b->data + b->len, p, n);
    b->len = need;
}
void wire_get(wire_buf *b, void *p, size_t n) {
    vm_check(b->pos <= b->len && n <= b->len - b->pos, "truncated payload");
    if (n)
        memcpy(p, b->data + b->pos, n);
    b->pos += n;
}
void wire_u32(wire_buf *b, uint32_t v) {
    unsigned char a[4];
    for (int i = 0; i < 4; i++)
        a[i] = (unsigned char)(v >> (24 - 8 * i));
    wire_put(b, a, 4);
}
uint32_t wire_get_u32(wire_buf *b) {
    unsigned char a[4];
    wire_get(b, a, 4);
    uint32_t v = 0;
    for (int i = 0; i < 4; i++)
        v = (v << 8) | a[i];
    return v;
}
void wire_u64(wire_buf *b, uint64_t v) {
    unsigned char a[8];
    for (int i = 0; i < 8; i++)
        a[i] = (unsigned char)(v >> (56 - 8 * i));
    wire_put(b, a, 8);
}
uint64_t wire_get_u64(wire_buf *b) {
    unsigned char a[8];
    wire_get(b, a, 8);
    uint64_t v = 0;
    for (int i = 0; i < 8; i++)
        v = (v << 8) | a[i];
    return v;
}
void wire_end(wire_buf *b) { vm_check(b->pos == b->len, "trailing payload bytes"); }
static void put_mpz(wire_buf *b, const fmpz_t v, size_t n) {
    unsigned char *bytes = pv_alloc(n, 1);
    mpz_t z;
    mpz_init(z);
    fmpz_get_mpz(z, v);
    vm_check(mpz_sgn(z) >= 0 && (mpz_sizeinbase(z, 2) + 7) / 8 <= n, "integer encoding overflow");
    mpz_export(bytes, NULL, -1, 1, 0, 0, z);
    wire_put(b, bytes, n);
    sodium_memzero(bytes, n);
    free(bytes);
    mpz_clear(z);
}
static void get_mpz(wire_buf *b, fmpz_t v, size_t n) {
    vm_check(b->pos <= b->len && n <= b->len - b->pos, "truncated integer");
    mpz_t z;
    mpz_init(z);
    mpz_import(z, n, -1, 1, 0, 0, b->data + b->pos);
    b->pos += n;
    fmpz_set_mpz(v, z);
    mpz_clear(z);
}
void wire_scalar(wire_buf *b, const fmpz_t v, const fmpz_t q) {
    vm_check(fmpz_sgn(v) >= 0 && fmpz_cmp(v, q) < 0, "noncanonical scalar");
    put_mpz(b, v, 32);
}
void wire_get_scalar(wire_buf *b, fmpz_t v, const fmpz_t q) {
    get_mpz(b, v, 32);
    vm_check(fmpz_cmp(v, q) < 0, "scalar outside field");
}
void wire_fq(wire_buf *b, const fq_t v, const pubpar *p) {
    fmpz_t a;
    fmpz_init(a);
    fq2fmpz(a, v, p->Fp);
    wire_scalar(b, a, p->p);
    fmpz_clear(a);
}
void wire_get_fq(wire_buf *b, fq_t v, const pubpar *p) {
    fmpz_t a;
    fmpz_init(a);
    wire_get_scalar(b, a, p->p);
    fq_set_fmpz(v, a, p->Fp);
    fmpz_clear(a);
}
void wire_ext(wire_buf *b, const ext_field_elem_t *v, const pubpar *p) {
    wire_scalar(b, v->a, p->p);
    wire_scalar(b, v->b, p->p);
}
void wire_get_ext(wire_buf *b, ext_field_elem_t *v, const pubpar *p) {
    wire_get_scalar(b, v->a, p->p);
    wire_get_scalar(b, v->b, p->p);
}
void wire_integer(wire_buf *b, const fmpz_t v) {
    vm_check(fmpz_sgn(v) >= 0, "negative wire integer");
    size_t n = (fmpz_bits(v) + 7) / 8;
    if (!n)
        n = 1;
    vm_check(n <= 4096, "oversized integer");
    wire_u32(b, (uint32_t)n);
    put_mpz(b, v, n);
}
void wire_get_integer(wire_buf *b, fmpz_t v) {
    size_t n = wire_get_u32(b);
    vm_check(n >= 1 && n <= 4096, "invalid integer length");
    get_mpz(b, v, n);
    vm_check(n == 1 || fmpz_bits(v) > 8 * (n - 1), "noncanonical integer");
}
static int read_all(int fd, unsigned char *p, size_t n, int allow_eof) {
    size_t done = 0;
    while (done < n) {
        ssize_t r = read(fd, p + done, n - done);
        if (r < 0 && errno == EINTR)
            continue;
        if (!r && allow_eof && !done)
            return 0;
        vm_check(r > 0, "connection closed or read failed");
        done += (size_t)r;
    }
    return 1;
}
void wire_frame_write(int fd, uint32_t type, const wire_buf *b) {
    wire_buf h = {0};
    wire_u32(&h, type);
    wire_u64(&h, b ? b->len : 0);
    size_t parts[2] = {h.len, b ? b->len : 0};
    const unsigned char *ptr[2] = {h.data, b ? b->data : NULL};
    for (int j = 0; j < 2; j++) {
        size_t done = 0;
        while (done < parts[j]) {
            ssize_t n = write(fd, ptr[j] + done, parts[j] - done);
            if (n < 0 && errno == EINTR)
                continue;
            vm_check(n > 0, "write failed");
            done += (size_t)n;
        }
    }
    wire_clear(&h);
}
int wire_frame_read(int fd, uint32_t *type, wire_buf *b) {
    unsigned char header[12];
    if (!read_all(fd, header, 12, 1))
        return 0;
    wire_buf h = {header, 12, 12, 0};
    *type = wire_get_u32(&h);
    uint64_t n = wire_get_u64(&h);
    vm_check(n <= VM_FRAME_LIMIT, "frame too large");
    wire_clear(b);
    b->data = pv_alloc((size_t)n, 1);
    b->len = b->cap = (size_t)n;
    read_all(fd, b->data, b->len, 0);
    return 1;
}
void vm_expect_frame(int fd, uint32_t type, wire_buf *b) {
    uint32_t got;
    vm_check(wire_frame_read(fd, &got, b) && got == type, "unexpected message type");
}
void vm_fixture(fq_t out, const char *kind, int m, int d, slong index, int trial, pubpar *p) {
    char label[128];
    int n = snprintf(label, sizeof(label), "pv-benchmark-v1:%s:%d:%d:%ld:%d", kind, m, d,
                     (long)index, trial);
    unsigned char digest[64];
    crypto_hash_sha512(digest, (unsigned char *)label, n);
    mpz_t z;
    mpz_init(z);
    mpz_import(z, 64, -1, 1, 0, 0, digest);
    fmpz_t a;
    fmpz_init(a);
    fmpz_set_mpz(a, z);
    fq_set_fmpz(out, a, p->Fp);
    fmpz_clear(a);
    mpz_clear(z);
}
void vm_fill_polynomial(fq_mat_t F, pubpar *p) {
    for (slong i = 0; i < fq_mat_nrows(F, p->Fp); i++)
        vm_fixture(fq_mat_entry(F, i, 0), "coefficient", p->m, p->d, i, 0, p);
}
slong vm_coefficient_count(int m, int d, uint64_t max_count) {
    slong n = mpoly_num_coefficients(m, d);
    vm_check(n > 0 && (uint64_t)n <= max_count, "dense polynomial exceeds coefficient limit; run "
                                                "preflight and size RAM before raising the limit");
    return n;
}
void wire_hss_setup(wire_buf *b, const gamma_ov_hss_keygen_out *k, int server) {
    const hss_public_key *p = &k->pk;
    const hss_evaluation_key *e = &k->ek[server];
    wire_u32(b, p->modulus_bits);
    wire_u32(b, p->base_bits);
    wire_u32(b, p->digit_count);
    wire_u64(b, p->conversion_budget);
    wire_integer(b, p->n);
    wire_integer(b, p->n2);
    wire_integer(b, p->field);
    wire_integer(b, p->value_bound);
    for (unsigned i = 0; i < p->digit_count; i++)
        wire_integer(b, p->encrypted_digits[i]);
    for (unsigned i = 0; i < e->count; i++)
        wire_integer(b, e->digits[i]);
    wire_put(b, e->prf_key, 32);
}
void wire_get_hss_setup(wire_buf *b, hss_public_key *p, hss_evaluation_key *e, pubpar *par) {
    memset(p, 0, sizeof(*p));
    memset(e, 0, sizeof(*e));
    p->modulus_bits = wire_get_u32(b);
    p->base_bits = wire_get_u32(b);
    p->digit_count = wire_get_u32(b);
    p->conversion_budget = wire_get_u64(b);
    vm_check(p->modulus_bits >= 3072 && p->modulus_bits <= 16384 && p->modulus_bits % 256 == 0 &&
                 p->base_bits == p->modulus_bits / 4 && p->digit_count >= 1 && p->digit_count <= 15,
             "invalid HSS dimensions");
    p->max_m = par->m;
    p->max_d = par->d;
    fmpz_init(p->n);
    fmpz_init(p->n2);
    fmpz_init(p->field);
    fmpz_init(p->value_bound);
    wire_get_integer(b, p->n);
    wire_get_integer(b, p->n2);
    wire_get_integer(b, p->field);
    wire_get_integer(b, p->value_bound);
    vm_check(fmpz_bits(p->n) >= p->modulus_bits - 1 && fmpz_bits(p->n) <= p->modulus_bits &&
                 fmpz_is_odd(p->n) && fmpz_equal(p->field, par->p),
             "invalid HSS modulus or field");
    fmpz_t expected;
    fmpz_init(expected);
    fmpz_mul(expected, p->n, p->n);
    vm_check(fmpz_equal(expected, p->n2), "HSS modulus square mismatch");
    fmpz_pow_ui(expected, par->p, (ulong)par->d + 1);
    fmpz_mul_ui(expected, expected, (ulong)mpoly_num_coefficients(par->m, par->d));
    vm_check(fmpz_equal(expected, p->value_bound), "HSS value bound mismatch");
    uint64_t budget = 16 * (uint64_t)mpoly_num_coefficients(par->m, par->d) *
                      (uint64_t)par->d * ((uint64_t)par->d + 1);
    vm_check(p->conversion_budget == budget, "HSS conversion budget mismatch");
    fmpz_clear(expected);
    p->encrypted_digits = pv_alloc(p->digit_count, sizeof(fmpz_t));
    e->count = p->digit_count;
    e->digits = pv_alloc(e->count, sizeof(fmpz_t));
    for (unsigned i = 0; i < p->digit_count; i++) {
        fmpz_init(p->encrypted_digits[i]);
        wire_get_integer(b, p->encrypted_digits[i]);
        vm_check(fmpz_sgn(p->encrypted_digits[i]) > 0 &&
                     fmpz_cmp(p->encrypted_digits[i], p->n2) < 0,
                 "invalid encrypted key digit");
    }
    for (unsigned i = 0; i < e->count; i++) {
        fmpz_init(e->digits[i]);
        wire_get_integer(b, e->digits[i]);
        vm_check(fmpz_bits(e->digits[i]) <= p->base_bits + PV_STAT_BITS + 8,
                 "oversized evaluation-key share");
    }
    wire_get(b, e->prf_key, 32);
}
void wire_hss_inputs(wire_buf *b, const hss_ciphertext *ct, const hss_public_key *pk, pubpar *p) {
    for (size_t i = 0; i < (size_t)p->m * (p->d + 1); i++) {
        vm_check(ct[i].count == pk->digit_count + 1, "ciphertext component count");
        for (unsigned j = 0; j < ct[i].count; j++)
            wire_integer(b, ct[i].c[j]);
    }
}
hss_ciphertext *wire_get_hss_inputs(wire_buf *b, const hss_public_key *pk, pubpar *p) {
    size_t count = (size_t)p->m * (p->d + 1);
    vm_check(count <= b->len / 5, "impossible ciphertext count");
    hss_ciphertext *ct = pv_alloc(count, sizeof(*ct));
    for (size_t i = 0; i < count; i++) {
        ct[i].count = pk->digit_count + 1;
        ct[i].c = pv_alloc(ct[i].count, sizeof(fmpz_t));
        for (unsigned j = 0; j < ct[i].count; j++) {
            fmpz_init(ct[i].c[j]);
            wire_get_integer(b, ct[i].c[j]);
            vm_check(fmpz_sgn(ct[i].c[j]) > 0 && fmpz_cmp(ct[i].c[j], pk->n2) < 0,
                     "ciphertext outside modulus");
        }
    }
    return ct;
}
