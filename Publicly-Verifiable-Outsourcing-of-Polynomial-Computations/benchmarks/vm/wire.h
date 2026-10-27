#ifndef PV_VM_WIRE_H
#define PV_VM_WIRE_H
#include "gamma_ov_hss.h"
#include "gamma_ov_sh.h"
#include "pi3.h"
#include <stddef.h>
#include <stdint.h>
#define VM_WIRE_VERSION 1
#define VM_FRAME_LIMIT ((size_t)512 * 1024 * 1024)
#define VM_MAX_SERVERS 128
#define VM_MAX_COEFFICIENTS_DEFAULT 5000000
/* Worker protocol; controller protocol uses the batch and row message types. */
enum {
    VM_SETUP = 1,
    VM_READY = 2,
    VM_QUERY = 3,
    VM_LOADED = 4,
    VM_GO = 5,
    VM_RESULT = 6,
    VM_BYE = 7,
    VM_BATCH_SETUP = 20,
    VM_BATCH_READY = 21,
    VM_BATCH_QUERY = 22,
    VM_BATCH_RESULT = 23,
    VM_ROW = 24,
    VM_DONE = 25
};
typedef struct {
    unsigned char *data;
    size_t len, cap, pos;
} wire_buf;
void vm_fail(const char *why);
void vm_check(int ok, const char *why);
uint64_t vm_ns(void);
int vm_number(const char *s, int min, int max);
int vm_scheme(const char *s);
const char *vm_scheme_name(int s);
void wire_clear(wire_buf *b);
void wire_put(wire_buf *b, const void *p, size_t n);
void wire_get(wire_buf *b, void *p, size_t n);
void wire_u32(wire_buf *b, uint32_t v);
uint32_t wire_get_u32(wire_buf *b);
void wire_u64(wire_buf *b, uint64_t v);
uint64_t wire_get_u64(wire_buf *b);
void wire_end(wire_buf *b);
void wire_scalar(wire_buf *b, const fmpz_t v, const fmpz_t q);
void wire_get_scalar(wire_buf *b, fmpz_t v, const fmpz_t q);
void wire_fq(wire_buf *b, const fq_t v, const pubpar *p);
void wire_get_fq(wire_buf *b, fq_t v, const pubpar *p);
void wire_ext(wire_buf *b, const ext_field_elem_t *v, const pubpar *p);
void wire_get_ext(wire_buf *b, ext_field_elem_t *v, const pubpar *p);
void wire_integer(wire_buf *b, const fmpz_t v);
void wire_get_integer(wire_buf *b, fmpz_t v);
void wire_frame_write(int fd, uint32_t type, const wire_buf *b);
int wire_frame_read(int fd, uint32_t *type, wire_buf *b);
void vm_expect_frame(int fd, uint32_t type, wire_buf *b);
void vm_fixture(fq_t out, const char *kind, int m, int d, slong index, int trial, pubpar *p);
void vm_fill_polynomial(fq_mat_t F, pubpar *p);
slong vm_coefficient_count(int m, int d, uint64_t max_count);
void wire_hss_setup(wire_buf *b, const gamma_ov_hss_keygen_out *key, int server);
void wire_get_hss_setup(wire_buf *b, hss_public_key *pk, hss_evaluation_key *ek, pubpar *p);
void wire_hss_inputs(wire_buf *b, const hss_ciphertext *inputs, const hss_public_key *pk,
                     pubpar *p);
hss_ciphertext *wire_get_hss_inputs(wire_buf *b, const hss_public_key *pk, pubpar *p);
#endif
