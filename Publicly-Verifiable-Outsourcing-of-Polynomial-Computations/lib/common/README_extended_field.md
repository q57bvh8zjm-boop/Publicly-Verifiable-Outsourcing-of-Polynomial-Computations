# 

## 

 $\mathbb{F}_p$  $\mathbb{F}_{p^2}$ 。 $a + b\omega$， $a, b \in \mathbb{F}_p$， $\omega$  $g(\omega) = \omega^2 + c$。

## 

- `extended_field.h` - ，
- `extended_field.c` - ，
- `extended_field_demo.c` - ，
- `README_extended_field.md` - 

## 

###  (`ext_field_elem_t`)

```c
typedef struct {
    fmpz_t a;              // 
    fmpz_t b;              // ω
} ext_field_elem_t;
```

 $a + b\omega$。

###  (`ext_field_ctx_t`)

```c
typedef struct {
    fmpz_t p;              // F_p
    fmpz_t c;              // g(ω) = ω^2 + c
    flint_rand_t state;    // 
} ext_field_ctx_t;
```

。

## 

### 
- `ext_field_ctx_init()` - 
- `ext_field_ctx_clear()` - 

### 
- `ext_field_elem_init()` - 
- `ext_field_elem_clear()` - 
- `ext_field_elem_zero()` - 
- `ext_field_elem_one()` - 
- `ext_field_elem_set()` - 
- `ext_field_elem_copy()` - 
- `ext_field_elem_random()` - 
- `ext_field_elem_print()` - 

### 
- `ext_field_elem_add()` - 
- `ext_field_elem_mul()` - 
- `ext_field_elem_inv()` - 
- `ext_field_elem_equal()` - 

## 

### 
$(a_1 + b_1\omega) + (a_2 + b_2\omega) = (a_1 + a_2) + (b_1 + b_2)\omega$

### 
$(a_1 + b_1\omega) \times (a_2 + b_2\omega) = (a_1a_2 - b_1b_2c) + (a_1b_2 + a_2b_1)\omega$

 $\omega^2 = -c$。

### 
 $\alpha = a + b\omega$，：
$$\alpha^{-1} = \frac{a - b\omega}{a^2 + b^2c}$$

## 

```c
#include "extended_field.h"

int main() {
    //  F_97^2
    ext_field_ctx_t ctx;
    fmpz_t p, c;
    fmpz_init(p);
    fmpz_init(c);
    fmpz_set_ui(p, 97);
    fmpz_set_ui(c, 2);
    ext_field_ctx_init(&ctx, p, c);
    
    // 
    ext_field_elem_t a, b, result;
    ext_field_elem_init(&a);
    ext_field_elem_init(&b);
    ext_field_elem_init(&result);
    
    // 
    fmpz_t val1, val2;
    fmpz_init(val1);
    fmpz_init(val2);
    
    fmpz_set_ui(val1, 3);
    fmpz_set_ui(val2, 5);
    ext_field_elem_set(&a, val1, val2, &ctx);  // a = 3 + 5ω
    
    fmpz_set_ui(val1, 7);
    fmpz_set_ui(val2, 2);
    ext_field_elem_set(&b, val1, val2, &ctx);  // b = 7 + 2ω
    
    // 
    ext_field_elem_add(&result, &a, &b, &ctx);
    
    // 
    ext_field_elem_clear(&a);
    ext_field_elem_clear(&b);
    ext_field_elem_clear(&result);
    fmpz_clear(val1);
    fmpz_clear(val2);
    ext_field_ctx_clear(&ctx);
    fmpz_clear(p);
    fmpz_clear(c);
    
    return 0;
}
```

## 

### 
```bash
mkdir build && cd build
CC=gcc cmake ..
make extended_field
```

### 
```bash
make extended_field_test
./test/correctness/extended_field_test
```

### 
```bash
gcc -I../lib/common -I/usr/include/flint -o extended_field_demo ../lib/common/extended_field_demo.c libextended_field.a -lflint -lgmp
./extended_field_demo
```




## 

- FLINT ：
- GMP ：
