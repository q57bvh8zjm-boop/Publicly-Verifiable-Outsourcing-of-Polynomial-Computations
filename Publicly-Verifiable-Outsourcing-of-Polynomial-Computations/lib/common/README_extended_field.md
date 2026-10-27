# Quadratic extension arithmetic

An `ext_field_elem_t` stores `a+b*w` over the prime field, with `w^2=-c`. A valid field requires `-c` to be a quadratic nonresidue. Gamma1 chooses such a `c` for the Ristretto scalar field during query generation.

Multiplication computes `(a1*a2-c*b1*b2, a1*b2+a2*b1)` into temporaries before writing the output, so either input may alias the output. Random elements use independent uniform full-field coordinates from `crypto_util.c`.

Initialize contexts/elements/matrices before use and call their corresponding clear functions. The end-to-end tests in `tests/test_protocols.c` and the benchmark demonstrate supported use. Historical standalone demo filenames and CMake targets are no longer part of this artifact.
