# Implemented HSS backend

`lib/common/hss_paillier.c` implements the digit-decomposed Paillier design of Orlandi, Scholl, and Yakoubov, *The Rise of Paillier: Homomorphic Secret Sharing and Public-Key Silent OT*, EUROCRYPT 2021, Construction 4.4. [Primary source](https://iacr.org/archive/eurocrypt2021/126960235/126960235.pdf); [ePrint record](https://eprint.iacr.org/2021/262).

This replaces the earlier Paillier–ElGamal-style prototype. It is a different concrete backend: its encrypted secret-key digits require a Paillier KDM assumption, as well as a PRF assumption. The current implementation therefore must not be described as an instantiation under DCR alone. A trusted setup generates fresh safe primes, the public modulus, encrypted key digits, and two evaluation keys. The factorization and combined decryption key are discarded after setup. The two evaluation keys must be distributed separately.

## Arithmetic invariant

Let `q` be the Ristretto scalar order, `N` the Paillier modulus, `B=2^base_bits`, and `s=sum_i s_i B^i` a decryption exponent with `s=0 mod phi(N)` and `s=1 mod N`. The public key contains independent encryptions of the digits `s_i`. An input ciphertext contains encryptions of `x` and `x*s_i`, with fresh rerandomization of every component.

A memory value represents subtractive integer shares of `(y, y*s_0, ..., y*s_(L-1))`. Multiplication reconstructs the exponent share `sum_i memory[i+1]*B^i`, raises each input component to that exponent, and applies the quotient/remainder distributed-log conversion. A common keyed mask indexed by the ciphertext/function digest, gate, and component is added modulo `N`. The two local evaluations use the same mask and are deterministic for a fixed key and query.

Each monomial starts from shares of one and multiplies by every variable in its tuple. Coefficients are canonical nonnegative lifts from `F_q`. The evaluator accumulates coefficient-weighted monomial shares over the integers. It does not reduce the accumulator modulo `N`. Server zero negates its accumulator and both servers reduce modulo `q`; their field outputs add to the polynomial value whenever the intermediate integer conversions are correct. Gamma3 subtracts `F(0)` on server one before the OV check and adds it back after verification.

## Implemented range calculation

The following bound is derived for this evaluator's nonnegative integer program. Let

```text
C = binomial(m+d,d)
V = C*q^(d+1)
Q = 16*C*d*(d+1)
margin = 128 + ceil(log2(Q)) + 4
```

`V` bounds each monomial, coefficient-weighted result, and accumulated unmasked value. A digit-bearing intermediate is below `B*V`. `Q` bounds the number of component conversions across all `d+1` OV columns; the implementation checks that the number of components is at most 16.

Starting at `n=3072`, setup increases `n` in steps of 256 until

```text
n-2 >= bitlength(V) + floor(n/4) + margin.
```

Both primes have `n/2` bits, hence `N >= 2^(n-2)`. The resulting parameter choice gives `Q*B*V/N < 2^-132`. This is a union bound for wrap failures with fresh ideal uniform conversion masks in a fixed honest computation. The digit-share interval is separately sized so its total statistical shift distance is at most `2^-130`. These terms describe two different effects; PRF security and primality-generation error are separate from these arithmetic bounds.

The application field is never replaced with a small test field to obtain favorable timings. A degree-10 computation over the full field needs a larger modulus than the degree-2 example. `paillier_bits` in the timing files records the requested key size; the actual product can have one fewer bit. Programs beyond the setup's variable/degree limits, wrong fields, malformed matrix dimensions, and nonunit ciphertext components are rejected.

## Scope of the security claim

The parameter calculation and correctness tests establish the implemented arithmetic model and detect concrete functional regressions. They do not prove the manuscript's stronger correctness condition for inputs selected after an evaluation key is disclosed. That key contains the PRF seed, so replacing its outputs by secret ideal masks in such an adaptive game requires a separate argument. Hashing the query into the mask domain prevents accidental reuse; it does not supply that argument.

Accordingly, Gamma3 remains an experimental conditional instantiation. The manuscript must either prove this stronger property for this exact setup, input formation, and program class, or retain it as an explicit assumption. A claim covering multiple evaluations must account for the total conversion budget and the chosen adaptive model. Passing high-degree tests cannot close this proof obligation. GMP/FLINT secret-key operations are not implemented as constant-time operations.
