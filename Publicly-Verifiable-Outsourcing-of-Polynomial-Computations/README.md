# Publicly Verifiable Outsourcing of Polynomial Computations

This research artifact implements three polynomial-outsourcing prototypes with a shared dense monomial order. The repaired implementation adds complete high-degree HSS evaluation, full-field randomness, public verification interfaces, end-to-end regression tests, and executable timing experiments.

Start with [REVIEWER_GUIDE.md](REVIEWER_GUIDE.md). The [revision notes](docs/REVISION_NOTES.md) map the changes to the review findings. The two-server construction uses an experimental Paillier HSS backend with explicit assumptions and bounds in [HSS_INSTANTIATION.md](docs/HSS_INSTANTIATION.md).

## Build and test

Ubuntu/Debian:

```sh
sudo apt-get install build-essential libflint-dev libgmp-dev libsodium-dev python3
cd Publicly-Verifiable-Outsourcing-of-Polynomial-Computations
make check
make test-hss
make test-hss-high
```

On macOS, install `flint`, `gmp`, and `libsodium` with Homebrew and run the same Make targets from this directory. `make check` runs the fast suite. The separate HSS suites generate fresh large keys and can take several minutes.

| Construction | Source | Current scope |
|---|---|---|
| Gamma1 | `protocols/scheme1.c` | Extension-field shares, encoded public verification key, client decoding |
| Gamma2 | `protocols/scheme2.c` | Shamir shares over the full Ristretto scalar field, OV verification |
| Gamma3 | `protocols/scheme3.c`, `lib/common/hss_paillier.c` | Two-server digit-decomposed Paillier HSS, all polynomial degrees, OV verification |
| Pi3 | `protocols/protocol2.c` | Historical baseline; polynomial evaluator regression only |

Gamma1 and Gamma2 use `k=d*t+1` servers. Gamma3 uses two servers and a threshold of one. The protocol APIs changed to separate public verification inputs from client secrets; the tests and benchmark are complete usage examples. Create public parameters with `pv_params_init` and clear every successful output with its matching cleanup function. Verification keys must come from the trusted query generator. These are in-process research APIs; they do not implement a network transport or wire-message parser.

## Run measurements

```sh
make benchmark
./build/benchmark 2 200 2 1 10 dense
./build/benchmark 2 4 2 1 10 var
python3 scripts/run_experiments.py --profile smoke --trials 3 --output results/measured/local-smoke
```

The benchmark arguments are `SCHEME M D T TRIALS dense|var`. Every timing row is emitted only after successful verification and equality with native evaluation. The VAR workload computes the quadratic score of a one-lag vector autoregression forecast and also checks a direct matrix-vector reference. It measures this polynomial kernel; signatures, receipts, binding proofs, transport, and application preprocessing are outside the measurement.

The runner stores each command's CSV and stderr, a combined `raw.csv`, compiler/library information, source hashes, and machine metadata. Public workload fixtures are deterministic; cryptographic masks, scaling points, keys, and encryption coins use operating-system randomness. See [MEASUREMENTS.md](docs/MEASUREMENTS.md) for the timing definitions and the full-grid command.

## Historical data

`results/reported/` contains archived coordinates transcribed from an earlier paper version. They are not measurements of the repaired implementation. Their provenance and known disagreement with the supplied manuscript are recorded in that directory. `make reported-plots` redraws them with an archival label. Use `results/measured/` for current runs.

## Polynomial representation

A polynomial of total degree at most `d` in `m` variables has `C(m+d,d)` coefficients. Index zero is the constant. Subsequent monomials are ordered by increasing degree, then lexicographically by their nondecreasing variable-index tuples:

```text
1; x0, x1; x0^2, x0*x1, x1^2; x0^3, x0^2*x1, x0*x1^2, x1^3; ...
```

Every current evaluator, including HSS, uses this order. All arithmetic exposed to verification uses the prime-order Ristretto255 scalar field. The extension-field representation is described in [lib/common/README_extended_field.md](lib/common/README_extended_field.md).
