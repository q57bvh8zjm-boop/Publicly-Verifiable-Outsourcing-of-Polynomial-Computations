# Publicly Verifiable Outsourcing of Polynomial Computations

This artifact implements three publicly verifiable polynomial-outsourcing schemes and a private classification application. It includes the protocol source, correctness and tampering checks, reproducible workloads, and recorded measurements.

## Build and run

Install a C compiler, FLINT, GMP, libsodium, and Python 3. On Ubuntu/Debian:

```sh
sudo apt-get install build-essential libflint-dev libgmp-dev libsodium-dev python3
cd Publicly-Verifiable-Outsourcing-of-Polynomial-Computations
make check
```

On macOS, install `flint`, `gmp`, and `libsodium` with Homebrew, then run the same Make targets. Start with [REVIEWER_GUIDE.md](REVIEWER_GUIDE.md) for the validation sequence.

## Private classification

A client delegates a fixed public polynomial SVM score, verifies and decodes the result, and applies the decision threshold locally. The Spambase experiment uses 57 features and a fixed split of 3,220 training and 1,381 test emails. At `t=1`, quadratic and cubic classifiers use three and four independent servers.

Use Python 3.9–3.12 with the pinned application dependencies:

```sh
python3 -m pip install -r applications/private-classification/requirements.txt
make classification-check
```

The command checks 48 protocol queries, rejects 96 single-response perturbations, and verifies that quantization preserves all 1,381 floating-point predictions for each model. Results go to `results/measured/local-classification/`. Use `CLASSIFICATION_OUTPUT=results/measured/local-classification-2` for a new run.

The recorded quadratic and cubic accuracies are 94.50% and 94.93%. [Application instructions](applications/private-classification/README.md) cover model regeneration, encoding, the measured costs in Appendix G/Table 7, and the protocol-feature comparison in Table 6. The public dataset, model fixtures, original measurements, and per-query predictions are included.

## Protocols

| Manuscript name | Code / CLI name | Source |
|---|---|---|
| $\Gamma_{\mathrm{FE}}^{\mathrm{Sh}}$ | `Gamma_FE_Sh` | `protocols/gamma_fe_sh.c` |
| $\Gamma_{\mathrm{OV}}^{\mathrm{Sh}}$ | `Gamma_OV_Sh` | `protocols/gamma_ov_sh.c` |
| $\Gamma_{\mathrm{OV}}^{\mathrm{HSS}}$ | `Gamma_OV_HSS` | `protocols/gamma_ov_hss.c`, `lib/common/hss_paillier.c` |
| $\Pi_3$ | `Pi3` | `protocols/protocol2.c` |

The Shamir constructions use `k=d*t+1` servers. The HSS construction uses two servers at threshold one. Its integer-lifting model, modulus bounds, and security and adaptive-correctness conditions are specified in [HSS_INSTANTIATION.md](docs/HSS_INSTANTIATION.md). The `Pi3` source supplies the historical polynomial-evaluator regression.

Create public parameters with `pv_params_init` and release successful outputs with their matching cleanup functions. Public verification consumes the verification key and response values; client decoding uses the client's secrets. Verification keys come from the trusted query generator. The application keeps complete response transcripts inside the client/auditor trust domain and gives each server its own input and response. The benchmark APIs run in one process; network transport is outside their timing boundary.

```sh
make test-hss
make test-hss-high
```

These suites generate fresh large Paillier keys and can take several minutes.

## Polynomial measurements

```sh
make benchmark
./build/benchmark Gamma_OV_Sh 200 2 1 10 dense
python3 scripts/run_experiments.py --profile smoke --trials 3 --output results/measured/local-smoke
```

Arguments are `SCHEME M D T TRIALS dense`. `SCHEME` accepts `Gamma_FE_Sh`, `Gamma_OV_Sh`, and `Gamma_OV_HSS`, with numeric aliases `1`, `2`, and `3`. Every timing row follows successful public verification and equality with native evaluation.

The runner saves CSVs, commands, source hashes, library versions, and machine metadata. Public workloads are deterministic; cryptographic masks, scaling points, keys, and encryption coins use operating-system randomness. [MEASUREMENTS.md](docs/MEASUREMENTS.md) defines the phases and full-grid commands. Earlier plot coordinates remain documented in [results/reported/](results/reported/README.md); new executions have their own raw measurements.

## Polynomial representation and names

A polynomial of total degree at most `d` in `m` variables has `C(m+d,d)` coefficients. Monomials are ordered by increasing degree, then lexicographically by their nondecreasing variable-index tuples:

```text
1; x0, x1; x0^2, x0*x1, x1^2; x0^3, x0^2*x1, x0*x1^2, x1^3; ...
```

Verification uses the prime-order Ristretto255 scalar field. The extension representation is described in [README_extended_field.md](lib/common/README_extended_field.md).

FE denotes field extensions, OV orthogonal vectors, Sh Shamir sharing, and HSS homomorphic secret sharing. The former `Gamma1`, `Gamma2`, and `Gamma3` labels map to the three canonical names above. Source/API prefixes are `gamma_fe_sh`, `gamma_ov_sh`, and `gamma_ov_hss`. Historical raw records retain their recorded names; plotting translates them to manuscript notation.
