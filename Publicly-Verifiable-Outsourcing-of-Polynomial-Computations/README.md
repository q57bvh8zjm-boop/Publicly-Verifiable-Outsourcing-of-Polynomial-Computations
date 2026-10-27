# Publicly Verifiable Outsourcing of Polynomial Computations

This research artifact implements three polynomial-outsourcing prototypes with a shared dense monomial order. The repaired implementation adds complete high-degree HSS evaluation, full-field randomness, public verification interfaces, end-to-end regression tests, and executable timing experiments.

For the paper's independent-VM experiments, start with [VM_EXPERIMENTS.md](docs/VM_EXPERIMENTS.md). It covers TCP deployment, Figures 2–5, the private Pi3 comparison, and the primitive microbenchmarks. [REVIEWER_GUIDE.md](REVIEWER_GUIDE.md) describes the correctness tests. The two-server construction's assumptions and bounds are specified in [HSS_INSTANTIATION.md](docs/HSS_INSTANTIATION.md).

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

| Manuscript name | Code / CLI name | Source | Current scope |
|---|---|---|---|
| $\Gamma_{\mathrm{FE}}^{\mathrm{Sh}}$ | `Gamma_FE_Sh` | `protocols/gamma_fe_sh.c` | Extension-field shares, encoded public verification key, client decoding |
| $\Gamma_{\mathrm{OV}}^{\mathrm{Sh}}$ | `Gamma_OV_Sh` | `protocols/gamma_ov_sh.c` | Shamir shares over the full Ristretto scalar field, OV verification |
| $\Gamma_{\mathrm{OV}}^{\mathrm{HSS}}$ | `Gamma_OV_HSS` | `protocols/gamma_ov_hss.c`, `lib/common/hss_paillier.c` | Two-server digit-decomposed Paillier HSS, all polynomial degrees, OV verification |
| $\Pi_3$ | `Pi3` | `protocols/pi3.c` | Private extension-field query, server computation, verification, and decoding |

Gamma_FE_Sh, Gamma_OV_Sh, and Pi3 use `k=d*t+1` servers. Gamma_OV_HSS uses two servers and a threshold of one. Public verification inputs and client secrets have separate interfaces. Create parameters with `pv_params_init` and clear successful outputs with the matching cleanup functions. The TCP harness in `benchmarks/vm/` serializes each server's inputs and runs the complete protocol across independent processes.

## Independent-VM measurements

```sh
make experiments
python3 scripts/run_vm_experiments.py --config experiments/vms.json --profile smoke --trials 3 --output results/measured/vm-smoke
python3 scripts/run_vm_experiments.py --config experiments/vms.json --profile paper --trials 30 --output results/measured/vm-paper
python3 scripts/run_microbenchmarks.py --output results/measured/vm-micro
python3 scripts/plot_vm_results.py results/measured/vm-paper --output results/plots/vm-paper
```

Copy `experiments/vms.example.json` to `experiments/vms.json`, set each server VM's address, and start its worker using the commands in the [VM guide](docs/VM_EXPERIMENTS.md). Every sample includes each server's computation time. The paper's `T_s` is their maximum; `T_c` is ProbGen + Verify + Decode. TCP wall times and setup are recorded separately. Figure 4 compares FE with Pi3; Figure 5 uses degrees 4, 8, and 10 with base-10 logarithms.

## Local kernel measurements

```sh
make benchmark
./build/benchmark Gamma_OV_Sh 200 2 1 10 dense
./build/benchmark Gamma_OV_Sh 4 2 1 10 var
python3 scripts/run_experiments.py --profile smoke --trials 3 --output results/measured/local-smoke
```

The benchmark arguments are `SCHEME M D T TRIALS dense|var`, with `SCHEME` set to `Gamma_FE_Sh`, `Gamma_OV_Sh`, or `Gamma_OV_HSS`. Legacy numeric selectors `1`, `2`, and `3` remain accepted; new output uses the canonical names. Every timing row follows successful verification and equality with native evaluation. The legacy `var` selector runs a vector-autoregression polynomial example; the manuscript's Verifiable Aggregate Receipts application has a separate implementation plan.

The runner stores each command's CSV and stderr, a combined `raw.csv`, compiler/library information, source hashes, and machine metadata. Public workload fixtures are deterministic; cryptographic masks, scaling points, keys, and encryption coins use operating-system randomness. See [MEASUREMENTS.md](docs/MEASUREMENTS.md) for the timing definitions and the full-grid command.

## Naming correspondence

FE denotes field extensions, OV orthogonal vectors, Sh Shamir sharing, and HSS homomorphic secret sharing. The former Gamma1, Gamma2, and Gamma3 names map to `Gamma_FE_Sh`, `Gamma_OV_Sh`, and `Gamma_OV_HSS`, respectively. Source files and C API prefixes are `gamma_fe_sh`, `gamma_ov_sh`, and `gamma_ov_hss`. Include the corresponding header and use that prefix when updating code written against the former `scheme1`, `scheme2`, or `scheme3` APIs.

Existing raw CSVs, execution logs, commands, and source hashes retain their original labels as recorded. Plotting translates those labels to the manuscript notation.

## Historical data

`results/reported/` contains archived coordinates transcribed from an earlier paper version. They are not measurements of the repaired implementation. Their provenance and known disagreement with the supplied manuscript are recorded in that directory. `make reported-plots` redraws them with an archival label. Use `results/measured/` for current runs.

## Polynomial representation

A polynomial of total degree at most `d` in `m` variables has `C(m+d,d)` coefficients. Index zero is the constant. Subsequent monomials are ordered by increasing degree, then lexicographically by their nondecreasing variable-index tuples:

```text
1; x0, x1; x0^2, x0*x1, x1^2; x0^3, x0^2*x1, x0*x1^2, x1^3; ...
```

Every current evaluator, including HSS, uses this order. All arithmetic exposed to verification uses the prime-order Ristretto255 scalar field. The extension-field representation is described in [lib/common/README_extended_field.md](lib/common/README_extended_field.md).
