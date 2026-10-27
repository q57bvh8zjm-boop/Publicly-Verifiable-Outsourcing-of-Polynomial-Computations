# Reviewer quick start

Run commands from `Publicly-Verifiable-Outsourcing-of-Polynomial-Computations/` after installing FLINT, GMP, libsodium, a C compiler, and Python 3.

```sh
make check            # build, archive-grid checks, fast end-to-end regressions
make test-hss         # x^2, dense mixed terms, full-field inputs, constants, tampering
make test-hss-high    # HSS powers 4, 6, 8, 10 at q-1 with a range-sized modulus
```

The fast suite runs Gamma_FE_Sh, Gamma_OV_Sh, and Pi3 at every degree 1–10 and threshold 1–3. It checks native equality, fresh query shares, changed outputs, malformed verification points, zero polynomials, and constants. The HSS suite also checks deterministic repeated evaluation and rejects malformed ciphertext components.

Check TCP execution and then deploy using [the VM guide](docs/VM_EXPERIMENTS.md):

```sh
make experiments
python3 tests/test_vm.py --hss
```

The VM runner covers Figures 2–5, records the maximum server time and all individual times, and supports resumed runs. The microbenchmark runner records the six primitive costs.

Generate new measurements:

```sh
python3 scripts/run_experiments.py --profile smoke --trials 3 --output results/measured/local-review
```

Inspect `raw.csv`, individual trial files, and `metadata.json`. The benchmark records setup, query generation, server work, verification, reconstruction, and native evaluation separately. It stops on a failed correctness check.

The [revision notes](docs/REVISION_NOTES.md) describe the changes. The [HSS instantiation note](docs/HSS_INSTANTIATION.md) gives the integer-lifting model, parameter bound, setup trust, and security assumptions. The [measurement note](docs/MEASUREMENTS.md) distinguishes VM execution from the local kernel harness.

`make reported-plots` only redraws historical coordinates. Those files do not establish performance of the repaired implementation. The former PIR estimators were removed from the active tree because they used the superseded API and estimated costs instead of measuring the complete protocol. Their source remains in Git history.
