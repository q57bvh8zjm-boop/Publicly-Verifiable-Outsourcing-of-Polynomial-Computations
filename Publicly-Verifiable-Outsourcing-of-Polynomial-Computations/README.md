# Publicly Verifiable Outsourcing of Polynomial Computations

For a one-page protocol map and the exact commands for generating Figures 2--5, start with **[`REVIEWER_GUIDE.md`](REVIEWER_GUIDE.md)**.

This repository contains the C artifact accompanying the paper *Publicly Verifiable Outsourcing of Polynomial Computations*. It implements the three proposed constructions, Gamma1, Gamma2, and Gamma3, and includes the Pi3 private-verification baseline used in the paper.

The repository also records the exact relationship between implementation files, experiment parameters, and paper figures. Numerical coordinates and plotting tools are included for Figures 2--5 and the reported tables.

## Quick start

On Ubuntu 22.04 or a compatible Debian-based system:

```bash
sudo apt-get update
sudo apt-get install build-essential libflint-dev libgmp-dev libsodium-dev python3
make check
```

`make check` validates the reported-data grids, checks source hygiene, builds all shipped programs, and runs the polynomial-evaluation tests.

To redraw the coordinates reported in the paper:

```bash
python3 -m pip install matplotlib
make reported-plots
```

The plots are written to `results/plots/`.

## Repository layout

```text
.
|-- protocols/
|   |-- scheme1.c, scheme1.h          Gamma1
|   |-- scheme2.c, scheme2.h          Gamma2
|   |-- scheme3.c, scheme3.h          Gamma3
|   |-- protocol2.c, protocol2.h      Pi3 comparison baseline
|   |-- scheme3_pir_optimize.c        exploratory PIR estimator
|   `-- protocol2_pir_optimize.c      exploratory baseline estimator
|-- lib/common/
|   |-- mpoly.c, mpoly.h              dense polynomial evaluation
|   `-- extended_field.c, .h          explicit F_(q^2) arithmetic
|-- tests/test_mpoly.c                 degree-3 evaluator regression tests
|-- results/reported/                  coordinates transcribed from the paper
|-- scripts/
|   |-- plot_reported_results.py       reported-coordinate plotting
|   `-- validate_artifact.py           grid and source validation
|-- REVIEWER_GUIDE.md                  one-page reviewer quick start
`-- Makefile
```

## Construction-to-source mapping

| Paper construction | Implementation | Main implementation path |
|---|---|---|
| Gamma1 | Publicly verifiable MSVC over an explicit extension field | `protocols/scheme1.c`, `lib/common/extended_field.c` |
| Gamma2 | Shamir-shared curves with orthogonal-vector verification | `protocols/scheme2.c`, `lib/common/mpoly.c` |
| Gamma3 | Two-server Paillier/HSS construction | `protocols/scheme3.c` |
| Pi3 baseline | Private-verification extension-field construction | `protocols/protocol2.c` |

## Paper experiments and included data

The paper defines the client improvement ratio

```text
Rt* = (Tp + Tv + Td) / Tn
```

where `Tn` is native polynomial-evaluation time, and `Tp`, `Tv`, and `Td` are the client-side ProbGen, Verify, and Decode/Reconstruct times. KeyGen is excluded because it can be amortized. Figures 2--4 plot `log2(Rt*)`.

| Paper result | Schemes and parameter grid | Data file |
|---|---|---|
| Figure 2 (`fig:Output1`) | Gamma1--3; `d=2`, `t=1`, `m=200,400,...,2000` | `results/reported/figure2_low_degree.csv` |
| Figure 3 (`fig:Rt`) | Gamma1--2; `d=2`, `t=2,3`, `m=200,400,...,2000` | `results/reported/figure3_thresholds.csv` |
| Figure 4 (`fig:bigd`) | Gamma1--3; `d=4,8,10`, `t=1`, `m=1,...,15` | `results/reported/figure4_high_degree.csv` |
| Table `tab:micro` | 100 trials per primitive | `results/reported/table_microbenchmark.csv` |
| Figure 5 (`fig:tc_gamma1_pi3`) | Gamma1/Pi3; `d=2,3`, `t=1,2`, `m=200,400,...,2000` | `results/reported/figure5_gamma1_vs_pi3.csv` |

The files under `results/reported/` provide the numerical coordinates used by the plotting script.

## Dense polynomial representation

`lib/common/mpoly.c` defines one coefficient order for a dense polynomial in `m` variables with total degree at most `d`:

```text
degree 0: 1
degree 1: x0, x1, ..., x(m-1)
degree 2: x0^2, x0*x1, ..., x(m-1)^2
degree 3: x0^3, x0^2*x1, ..., x(m-1)^3
...
```

Within each degree, a monomial is represented by a nondecreasing tuple of variable indices. The number of coefficients is `C(m + d, d)`.

The same graded order is used by the native evaluator, Gamma1 extension-field evaluator, Gamma2 through the shared evaluator, and the Pi3 baseline. The regression test covers a complete two-variable degree-3 polynomial in all three evaluator variants.

Run the tests directly with:

```bash
make test
```

## Build targets

| Command | Purpose |
|---|---|
| `make all` | Build core objects and the two exploratory PIR programs |
| `make core` | Build only the construction and common-library objects |
| `make test` | Build and run polynomial evaluator regression tests |
| `make validate` | Validate all Figure 2--5 grids and source hygiene |
| `make check` | Run validation, tests, and the complete build |
| `make var_benchmark` | Build the exploratory VAR estimator |
| `make reported-plots` | Replot the paper's reported coordinates |
| `make clean` | Remove generated native binaries and objects |

The Makefile detects common Homebrew FLINT locations on macOS. The artifact is also checked on Ubuntu through `.github/workflows/ci.yml`.

## Reported experimental platform

- Client VM: Ubuntu Desktop 20.04.2 LTS, 4 GB RAM, one Intel i7-6700 core limited to 800 MHz.
- Server VM: Ubuntu Desktop 20.04.2 LTS, 32 GB RAM, one 2.30 GHz Intel Xeon Gold 5218 core.
- The VMs used bridged networking on the same physical workstation and communicated over TCP.
- Default security parameter: 128 bits.
- Reported software versions: FLINT 2.8.0 and libsodium 1.0.18.
- Reported group: ristretto255 through libsodium.

Current systems may use newer compatible library versions. Record the actual versions and machine configuration with any newly generated timing results.

For the shortest path from protocol names to source files and Figure 2--5 commands, see `REVIEWER_GUIDE.md`.
