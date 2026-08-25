# Reviewer Quick Start

This is the only file needed for a quick artifact review. It explains which source files implement each protocol and gives one-command instructions for generating the paper figures.

## 1. Install dependencies

Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install build-essential libflint-dev libgmp-dev libsodium-dev python3 python3-pip
python3 -m pip install matplotlib
```

## 2. Build and verify the artifact

From the repository root, run:

```bash
make check
```

This command validates the result datasets, builds all protocol source files, and runs the polynomial-evaluation tests.

## 3. Protocol-to-file mapping

| Protocol in the paper | Main source files | Role |
|---|---|---|
| $\Gamma_1$ | `protocols/scheme1.c`, `protocols/scheme1.h` | Field-extension publicly verifiable construction |
| $\Gamma_2$ | `protocols/scheme2.c`, `protocols/scheme2.h` | Shamir/orthogonal-vector publicly verifiable construction |
| $\Gamma_3$ | `protocols/scheme3.c`, `protocols/scheme3.h` | Two-server Paillier/HSS construction |
| $\Pi_3$ | `protocols/protocol2.c`, `protocols/protocol2.h` | Private-verification baseline used in the $\Gamma_1/\Pi_3$ comparison |
| $\Pi_4$ PIR experiment | `protocols/protocol2_pir_optimize.c` | Field-extension PIR parameter and cost experiment |
| $\Pi_5$ PIR experiment | `protocols/scheme3_pir_optimize.c` | HSS/PIR parameter and cost experiment |
| Shared polynomial code | `lib/common/mpoly.c`, `lib/common/mpoly.h` | Dense multivariate evaluation and interpolation |
| Shared extension field | `lib/common/extended_field.c`, `lib/common/extended_field.h` | Arithmetic over the explicit quadratic extension field |

The filenames retain the implementation numbering used by the source code. The table above gives their paper notation directly.

## 4. Generate the paper figures

Run any of the following commands from the repository root. Each command creates the corresponding PNG file under `results/plots/`.

| Paper result | Protocols and parameters | Command | Generated file |
|---|---|---|---|
| Figure 2 (`fig:Output1`) | $\Gamma_1,\Gamma_2,\Gamma_3$; $d=2,t=1,m=200,400,\ldots,2000$ | `make figure2` | `results/plots/figure2_low_degree.png` |
| Figure 3 (`fig:Rt`) | $\Gamma_1,\Gamma_2$; $d=2,t=2,3,m=200,400,\ldots,2000$ | `make figure3` | `results/plots/figure3_thresholds.png` |
| Figure 4 (`fig:bigd`) | $\Gamma_1,\Gamma_2,\Gamma_3$; $d=4,8,10,t=1,m=1,\ldots,15$ | `make figure4` | `results/plots/figure4_d4.png`, `figure4_d8.png`, `figure4_d10.png` |
| Figure 5 (`fig:tc_gamma1_pi3`) | $\Gamma_1/\Pi_3$; $d=2,3,t=1,2,m=200,400,\ldots,2000$ | `make figure5` | `results/plots/figure5_gamma1_vs_pi3.png` |

Generate all figures at once:

```bash
make figures
```

The figure commands use the corresponding CSV datasets in `results/reported/`:

```text
Figure 2 -> results/reported/figure2_low_degree.csv
Figure 3 -> results/reported/figure3_thresholds.csv
Figure 4 -> results/reported/figure4_high_degree.csv
Figure 5 -> results/reported/figure5_gamma1_vs_pi3.csv
```

## 5. Run the PIR experiments

Build and run the two PIR experiment programs:

```bash
make all
./build/protocol2_pir_optimize
./build/scheme3_pir_optimize
```

The first command corresponds to the $\Pi_4$ field-extension PIR experiment. The second corresponds to the $\Pi_5$ HSS/PIR experiment.

## 6. Other useful commands

```bash
make test             # Run polynomial evaluator tests
make validate         # Validate all Figure 2--5 datasets
make var_benchmark    # Build the VAR application benchmark
./build/var_benchmark # Run the VAR application benchmark
make clean            # Remove generated binaries and object files
```
