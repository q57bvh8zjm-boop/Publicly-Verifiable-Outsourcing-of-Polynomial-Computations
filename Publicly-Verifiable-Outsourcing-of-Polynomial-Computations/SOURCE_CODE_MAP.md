# Source-code to paper-experiment map

This document maps each implementation file and result dataset to the corresponding construction, figure, or table in the paper.

## Core constructions

| Repository file | Paper object | Purpose |
|---|---|---|
| `protocols/scheme1.c`, `scheme1.h` | $\Gamma_1$ | Publicly verifiable MSVC using field-extension verification over $\mathbb F_{q^2}$ |
| `protocols/scheme2.c`, `scheme2.h` | $\Gamma_2$ | Publicly verifiable MSVC using $d+1$ Shamir-shared curves and orthogonal-vector verification |
| `protocols/scheme3.c`, `scheme3.h` | $\Gamma_3$ | Two-server HSS-based publicly verifiable MSVC |
| `protocols/protocol2.c`, `protocol2.h` | $\Pi_3$ | Private-verification field-extension baseline used in the $\Gamma_1$ comparison |
| `lib/common/mpoly.c`, `mpoly.h` | Polynomial utilities | Dense graded polynomial evaluation, coefficient counting, and interpolation |
| `lib/common/extended_field.c`, `.h` | Extension-field utilities | $\mathbb F_{q^2}$ arithmetic, matrices, polynomial helpers, and operation counters |

## Main experimental figures

| Paper result | Parameters | Result data |
|---|---|---|
| Figure 2 (`fig:Output1`) | $d=2,t=1,m=200,400,\ldots,2000$; $\Gamma_1,\Gamma_2,\Gamma_3$ | `results/reported/figure2_low_degree.csv` |
| Figure 3 (`fig:Rt`) | $d=2,t\in\{2,3\},m=200,\ldots,2000$; $\Gamma_1,\Gamma_2$ | `results/reported/figure3_thresholds.csv` |
| Figure 4 (`fig:bigd`) | $d\in\{4,8,10\},t=1,m=1,\ldots,15$; $\Gamma_1,\Gamma_2,\Gamma_3$ | `results/reported/figure4_high_degree.csv` |
| Table `tab:micro` | 100 trials per primitive | `results/reported/table_microbenchmark.csv` |
| Figure 5 (`fig:tc_gamma1_pi3`) | $d\in\{2,3\},t\in\{1,2\},m=200,\ldots,2000$ | `results/reported/figure5_gamma1_vs_pi3.csv` |

For Figures 2--4, the plotted quantity is

$$\operatorname{Rt}^*_{\mathcal{VC}}(d,t,m)=\frac{T_p+T_v+T_d}{T_n},$$

and the vertical axis reports $\lg(\operatorname{Rt}^*)$.

## Application programs

| File | Purpose | Paper relation |
|---|---|---|
| `var_benchmark.c` | Auditable-aggregation/VAR operation estimator | Table `VAR` |
| `protocols/scheme3_pir_optimize.c` | HSS/PIR parameter-cost estimator | PIR application evaluation |
| `protocols/protocol2_pir_optimize.c` | Field-extension PIR baseline estimator | PIR baseline evaluation |

## Plot generation

```bash
make reported-plots
```

Generated images are written to `results/plots/`.
