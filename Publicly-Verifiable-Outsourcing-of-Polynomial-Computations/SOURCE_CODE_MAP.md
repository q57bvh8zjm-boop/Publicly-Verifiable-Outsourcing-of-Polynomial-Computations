# Source-code to paper-experiment map

This file answers the practical artifact-evaluation question: **which source file is related to which construction, figure, or table in the paper?**

## Core constructions

| Repository file | Paper object | What it implements | Status |
|---|---|---|---|
| `protocols/scheme1.c`, `scheme1.h` | $\Gamma_1$ | Publicly verifiable MSVC using field-extension verification over $\mathbb F_{q^2}$ | Core implementation present |
| `protocols/scheme2.c`, `scheme2.h` | $\Gamma_2$ | Publicly verifiable MSVC using $d+1$ Shamir-shared curves and orthogonal-vector verification | Core implementation present |
| `protocols/scheme3.c`, `scheme3.h` | $\Gamma_3$ | Two-server HSS-based publicly verifiable MSVC | Core implementation present; HSS port contains assumptions documented in source comments |
| `protocols/protocol2.c`, `protocol2.h` | $\Pi_3$ of Zhang et al. | Private-verification field-extension baseline used for the $\Gamma_1$ comparison | Baseline evaluator supports arbitrary total degree, including Figure 5's $d=3$ case |
| `lib/common/mpoly.c`, `mpoly.h` | Shared polynomial utilities | Base-field evaluation and interpolation | Dense graded evaluation supports arbitrary total degree and exposes the coefficient count |
| `lib/common/extended_field.c`, `.h` | Shared $\mathbb F_{q^2}$ utilities | Extension-field arithmetic, matrices, polynomial helpers and operation counters | Present |

## Main experimental figures

The current C tree does **not** contain the original end-to-end benchmark drivers that produced Figures 2--5. The numerical coordinates printed in the paper are therefore included under `results/reported/` and can be replotted with `scripts/plot_reported_results.py`.

| Paper result | Parameters | Reported-data file | Raw benchmark driver in current repository? |
|---|---|---|---|
| Figure 2 (`fig:Output1`) | $d=2,t=1,m=200,400,\ldots,2000$; $\Gamma_1,\Gamma_2,\Gamma_3$ | `results/reported/figure2_low_degree.csv` | **No** |
| Figure 3 (`fig:Rt`) | $d=2,t\in\{2,3\},m=200,\ldots,2000$; $\Gamma_1,\Gamma_2$ | `results/reported/figure3_thresholds.csv` | **No** |
| Figure 4 (`fig:bigd`) | $d\in\{4,8,10\},t=1,m=1,\ldots,15$; all three schemes | `results/reported/figure4_high_degree.csv` | **No; Gamma3's HSS evaluator remains linear-only and the original driver is absent** |
| Table `tab:micro` | 100-trial client-side microbenchmarks | `results/reported/table_microbenchmark.csv` | **No exact driver** |
| Figure 5 (`fig:tc_gamma1_pi3`) | $d\in\{2,3\},t\in\{1,2\},m=200,\ldots,2000$ | `results/reported/figure5_gamma1_vs_pi3.csv` | **No original timing driver; both evaluators now support $d=3$** |

For Figures 2--4 the plotted quantity is

$$\operatorname{Rt}^*_{\mathcal{VC}}(d,t,m)=\frac{T_p+T_v+T_d}{T_n},$$

and the paper plots $\lg(\operatorname{Rt}^*)$.

## Application-oriented / exploratory programs

| File | Intended role | Relation to paper | Important caveat |
|---|---|---|---|
| `var_benchmark.c` | Atomic-operation estimator for the auditable-aggregation/VAR application | Related to Table `VAR` | **Not** an end-to-end execution of $\Gamma_1$--$\Gamma_3$; the $\Gamma_3$ issuance number is extrapolated |
| `protocols/scheme3_pir_optimize.c` | Cost estimator for HSS/PIR parameter choices | Related to the PIR application discussion | Uses database size $2^{20}$ and analytical scaling; it does not directly reproduce Table `PIR_i` (paper: $n=2^{26}$) |
| `protocols/protocol2_pir_optimize.c` | Cost estimator for the field-extension baseline | Baseline PIR exploration | Not a direct reproduction of the $\Pi_4/\Pi_5$ rows in Table `PIR_i` |

