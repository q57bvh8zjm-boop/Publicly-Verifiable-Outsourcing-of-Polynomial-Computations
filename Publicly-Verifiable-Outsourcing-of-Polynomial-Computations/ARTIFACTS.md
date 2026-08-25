# Artifacts for paper evaluation

## Core scheme implementations

- **$\Gamma_1$ -- field-extension public verification:** `protocols/scheme1.c`, `protocols/scheme1.h`
- **$\Gamma_2$ -- Shamir/orthogonal-vector public verification:** `protocols/scheme2.c`, `protocols/scheme2.h`
- **$\Gamma_3$ -- two-server HSS construction:** `protocols/scheme3.c`, `protocols/scheme3.h`
- **$\Pi_3$ private-verification baseline:** `protocols/protocol2.c`, `protocols/protocol2.h`

## Shared libraries

- Multivariate-polynomial utilities: `lib/common/mpoly.c`, `lib/common/mpoly.h`
- Explicit quadratic extension-field arithmetic: `lib/common/extended_field.c`, `lib/common/extended_field.h`

## Experimental results

The Figure 2--5 and table datasets are stored in `results/reported/`. Use the plotting script to generate the corresponding images:

```bash
make reported-plots
```

See `SOURCE_CODE_MAP.md` for the source-file, parameter, result-file, and paper-figure correspondence.

## Application programs

- `var_benchmark.c`: auditable-aggregation/VAR operation estimator
- `protocols/scheme3_pir_optimize.c`: HSS/PIR parameter-cost estimator
- `protocols/protocol2_pir_optimize.c`: field-extension PIR baseline estimator

## Build and validation

```bash
make check
```
