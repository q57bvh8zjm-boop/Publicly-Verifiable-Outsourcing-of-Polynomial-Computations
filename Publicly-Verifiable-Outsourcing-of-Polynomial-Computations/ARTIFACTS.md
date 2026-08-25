# Artifacts for paper evaluation

## Core scheme implementations

- **$\Gamma_1$ — field-extension public verification**: `protocols/scheme1.c`, `protocols/scheme1.h`
- **$\Gamma_2$ — Shamir/orthogonal-vector public verification**: `protocols/scheme2.c`, `protocols/scheme2.h`
- **$\Gamma_3$ — two-server HSS construction**: `protocols/scheme3.c`, `protocols/scheme3.h`
- **$\Pi_3$ private-verification baseline** used in the comparison with $\Gamma_1$: `protocols/protocol2.c`, `protocols/protocol2.h`

The paper proposes three new schemes. A previous version of this artifact file incorrectly listed nonexistent `scheme4.c/h`; that stale entry has been removed.

## Shared libraries

- Multivariate-polynomial utilities: `lib/common/mpoly.c`, `lib/common/mpoly.h`
- Explicit quadratic extension-field arithmetic: `lib/common/extended_field.c`, `lib/common/extended_field.h`

A previous version listed `lib/common/he.c/h`; those files are not present. HSS/Paillier-related code is inside `protocols/scheme3.c`.

## Experimental-result mapping

See `SOURCE_CODE_MAP.md`. The paper's printed numerical coordinates are under `results/reported/` and can be replotted with `python3 scripts/plot_reported_results.py`.

These files reproduce the **reported figures** and are deliberately kept separate from raw benchmark data. The original end-to-end timing drivers for Figures 2--5 are not present; see `REPRODUCIBILITY_STATUS.md`.

## Exploratory application programs

- `var_benchmark.c`: auditable-aggregation/VAR atomic-operation estimator; not an end-to-end execution of all three schemes.
- `protocols/scheme3_pir_optimize.c`: HSS/PIR parameter-cost estimator.
- `protocols/protocol2_pir_optimize.c`: field-extension baseline PIR parameter-cost estimator.

## Build

```bash
make all
```
