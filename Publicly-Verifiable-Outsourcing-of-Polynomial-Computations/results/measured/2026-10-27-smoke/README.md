# Measured smoke run — 2026-10-27

This run contains 36 successful trials from 12 configurations on macOS arm64. Each configuration has three trials. `metadata.json` records the precise platform, compiler, commands, library output files, and measured source hashes.

The dense and VAR kernels at `m=2,d=2,t=1` cover Gamma_FE_Sh, Gamma_OV_Sh, and Gamma_OV_HSS. Gamma_FE_Sh and Gamma_OV_Sh also cover `m=2,t=1,d=4,8,10`. All rows passed verification and native equality. The VAR rows additionally matched a direct matrix-vector reference.

The sample is too small and the workload too small to establish paper-scale performance or crossover points. Gamma_OV_HSS has substantial query-generation and server cost in this backend. Setup time varies because safe primes are generated freshly; it is recorded separately. HSS regression tests and sanitizer runs finished before these measurements began.

`raw.csv` combines the per-command CSVs without aggregation. Each `.stderr` contains library versions. `validation.txt` records the local test outcomes. See `../../../docs/MEASUREMENTS.md` for phase and payload definitions.

Original data and run metadata retain the labels `Gamma1`, `Gamma2`, and `Gamma3`, corresponding to `Gamma_FE_Sh`, `Gamma_OV_Sh`, and `Gamma_OV_HSS` in the current manuscript. These records predate the naming update; their recorded commands remain usable through numeric CLI aliases.
