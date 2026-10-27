# Archived reported coordinates

These CSVs were inherited from the earlier repository version and transcribed from paper plots/tables. They are preserved for historical comparison. They are not raw timing logs and do not measure the repaired implementation.

Files: `figure2_low_degree.csv`, `figure3_thresholds.csv`, `figure4_high_degree.csv`, `figure5_gamma1_vs_pi3.csv`, `table_microbenchmark.csv`, and `table_var.csv`.

Known mismatch: the archived Gamma_OV_HSS `d=4,m=14` point in Figure 4 has a negative `log2_Rt_star`, whereas the supplied manuscript says it remains slower than native evaluation for `m<=15`. The original hardware runs and provenance needed to resolve that mismatch are not present. The VAR table also includes application components that this repository does not implement.

`make reported-plots` redraws these coordinates with an archival label. New measurements are stored separately under `results/measured/`; see `docs/MEASUREMENTS.md`.

Original data and run metadata retain the labels `Gamma1`, `Gamma2`, and `Gamma3`, corresponding to `Gamma_FE_Sh`, `Gamma_OV_Sh`, and `Gamma_OV_HSS` in the current manuscript. These records predate the naming update; their recorded commands remain usable through numeric CLI aliases.
