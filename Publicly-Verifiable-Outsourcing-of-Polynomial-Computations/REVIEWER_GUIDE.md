# Reviewer quick start

Run these commands from `Publicly-Verifiable-Outsourcing-of-Polynomial-Computations/` with FLINT, GMP, libsodium, a C compiler, and Python 3 installed.

## Core protocols

```sh
make check
make test-hss
make test-hss-high
```

The fast suite checks Gamma_FE_Sh and Gamma_OV_Sh at degrees 1–10 and thresholds 1–3: native equality, fresh shares, altered outputs, malformed verification points, zero polynomials, and constants. HSS tests add deterministic repeated evaluation, malformed ciphertext rejection, and high-degree inputs at `q-1`. HSS key generation and evaluation can take several minutes.

## Private classification: Appendix G and Table 7

With Python 3.9–3.12:

```sh
python3 -m pip install -r applications/private-classification/requirements.txt
make classification-check
```

Expected checks are 48 correct protocol queries, 96 rejected response perturbations, and zero changed labels after quantization on all 1,381 test samples for each model. Inspect `summary.json`, the four CSVs, prediction TSVs, and `run-metadata.json` in `results/measured/local-classification/`. To rerun, choose a fresh `CLASSIFICATION_OUTPUT` directory.

The [application README](applications/private-classification/README.md) gives the fixed split, original measurements, model regeneration command, and timing boundaries. Table 6 compares protocol properties and server counts; the application measurements are reported in Appendix G/Table 7.

## Dense polynomial measurements

```sh
python3 scripts/run_experiments.py --profile smoke --trials 3 --output results/measured/local-review
```

Inspect `raw.csv`, per-command files, and `metadata.json`. The benchmark records setup, query generation, server work, verification, reconstruction, and native evaluation, and stops on a failed correctness check.

[MEASUREMENTS.md](docs/MEASUREMENTS.md) defines the measured phases and full parameter grid. [HSS_INSTANTIATION.md](docs/HSS_INSTANTIATION.md) states the backend's integer-lifting model, parameter bound, setup trust, and security and correctness conditions. [REVISION_NOTES.md](docs/REVISION_NOTES.md) maps each implemented capability to its validation.

`make reported-plots` redraws the archived paper coordinates described in [results/reported/README.md](results/reported/README.md). Raw execution measurements are stored separately from these transcriptions.
