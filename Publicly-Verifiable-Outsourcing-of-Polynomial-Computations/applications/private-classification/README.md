# Verifiable private polynomial classification

This application evaluates fixed public quadratic and cubic SVM spam classifiers with the existing Gamma_FE_Sh and Gamma_OV_Sh protocol APIs. One client holds the private query and applies the decision threshold after score verification and signed decoding. Servers evaluate their own shares independently. At threshold one, the quadratic and cubic workloads use three and four servers.

## Reproduce

From the artifact root, use Python 3.9–3.12 with the pinned packages and the FLINT, GMP, and libsodium C dependencies:

```sh
python3 -m pip install -r applications/private-classification/requirements.txt
make classification-check
```

The default output directory is `results/measured/local-classification/`. Set `CLASSIFICATION_OUTPUT=results/measured/local-classification-2` for another run. Each run uses a new output directory. The runner checks all 48 timed queries against an exact integer oracle, tests 96 single-server response perturbations, and checks all 1,381 native test predictions for each model against the floating-point labels. Results and their metadata are written to the new directory. The original paper measurements remain in `reported/`; the raw CSVs retain all measured fields, including native arithmetic references.

To rebuild the model fixtures from the included public dataset, use Python 3.9 with the pinned requirements:

```sh
python3 applications/private-classification/prepare.py --output results/measured/classification-models
python3 applications/private-classification/run.py --fixtures results/measured/classification-models/data --models results/measured/classification-models/results/models.json --output results/measured/classification-rebuilt
```

The included fixtures fix the exact reported model. Rebuilding may change the least significant coefficient bits through floating-point training; the runner verifies the rebuilt integer oracle and prediction labels.

The fixed split uses seed 20261028: 3,220 training and 1,381 test emails, each with 57 features. Preprocessing is log1p, training-fitted standardization, and clipping to [-4,4]. SVM parameters are C=1, gamma=1/57, coef0=1, and degree 2 or 3. Linear and RBF models provide task-quality references using the same split and preprocessing.

## Encoding and measured results: Appendix G and Table 7

The quadratic score is compiled into 1,711 polynomial coefficients. Inputs use scale 4096; degree-r coefficients use scale 2^24 * 4096^(2-r). The cubic score retains its 600 support vectors, with input/support-vector scale 1024, weight scale 2^20, and score scale 2^20 * (57*1024^2)^3. Both run exactly over the Ristretto scalar field, with no intermediate rounding. The recorded integer score bounds are less than q/2. Test accuracy is 94.50% and 94.93%, with no label changes after quantization.

| Degree | Scheme | Servers | Client (ms) | Verify (ms) | Request bytes | Reply bytes | Verification-key bytes |
|---:|---|---:|---:|---:|---:|---:|---:|
| 2 | Gamma_FE_Sh | 3 | 0.558 | 0.184 | 11136 | 192 | 192 |
| 2 | Gamma_OV_Sh | 3 | 0.634 | 0.122 | 16416 | 288 | 96 |
| 3 | Gamma_FE_Sh | 4 | 0.663 | 0.263 | 14848 | 256 | 256 |
| 3 | Gamma_OV_Sh | 4 | 0.924 | 0.169 | 29184 | 512 | 128 |

Times are medians over 12 queries after one warmup on a macOS 15.6.1 ARM64 host. Client time includes field conversion, ProbGen, request serialization, Verify, Decode, and client-object cleanup. Verify is part of client time. It starts with already quantized features; training, model loading, feature extraction, network transport, and response parsing are outside the interval. Parties execute sequentially. Request/reply counts use 32-byte field scalars, exclude headers, and report verification-key traffic separately.

The adapter links the main protocol source files unchanged and implements application-specific evaluation on shares. Model parameters are public. Input privacy is server-local; the complete response transcript stays with the client and its authorized auditor. The existing trusted-query-generation and protocol security assumptions apply. The scope is exact score computation under the chosen fixed model.

## Data and provenance

Spambase: Mark Hopkins, Erik Reeber, George Forman, and Jaap Suermondt (1999), UCI Machine Learning Repository, [DOI 10.24432/C53G6X](https://doi.org/10.24432/C53G6X), licensed CC BY 4.0. The archive is included unchanged; fixtures and model files are derived from it through the documented preprocessing and training. `reported/models.json` records the split and integer oracle; `reported/run-metadata.json` records the original executable and core source hashes. Table 6 compares protocol properties and server counts across EPIC, vPIN, Conan, Pi3/Pi4/Pi5, and the Shamir constructions. The measured application costs above correspond to Appendix G/Table 7.
