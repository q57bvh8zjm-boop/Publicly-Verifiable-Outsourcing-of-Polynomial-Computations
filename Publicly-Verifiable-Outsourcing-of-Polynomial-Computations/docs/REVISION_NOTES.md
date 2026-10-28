# Implementation and validation

| Component | Implementation | Validation |
|---|---|---|
| Polynomial evaluation | Shared graded monomial order, including constant and mixed high-degree terms. | `make check` exercises degrees 1–10 and thresholds 1–3. |
| HSS evaluation | Digit-decomposed Paillier backend with integer lifts and complete monomial traversal. | `make test-hss` covers quadratic and mixed cubic terms, full-field inputs, constants, and malformed ciphertexts; `make test-hss-high` covers powers 4–10 at `q-1`. |
| Randomness | OS-backed rejection sampling across the full field; distinct nonzero scaling points. | Protocol tests check fresh query shares; sampler construction supplies uniformity. |
| Public verification | FE uses an encoded public key; OV verification consumes the public OV key and field responses. | End-to-end tests check decoded equality, altered responses, and malformed points. |
| Extension arithmetic | Multiplication computes both coordinates before writing the output. | FE tests exercise aliased operations through complete evaluations. |
| Private classification | Fixed quadratic and cubic SVM scores over quantized Spambase features. | `make classification-check` verifies 48 protocol queries, 96 rejected response perturbations, and all test labels. |
| Measurements | Executable benchmarks record phase times, raw data, source hashes, and environment metadata. | Timing rows are emitted after verification and native equality checks. |

The private classification implementation corresponds to Appendix G/Table 7. Table 6 presents the protocol-feature comparison, including EPIC, vPIN, Conan, and the polynomial-outsourcing schemes. The paper's other applications describe polynomial mappings; the executable application in this release is private classification.

The HSS backend's assumptions, parameter bounds, and post-evaluation-key adaptive-correctness requirement are documented in [HSS_INSTANTIATION.md](HSS_INSTANTIATION.md). Historical plot transcriptions and their provenance are documented in [results/reported/README.md](../results/reported/README.md). Current measurement boundaries are in [MEASUREMENTS.md](MEASUREMENTS.md).
