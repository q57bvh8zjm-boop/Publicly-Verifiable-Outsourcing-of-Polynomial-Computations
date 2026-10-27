# Revision responding to artifact review

| Finding | Change | Evidence / remaining work |
|---|---|---|
| HSS ignored quadratic and higher terms | Replace the linear-only implementation with a digit-decomposed Paillier backend and complete graded monomial traversal. | `test-hss` includes `x^2` at 2, mixed cubic terms, full-field inputs/coefficients, zero and constant polynomials. `test-hss-high` checks degrees 4–10 at `q-1`. |
| Shamir masks came from 1–500; scaling points came from a small range | Sample coefficients uniformly by OS-backed rejection sampling over the full field; sample distinct nonzero scaling points over the full field. | No small-range sampler or unseeded FLINT RNG remains in Gamma_FE_Sh, Gamma_OV_Sh, and Gamma_OV_HSS. Tests exercise fresh query shares. Distribution follows from the sampler, not a frequency-test claim. |
| Verification depended on client secrets | Give Gamma_FE_Sh a per-query public key containing encoded powers; give Gamma_OV_Sh and Gamma_OV_HSS verification functions only the public OV key and returned field values. | End-to-end tests use those interfaces, including rejected changes and malformed points. |
| Extension multiplication failed for aliased operands | Compute both coordinates before writing either output. | End-to-end FE tests repeatedly use in-place multiplication at degrees 1–10 and thresholds 1–3. |
| Old figure CSVs were presented as reproducibility evidence | Preserve them as explicitly archived transcriptions; add executable benchmarks and raw-log collection. | Measured data lives in `results/measured/`. The full paper grid and original VM platform have not been rerun by the smoke experiment. |
| VAR source was absent; authentication/proof costs were unsupported | Replace the broken target with a measured VAR polynomial kernel and an independent direct reference. Remove obsolete PIR cost estimators. | No complete receipt/signature/binding-proof implementation is claimed. Those costs require implementations and fresh application benchmarks. |
| HSS assumptions and range were not tied to code | Document the exact backend, integer lifts, setup trust, component budget, modulus sizing, and KDM/PRF assumptions. | Post-evaluation-key adaptive correctness remains a separate proof obligation. |

The historical Figure 4 data has Gamma_OV_HSS's `d=4,m=14` point below zero, while the supplied manuscript says Gamma_OV_HSS remains slower through `m=15`. The archived coordinates are preserved verbatim so the disagreement remains auditable. They must be replaced with measured results before using the figure to support the revised implementation.

## Manuscript changes still required

The repository contains no TeX source for the supplied manuscript. Update its implementation section to identify the new HSS backend and its KDM assumption. Use the range calculation from the instantiation note with its stated model. Keep the stronger adaptive condition explicit until it has a proof. Replace performance claims based on historical data with new runs on a documented platform, and describe the current VAR result as a polynomial-kernel experiment. Full application overhead requires a concrete authentication and binding-proof system.
