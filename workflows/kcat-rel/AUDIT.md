# Frozen evidence audit

Audit date: 2026-09-16

| Gate | Result | Evidence |
| --- | --- | --- |
| Frozen primary cohort | PASS | Exactly 17 unique canonical variants; D382Y and P710R absent |
| One label per variant | PASS | 17 canonical rows from 22 linked measurement rows |
| Eligible biological assay | PASS | Pure-actin steady-state actin-activated `kcat` only |
| Accepted construct rule | PASS | 10 sS1/S1 labels; 7 fallback 2-hep labels |
| Matched WT context | PASS | Every linked row has identical WT/mutant context IDs |
| Provenance | PASS | Source study, DOI, source location, construct, chemistry, and evidence epoch retained |
| Exact label math | PASS | `kcat_rel` and `ln(kcat_rel)` recomputed at absolute tolerance `1e-12` |
| Uncertainty | PASS with stated limitations | Reported error types retained; propagated columns audited |
| Pseudoreplication | PASS | Morck replicates remain separate and collapse to five variant-level labels |
| Study grouping | PASS | Nine source-study IDs retained for later grouped validation |

The machine-readable gate is `src/kcat_rel/evidence.py`. The retained
limitations are listed in [`README.md`](README.md); this audit does not approve
features, fit a diagnostic model, or authorize forward prediction.

## Frozen fallback-feature diagnostic (2026-09-16)

| Gate | Result | Evidence |
| --- | --- | --- |
| P12883 input pin | PASS | Raw FASTA SHA-256 `45c96586dd50e0ede770a3bb4fb0aab125d103d53b2bc0b0fed89619e07a16ad`; normalized sequence SHA-256 `e98e5d01a359820bcd1a00c3422025d9fedf1fdbc46ab586e3f5b731f4585be6` |
| 8ACT input pin | PASS | Compressed mmCIF SHA-256 `69bc81d3b6ee01a09e300a0e9b65f5a02837bbda9394e5566b248bf612a1b5e0` |
| Label-blind fallback features | PASS | 17 ordered rows and three finite features; audited before labels were opened |
| Ridge primary gate | FAIL | Study macro MAE +12.52% worse, variant MAE +11.09% worse; both RMSE limits and deletion influence also fail |
| GPR primary gate | FAIL | Study macro MAE +7.85% worse, variant MAE +9.52% worse; both RMSE limits and deletion influence also fail |
| Forward prediction | NOT AUTHORIZED | Neither candidate passed every predeclared gate; no winner was selected |

The complete deterministic artifact set, including fold-local tuning, OOF
predictions, a per-split inner-fold audit, secondary uncertainty sensitivities, input hashes, package
versions, command, fallback reason, and output hashes is in
[`results/diagnostic-v0`](results/diagnostic-v0). `SUCCESS` was written last to
mark completed reporting, not a positive model decision.

The frozen run's `inner-fold-audit.json` is an artifact-only
`post_execution_reconstruction` anchored to execution revision
`91affadb0f5bdb3fbf8f09b22209836f38ddf9c5`. It reconstructs every nested
source holdout and its scaler from inner-training rows only; the diagnostic was
not rerun and the six frozen scientific result files were not changed.
