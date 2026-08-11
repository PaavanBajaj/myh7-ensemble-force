# Continuous LSAR (`lsar-na-rel`) — method

First vertical slice. Modeled biochemical quantity: continuous LSAR
\(\mathrm{LSAR}=k_{\mathrm{cat,long-tail}}/k_{\mathrm{cat,short-tail}}\) on paired
two-headed human β-cardiac HMM constructs.

Canonical force-facing relative-head label:

\[
N_{a,\mathrm{rel}} = \frac{\mathrm{LSAR}_{mut}}{\mathrm{LSAR}_{WT,\mathrm{study}}},
\qquad
y = \ln(N_{a,\mathrm{rel}}).
\]

\(y\) is the intended GPR target. Global \(\mathrm{LSAR}/0.57\) and the release
transform \(a=(\mathrm{LSAR}_{mut}-\mathrm{LSAR}_{WT})/(1-\mathrm{LSAR}_{WT})\)
are retained only as separate derived/audit fields.

## Status (S1)

Provenance-ready public derivatives are published under
`data/public/lsar-na-rel/` with source identifiers in
`data/public/manifests/lsar-na-rel-sources.json`. The primary pipeline table is
`na-rel-primary.csv`. No fitted model and no pinned regeneration pipeline yet
(S2+).

## Inclusion criteria

Exact continuous-LSAR candidates require:

- human β-cardiac myosin / MYH7
- the evaluated amino-acid variant
- paired short-tail and long-tail two-headed HMM from an accepted family
  (25-hep/2-hep or 15-hep/8-hep)
- actin-activated steady-state ATPase fitted \(k_{cat}\)
- exact numeric primary values (text, table, SI, or first-party numeric source data)

## Input LSAR policy

1. Prefer the paper’s reported LSAR (and printed uncertainty) when available.
2. Otherwise use the exact quotient of printed mean \(k_{cat}\) values
   (long/short), with uncertainty propagated from printed \(k_{cat}\) errors.
3. Study WT LSAR: prefer the paper’s reported WT LSAR ± error; otherwise the
   exact WT long/short quotient from printed means ± propagated error. Morck WT
   uses the paper’s stated WT 2-/25-hep ratio \(0.61\pm0.10\).
4. Morck mutants use the paper’s reported average LSAR across two biological
   replicates.

## Derived quantities (kept distinct)

1. Continuous LSAR (biochemical input).
2. Within-study \(N_{a,\mathrm{rel}}=\mathrm{LSAR}_{mut}/\mathrm{LSAR}_{WT}\)
   and \(y=\ln(N_{a,\mathrm{rel}})\).
3. Global Spudich-baseline ratio \(\mathrm{LSAR}/0.57\)
   (`na_rel_global_spudich057`).
4. Spudich additional-\(N_a\) percent: \(100(\mathrm{LSAR}-0.57)/0.43\).
5. Release (audit-only):
   \(a=(\mathrm{LSAR}_{mut}-\mathrm{LSAR}_{WT})/(1-\mathrm{LSAR}_{WT})\).
6. Table-implied LSAR from rounded Table 1 percent:
   \(0.57+0.43(p/100)\) — audit-only; never a training label.

Raw LSAR and within-study transforms are not capped. Apparent Table 1 display
caps (and near-WT coding of negatives as 0) are recorded only as display fields.

## Uncertainty propagation

For a ratio \(R=X/Y\) with independent relative errors:

\[
\frac{\sigma_R}{R}=\sqrt{\left(\frac{\sigma_X}{X}\right)^2+\left(\frac{\sigma_Y}{Y}\right)^2},
\qquad
\sigma_{\ln R}=\frac{\sigma_R}{R}.
\]

Error types remain as reported (SEM, SE of fit, or unspecified) and are not
re-labeled as a common SEM.

## Evidence handling

- Preserve reported LSAR separately from recomputed long/short.
- One canonical label row per variant; biological replicates live only in
  `measurements.csv`.
- Primary pipeline table (`na-rel-primary.csv`): 13 Primary + 3 Risky allowed.
- G256E: approximate relative/figure evidence; excluded from primary Na_rel.
- Q222K, E536D, V606M: unresolved unpublished; `exclude_unresolved`.
- D239N, R453C, G741R: assay absent, not zero; `exclude_blank`.
- G768R: later Pathak published epoch is canonical; Spudich 2024 unpublished
  epoch remains unresolved and separate.
- R403Q / R663H: exact pairs with low biological \(n\) (`review_low_n`; Risky
  allowed in the primary Na_rel table).
- Preprint and final versions of the same summary are one evidence identity.

## Construct / assay covariates (pooling review before S2+)

Do not silently pool across:

- 2/25-hep versus 8/15-hep constructs
- colorimetric versus NADH-coupled assays
- temperature/buffer differences
- study-level WT baselines
- preparation-level averaging versus ratio-of-means
- low biological replication

## Public pack

See `data/public/lsar-na-rel/README.md` and `schema.json` for file roles,
eligibility vocabulary, checksums, and the licensing limitation.
