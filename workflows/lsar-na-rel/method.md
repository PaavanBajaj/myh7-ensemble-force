# Continuous LSAR (`lsar-na-rel`) — method

First vertical slice. Modeled label: continuous LSAR
\(\mathrm{LSAR}=k_{\mathrm{cat,long-tail}}/k_{\mathrm{cat,short-tail}}\) on paired
two-headed human β-cardiac HMM constructs. Force-facing
\(N_{a,\mathrm{rel}}=\mathrm{LSAR}/0.57\) is recorded as a derived field and owned
in detail by the later `ensemble-force` method.

## Status (S1)

Provenance-ready public derivatives are published under
`data/public/lsar-na-rel/` with source identifiers in
`data/public/manifests/lsar-na-rel-sources.json`. No fitted model and no pinned
regeneration pipeline yet (S2+).

## Inclusion criteria

Exact continuous-LSAR candidates require:

- human β-cardiac myosin / MYH7
- the evaluated amino-acid variant
- paired short-tail and long-tail two-headed HMM from an accepted family
  (25-hep/2-hep or 15-hep/8-hep)
- actin-activated steady-state ATPase fitted \(k_{cat}\)
- exact numeric primary values (text, table, SI, or first-party numeric source data)

## Derived quantities (kept distinct)

1. Continuous LSAR (modeled biochemical label).
2. Spudich additional-\(N_a\) percent: \(100(\mathrm{LSAR}-0.57)/0.43\).
3. \(N_{a,\mathrm{rel}}=\mathrm{LSAR}/0.57\).
4. Within-study normalization:
   \(100(\mathrm{LSAR}_{mut}-\mathrm{LSAR}_{WT})/(1-\mathrm{LSAR}_{WT})\).
5. Table-implied LSAR from rounded Table 1 percent:
   \(0.57+0.43(p/100)\) — audit-only; never a training label.

Raw LSAR and within-study normalizations are not capped. Apparent Table 1
display caps (and near-WT coding of negatives as 0) are recorded only as display
fields.

## Evidence handling

- Preserve reported LSAR separately from recomputed long/short.
- One canonical label row per variant; biological replicates live only in
  `measurements.csv`.
- G256E: approximate relative/figure evidence; `exclude_approximate`.
- Q222K, E536D, V606M: unresolved unpublished; `exclude_unresolved`.
- D239N, R453C, G741R: assay absent, not zero; `exclude_blank`.
- G768R: later Pathak published epoch is canonical; Spudich 2024 unpublished
  epoch remains unresolved and separate.
- R403Q / R663H: exact pairs with low biological \(n\) (`review_low_n`).
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
