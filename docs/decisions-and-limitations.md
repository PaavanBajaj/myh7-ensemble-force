# Decisions and limitations

Living log for promoted major methodology decisions and standing limitations.
No ADR folder for v1 — updates land here after explicit promotion OK.

## Status (S1)

Promoted decisions for the continuous LSAR public pack
(`data/public/lsar-na-rel/`):

1. Modeled biochemical label is continuous LSAR
   (\(k_{\mathrm{cat,long}}/k_{\mathrm{cat,short}}\)), not the Spudich display
   percent and not \(N_{a,\mathrm{rel}}\).
2. Spudich additional-\(N_a\), \(N_{a,\mathrm{rel}}\), within-study
   normalization, and Table-implied reverse LSAR are stored as distinct derived
   fields; reverse LSAR is audit-only.
3. One canonical label row per variant; replicates belong only in measurements.
4. Reported LSAR is preserved separately from recomputed long/short; neither is
   silently replaced to force Table 1 agreement.
5. Raw LSAR and within-study labels are not capped; display caps are recorded
   separately.
6. G256E remains approximate and excluded from exact-label modeling.
7. Q222K, E536D, and V606M remain unresolved unpublished.
8. D239N, R453C, and G741R are assay-absent blanks, not zeroes.
9. G768R’s later Pathak measurement is canonical and separate from Spudich’s
   unpublished 2024 epoch.
10. R403Q / R663H remain low-biological-\(n\) evidence (`review_low_n`).
11. Preprint/final versions and repeated presentations of the same summary are
    one evidence identity.
12. No new outside-Table-1 exact LSAR variant was recovered in the audited
    expansion screen summarized with this pack.
13. No separate public data license is declared for the compiled derivative
    beyond the repository MIT license’s ordinary scope (see pack README).

## Limitations (standing)

- Research-use only; not a diagnostic test; not clinical advice.
- Vault evidence is cite-not-contain — literature PDFs and full-text extracts
  do not live in this repository.
- S1 is provenance-ready only; regenerating the pack from a pinned public
  pipeline requires S2.
- Construct family, assay chemistry, study WT baseline, and low-\(n\) rows
  still need pooling review before any fitted model.
- Later label trees (`kcat-rel`, `v-rel`, `ensemble-force`) are created only
  when those slices start.
