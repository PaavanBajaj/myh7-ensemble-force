# Public LSAR / \(N_{a,\mathrm{rel}}\) release (`lsar-na-rel`)

S1 provenance-ready scrubbed derivatives for continuous LSAR on the 23 MYH7
variants compiled in Spudich et al. 2024 Table 1 (pink additional-\(N_a\) by
LSAR column), plus within-study \(N_{a,\mathrm{rel}}\) labels for the 16 validated
rows.

Research use only. Not a diagnostic test. Not clinical advice.

## What is published here

| File | Role |
| --- | --- |
| `na-rel-primary.csv` | **Primary lsar_na modeling table** (S1 data artifact; not an S2 regeneration pipeline): `variant`, `LSAR_mut`, `LSAR_WT`, `Na_rel`, `ln(Na_rel)` |
| `schema.json` | Field definitions, formulas, and eligibility vocabulary |
| `measurements.csv` | Source-level / biological-replicate evidence |
| `labels.csv` | One canonical row per variant (uncertainties + audit transforms) |
| `exclusions.csv` | Approximate, unresolved, assay-absent, and unresolved-prior-epoch rows |
| `SHA256SUMS` | SHA-256 digests of the derivative files above |
| `../manifests/lsar-na-rel-sources.json` | Public source identifiers (DOI / PMID / PMCID) and derivative checksums |

## Scientific summary

- **Modeled biochemical label:** \(\mathrm{LSAR}=k_{\mathrm{cat,long}}/k_{\mathrm{cat,short}}\).
- **Canonical \(N_{a,\mathrm{rel}}\):** within-study
  \(\mathrm{LSAR}_{mut}/\mathrm{LSAR}_{WT,\mathrm{study}}\).
- **GPR target:** \(y=\ln(N_{a,\mathrm{rel}})\).
- **Kept separate:** Spudich additional-\(N_a\) percent
  \(100(\mathrm{LSAR}-0.57)/0.43\), global \(\mathrm{LSAR}/0.57\)
  (`na_rel_global_spudich057`), and release
  \(a=(\mathrm{LSAR}_{mut}-\mathrm{LSAR}_{WT})/(1-\mathrm{LSAR}_{WT})\)
  (audit-only).
- **Primary pack rows:** 13 Primary + 3 Risky allowed (R403Q, R663H, G768R).
  G256E is excluded from `na-rel-primary.csv`.
- **Reported vs recomputed:** author-reported LSAR is preserved separately from
  long/short recomputed from printed means.
- **No capping** of raw LSAR or within-study transforms; Table 1 display caps
  are recorded only as display fields.
- **G256E** is approximate and not exact-label modeling-ready.
- **Q222K, E536D, V606M** remain unresolved unpublished.
- **D239N, R453C, G741R** are assay-absent blanks, not zeroes.
- **G768R** canonical label uses the later Pathak published epoch; the Spudich
  2024 unpublished epoch stays unresolved and separate.
- Preprint/final versions of the same summary are one evidence identity.

## Stage claim

This pack is **S1 (provenance-ready)**. It is not S2: there is no pinned,
runnable regeneration pipeline in this release.

## Licensing limitation

No separate public data license is declared for this compiled derivative.
Source articles retain their publishers’ copyrights. The repository MIT license
covers repository software and documentation text, not third-party paper
copyright or redistributable full text. Redistribute only these scrubbed numeric
derivatives and public identifiers under cite-not-contain rules.

## Cite, do not contain

This directory does not include literature PDFs, screenshots, copyrighted full
text, vault paths, or temporary research notes. Locate primary evidence from the
stable identifiers in `../manifests/lsar-na-rel-sources.json`.
