# Public data

Scrubbed, redistribute-safe derivatives and manifests only.

- Never place literature PDFs, full-text extracts, or vault binaries here.
- Cite-not-contain: evidence stays in the private vault; this tree holds pointers
  and shareable scrubbed tables only.

## Current packs

| Pack | Stage | Path | Manifest |
| --- | --- | --- | --- |
| Continuous LSAR / within-study $N_{a,\mathrm{rel}}$ (`lsar-na-rel`) | S3 | [`lsar-na-rel/`](lsar-na-rel/) (primary: [`na-rel-primary.csv`](lsar-na-rel/na-rel-primary.csv), features: [`features.csv`](lsar-na-rel/features.csv)) | [`manifests/lsar-na-rel-sources.json`](manifests/lsar-na-rel-sources.json) |
| Actin-activated ATPase ratio (`kcat-rel`) | Frozen diagnostic | [`kcat-rel/`](kcat-rel/) (labels: [`canonical-labels.csv`](kcat-rel/canonical-labels.csv), features: [`derived/kcat-rel-v0-no-msa-features.csv`](kcat-rel/derived/kcat-rel-v0-no-msa-features.csv)) | [`derived/kcat-rel-v0-no-msa-feature-manifest.json`](kcat-rel/derived/kcat-rel-v0-no-msa-feature-manifest.json) |
| MYH7 correlation screen | Descriptive analysis | [`myh7-correlation-screen/`](myh7-correlation-screen/) (results: [`comparison-results-v1.csv`](../../results/myh7-correlation-screen/results/screen-v1/comparison-results-v1.csv), ranking: [`spearman-rankings-no-ihm-v1.csv`](../../results/myh7-correlation-screen/results/screen-v1/spearman-rankings-no-ihm-v1.csv)) | [`analysis-manifest-v1.json`](../../results/myh7-correlation-screen/results/screen-v1/analysis-manifest-v1.json) |
| MYH7 nonlinear dependence | Descriptive extension | [`myh7-nonlinear-dependence/`](myh7-nonlinear-dependence/) (reuses published correlation inputs; [combined CSV](../../results/myh7-nonlinear-dependence/results/screen-v1/combined-correlations-v1.csv)) | [`dependence-plan-v1.json`](myh7-nonlinear-dependence/dependence-plan-v1.json) |

Code for these packs lives in `../../src/`; tests live in `../../tests/`;
workflow guides and frozen result packs live in `../../results/`.
