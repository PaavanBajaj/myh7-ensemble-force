# Predicting Hypertrophic Cardiomyopathy gene MYH7 relative ensemble force in 30+ variants

This project estimates how much ensemble force HCM-linked MYH7 variants may
change in the heart muscle motor, relative to normal MYH7. It does that by
modeling a few measurable biochemical labels separately, then combining those
predictions into one relative ensemble force estimate (with uncertainty).

We start with one label, continuous LSAR. Later labels are relative ATPase rate
(kcat_rel) and relative velocity (v_rel). After those three labels exist, we
combine them into the relative ensemble force estimate. Each model also reports
uncertainty, so a guess is never shown as more certain than the evidence supports.

Research use only. This is not a diagnostic test and not clinical advice.

## Lineage

Initial project idea archived under `hcm-mava-response` (now private). That
workflow aimed to predict mavacamten treatment response in HCM patients with
machine learning. It was discontinued because usable labeled data were severely
insufficient. This repository is a clean start for the MYH7 relative
ensemble force workflow. It does not continue that Git history and does not
redistribute literature PDFs or other restricted content.

## Current stage

**S3 Label results shipped.** The continuous LSAR / within-study
$N_{a,\mathrm{rel}}$ workflow is runnable from repository-local inputs. It
includes outcome-free feature generation, 110 tests, an `osx-arm64` explicit
environment lock, the registered three-run LOSO ladder, and regenerated result
packs under `workflows/lsar-na-rel/results/`.

## Start here

1. `environment.yml` and `environment-osx-arm64.lock`: direct and exact pins
2. `docs/reproducibility.md`: stage ladder and what "reproducible" means
3. `docs/workflows.md`: index of live workflow trees
4. `data/public/lsar-na-rel/`: labels, generated features, provenance, and audits
5. `workflows/lsar-na-rel/`: commands, configuration, method, and results

## Documentation map

- `docs/decisions-and-limitations.md`: living major-decision and limitations log
- `docs/provenance.md`: cite-not-contain evidence account
- `docs/reproducibility.md`: stage contract
- `docs/workflows.md`: live slice index
- `data/public/`: scrubbed shareable derivatives and manifests only (never vault contents)
