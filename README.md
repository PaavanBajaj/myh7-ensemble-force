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

**S0 Scaffold.** Clone this repo and read the intent. There are no fitted
models yet, no public label tables, and no results packs.

## Start here

1. `environment.yml`: env pin stub (fill before any results pack)
2. `docs/reproducibility.md`: stage ladder and what "reproducible" means
3. `docs/workflows.md`: index of live workflow trees
4. `workflows/lsar-na-rel/`: first vertical slice (stubs only at S0)

## Documentation map

- `docs/decisions-and-limitations.md`: living major-decision and limitations log
- `docs/provenance.md`: cite-not-contain evidence account
- `docs/reproducibility.md`: stage contract
- `docs/workflows.md`: live slice index
- `data/public/`: scrubbed shareable derivatives and manifests only (never vault contents)
