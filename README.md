# MYH7 relative ensemble force

Uncertainty-aware relative ensemble-force estimates for HCM-associated MYH7
variants from separate label models (continuous LSAR first; later kcat_rel,
v_rel, and the ensemble combination).

Research-use only. Not a diagnostic test. Not clinical advice.

## Lineage

Initial project idea archived under `hcm-mava-response` (now private). That
workflow aimed to predict mavacamten treatment response in HCM patients with
machine learning; it was discontinued because usable labeled data were severely
insufficient. This repository is a clean start for the MYH7 relative
ensemble-force workflow. It does not continue that Git history and does not
redistribute literature PDFs or other restricted content.

## Current stage

**S0 — Scaffold.** Clone and read intent. No fitted models, no public label
tables, and no results packs yet.

## Start here

1. `environment.yml` — env pin stub (fill before any results pack)
2. `docs/reproducibility.md` — stage ladder and what “reproducible” means
3. `docs/workflows.md` — index of live workflow trees
4. `workflows/lsar-na-rel/` — first vertical slice (stubs only at S0)

## Documentation map

- `docs/decisions-and-limitations.md` — living major-decision + limitations log
- `docs/provenance.md` — cite-not-contain evidence account
- `docs/reproducibility.md` — stage contract
- `docs/workflows.md` — live slice index
- `data/public/` — scrubbed shareable derivatives and manifests only (never vault contents)
