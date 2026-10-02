# MYH7 variant biochemistry and mechanism

This repository publishes the reproducible LSAR / within-study
$N_{a,\mathrm{rel}}$ workflow from an earlier effort to predict relative
ensemble force for MYH7 variants. The current research direction is to examine
measured biochemical effects and evidence for variant mechanisms. The public
LSAR analysis is a small, descriptive input to that work, not a validated
predictor of head-state occupancy or force.

The original plan combined predicted ATPase turnover, unloaded velocity, and
available myosin heads into one force estimate. It did not reach that stage:
velocity modeling, forward prediction, and the combined calculation were not
implemented. A separate, private 17-variant ATPase diagnostic also failed its
predeclared prediction gates. Its data and results are not published here. The
scope and validation plan for new mechanism analyses are still being defined.

Research use only. This is not a diagnostic test and not clinical advice.

## Lineage

Initial project idea archived under `hcm-mava-response` (now private). That
workflow aimed to predict mavacamten treatment response in HCM patients with
machine learning. It was discontinued because usable labeled data were severely
insufficient. This repository is a clean start for the MYH7 relative
ensemble force workflow, which has since become a historical research direction.
It does not continue that Git history and does not redistribute literature PDFs
or other restricted content.

## Current stage

**S3 LSAR results shipped.** The continuous LSAR / within-study
$N_{a,\mathrm{rel}}$ workflow is runnable from repository-local inputs. It
includes outcome-free feature generation, 110 tests, an `osx-arm64` explicit
environment lock, the registered three-run LOSO ladder, and regenerated result
packs under `workflows/lsar-na-rel/results/`. Its best pooled out-of-fold MAE
was 0.130886 versus 0.131202 for the training-fold mean baseline. That small
difference does not establish useful prediction, statistical significance, or a
mutation mechanism. S3 describes release reproducibility, not model success.

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
