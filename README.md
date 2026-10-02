# MYH7 variant biochemistry and mechanism

This repository publishes the completed LSAR / within-study
$N_{a,\mathrm{rel}}$ workflow, the frozen `kcat_rel` diagnostic, and the MYH7
correlation screen. The current research direction is to examine measured
biochemical effects and validate associations with matched preparations.
These analyses do not establish head-state occupancy, a force predictor, or a
mutation mechanism.

The original plan combined predicted ATPase turnover, unloaded velocity, and
available myosin heads into one force estimate. It did not reach that stage:
velocity modeling, forward prediction, and the combined calculation were not
implemented. The 17-variant ATPase diagnostic failed its predeclared prediction
gates. Its data, code, and negative result are now published. The subsequent
correlation screen compares measured ATPase, LSAR-derived available-head proxy,
unloaded velocity, and structural features without training a predictor.

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

**Completed public workflows: LSAR, `kcat_rel`, and correlation screen.** The continuous LSAR / within-study
$N_{a,\mathrm{rel}}$ workflow is runnable from repository-local inputs. It
includes outcome-free feature generation, 110 tests, an `osx-arm64` explicit
environment lock, the registered three-run LOSO ladder, and regenerated result
packs under `results/lsar-na-rel/results/`. Its best pooled out-of-fold MAE
was 0.130886 versus 0.131202 for the training-fold mean baseline. That small
difference does not establish useful prediction, statistical significance, or a
mutation mechanism. The `kcat_rel` Ridge and GPR diagnostics also failed their
predeclared gates. The exploratory screen reports all 78 eligible comparisons;
its strongest outcome relationship is `kcat_rel` versus unloaded velocity
(Spearman rho +0.670, 14 paired variants), which motivates matched-measurement
validation, not a causal claim.

## Start here

1. `environment.yml` and `environment-osx-arm64.lock`: direct and exact pins
2. `docs/reproducibility.md`: stage ladder and what "reproducible" means
3. `docs/workflows.md`: index of completed workflow trees
4. `results/kcat-rel/`: ATPase method, config, and frozen negative diagnostic
5. `results/myh7-correlation-screen/`: correlation method and runbook
6. `data/public/`: published inputs and derived tables; frozen result packs and plots are under `results/`

Executable code lives in `src/lsar_na_rel/`, `src/kcat_rel/`, and
`src/myh7_correlation_screen/`. Their tests live in the matching `tests/`
directories. Each `results/<workflow>/` contains a runbook, configuration,
and its frozen result pack. Downloaded structures use ignored `.cache/` paths.
Run commands from this repository's root.

## Documentation map

- `docs/decisions-and-limitations.md`: living major-decision and limitations log
- `docs/provenance.md`: cite-not-contain evidence account
- `docs/reproducibility.md`: stage contract
- `docs/workflows.md`: live slice index
- `data/public/`: scrubbed shareable inputs, derivatives, and manifests only
  (never vault contents)
