# MYH7 Variant Biochemistry

This repository shares four completed MYH7 research workflows: an LSAR-derived
proxy for available myosin heads, a `kcat_rel` ATPase diagnostic, a screen
of relationships among ATPase, LSAR, unloaded velocity, and structural
features, and a nonlinear and joint dependence extension. The current focus
is to describe measured effects and check promising
relationships with matched experiments. These results do not establish a
mutation mechanism or support clinical decisions.

## Lineage

The project first aimed to estimate MYH7 ensemble force by combining predicted
ATPase turnover, unloaded velocity, and available myosin heads. The ATPase
diagnostic failed its prediction gates, and the velocity model and combined
force calculation were never built. That ensemble force plan is now historical.

Initial project idea archived under `hcm-mava-response` (now private). That
workflow aimed to predict mavacamten treatment response in HCM patients with
machine learning. It was discontinued because usable labeled data were severely
insufficient. This repository is a clean start for the MYH7 relative
ensemble force workflow, which has since become a historical research direction.
It does not continue that Git history and does not redistribute literature PDFs
or other restricted content.

## Current stage

The LSAR workflow is runnable, and the completed `kcat_rel` diagnostic is
preserved as a frozen record; neither produced a useful predictor. The
correlation screen reports all 78 eligible comparisons; its strongest outcome
association is between `kcat_rel` and
unloaded velocity (Spearman rho +0.670, 14 paired variants). The next step is to
test that relationship with matched measurements before drawing broader
conclusions. The nonlinear extension publishes 90 comparisons and 1,761
sensitivity rows. Its fixed feature blocks do not provide a stable new lead;
ATPase versus velocity remains the follow-up priority.

## Start here

1. `environment.yml` and `environment-osx-arm64.lock`: direct and exact pins
2. `docs/reproducibility.md`: stage ladder and what "reproducible" means
3. `docs/workflows.md`: index of completed workflow trees
4. `results/kcat-rel/`: ATPase method, config, and frozen negative diagnostic
5. `results/myh7-correlation-screen/`: correlation method and runbook
6. `results/myh7-nonlinear-dependence/`: joint feature-block scores and combined CSV
7. `data/public/`: published inputs and derived tables; frozen result packs and plots are under `results/`

Executable code lives in `src/lsar_na_rel/`, `src/kcat_rel/`,
`src/myh7_correlation_screen/`, and `src/myh7_nonlinear_dependence/`.
Their tests live in the matching `tests/`
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

## Additional Notes

This repository serves as the public research home for my MYH7 project. More detailed research decisions, experimental code, intermediate analyses, and other work-in-progress materials are developed and iterated on in the private `research-workhorse` repository before being deliberately refined and published here.
