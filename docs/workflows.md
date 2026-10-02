# Workflows

Thin index of **live** public workflow trees. The original predictive plan's
later label and ensemble-force trees were not implemented.

| Workflow guide and results | Code | Tests | Public data | Status |
| --- | --- | --- | --- | --- |
| [`lsar-na-rel`](../results/lsar-na-rel/) | [`src/lsar_na_rel`](../src/lsar_na_rel/) | [`tests/lsar_na_rel`](../tests/lsar_na_rel/) | [`data/public/lsar-na-rel`](../data/public/lsar-na-rel/) | S3 results shipped; three registered LOSO runs |
| [`kcat-rel`](../results/kcat-rel/) | [`src/kcat_rel`](../src/kcat_rel/) | [`tests/kcat_rel`](../tests/kcat_rel/) | [`data/public/kcat-rel`](../data/public/kcat-rel/) | 17-variant evidence and frozen negative result published |
| [`myh7-correlation-screen`](../results/myh7-correlation-screen/) | [`src/myh7_correlation_screen`](../src/myh7_correlation_screen/) | [`tests/myh7_correlation_screen`](../tests/myh7_correlation_screen/) | [`data/public/myh7-correlation-screen`](../data/public/myh7-correlation-screen/) | 78 comparisons, point data, plots, and sensitivities published |

Each `results/<workflow>/` directory holds its README runbook, `config/`, and
frozen `results/`. `src/` holds executable code, `tests/` holds command and
contract checks, `data/public/` holds shareable inputs and intermediate tables,
and `docs/` holds method and evidence accounts. Downloaded structures use the
ignored repository-level `.cache/` directory.

Public pack: [`data/public/lsar-na-rel/`](../data/public/lsar-na-rel/) ·
primary pipeline table: [`na-rel-primary.csv`](../data/public/lsar-na-rel/na-rel-primary.csv) ·
manifest: [`data/public/manifests/lsar-na-rel-sources.json`](../data/public/manifests/lsar-na-rel-sources.json)

Runbook: [`results/lsar-na-rel/README.md`](../results/lsar-na-rel/README.md) ·
results: [`results/lsar-na-rel/results/`](../results/lsar-na-rel/results/)

The [correlation decision](myh7-correlation-screen-decision.md) selects a
[matched-preparation validation question](matched-atpase-velocity-protocol.md).
There is no velocity predictor or combined-force workflow.
