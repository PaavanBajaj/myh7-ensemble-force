# Workflows

Thin index of **live** public workflow trees. The original predictive plan's
later label and ensemble-force trees were not implemented.

| Workflow guide and results | Code | Tests | Public data | Status |
| --- | --- | --- | --- | --- |
| [`lsar-na-rel`](../workflows/lsar-na-rel/) | [`src/lsar_na_rel`](../src/lsar_na_rel/) | [`tests/lsar_na_rel`](../tests/lsar_na_rel/) | [`data/public/lsar-na-rel`](../data/public/lsar-na-rel/) | S3 results shipped; three registered LOSO runs |
| [`kcat-rel`](../workflows/kcat-rel/) | [`src/kcat_rel`](../src/kcat_rel/) | [`tests/kcat_rel`](../tests/kcat_rel/) | [`data/public/kcat-rel`](../data/public/kcat-rel/) | 17-variant evidence and frozen negative result published |
| [`myh7-correlation-screen`](../workflows/myh7-correlation-screen/) | [`src/myh7_correlation_screen`](../src/myh7_correlation_screen/) | [`tests/myh7_correlation_screen`](../tests/myh7_correlation_screen/) | [`data/public/myh7-correlation-screen`](../data/public/myh7-correlation-screen/) | 78 comparisons, point data, plots, and sensitivities published |

Each `workflows/` directory holds its runbook, configuration, and frozen run
records. `src/` holds executable code, `tests/` holds command and contract
checks, and `data/public/` holds shareable inputs and generated tables.

Public pack: [`data/public/lsar-na-rel/`](../data/public/lsar-na-rel/) ·
primary pipeline table: [`na-rel-primary.csv`](../data/public/lsar-na-rel/na-rel-primary.csv) ·
manifest: [`data/public/manifests/lsar-na-rel-sources.json`](../data/public/manifests/lsar-na-rel-sources.json)

Runbook: [`workflows/lsar-na-rel/README.md`](../workflows/lsar-na-rel/README.md) ·
results: [`workflows/lsar-na-rel/results/`](../workflows/lsar-na-rel/results/)

The [correlation decision](myh7-correlation-screen-decision.md) selects a
[matched-preparation validation question](matched-atpase-velocity-protocol.md).
There is no velocity predictor or combined-force workflow.
