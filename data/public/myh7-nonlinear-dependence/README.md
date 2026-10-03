# Nonlinear dependence input plan

This pack adds a frozen public plan while reusing the existing public correlation derivatives. Numerical inputs are not duplicated. The [runbook](../../../results/myh7-nonlinear-dependence/README.md) explains regeneration; output CSVs and the run manifest are in [the frozen result pack](../../../results/myh7-nonlinear-dependence/results/screen-v1/).

Seven inputs live in `data/public/myh7-correlation-screen/`: `analysis-table-v1.csv`, `comparison-inventory-v1.csv`, `comparison-plan-manifest-v1.json`, `outcomes-v1.csv`, `dependence-ledger-v1.csv`, `velocity-dependence-v1.csv`, and `feature-coverage-v1.csv`. Two inputs live in `results/myh7-correlation-screen/results/screen-v1/`: `comparison-results-v1.csv` and `sensitivity-results-v1.csv`. A custom input directory must contain all nine files.

The [plan](dependence-plan-v1.json) hashes those public derivatives, the fixed configuration and all four executable package files. It was regenerated after packaging; the private development plan is not reused. The plan declares 78 historical pairs, nine compact primary feature-block comparisons and three all-feature sensitivities. This extension is post hoc relative to the Spearman screen, not a prospective preregistration.

Only source-linked numerical derivatives are redistributed. Literature, original workbooks and restricted evidence stay outside this tree. No additional data license is declared beyond the repository license's ordinary scope.
