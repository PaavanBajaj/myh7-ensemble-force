# MYH7 nonlinear and joint dependence screen

This published descriptive extension follows the [approved spec](../../docs/myh7-nonlinear-dependence-spec.md) and [method review](../../docs/myh7-nonlinear-dependence/ml-alternatives.md). It adds conventional distance correlation and a signed normalized U-centered squared distance-correlation score to the existing 78 comparisons, then evaluates 12 fixed feature-block/outcome comparisons. It does not train a predictor or calculate actual-cohort p-values.

The main deliverable is [combined-correlations-v1.csv](results/screen-v1/combined-correlations-v1.csv): 90 rows with historical effects, new scores, coverage, sensitivity ranges and three plain-language sentences per row. The [data dictionary](data-dictionary.md) explains the columns. [dependence-sensitivities-v1.csv](results/screen-v1/dependence-sensitivities-v1.csv) lists every variant and source-study deletion plus the exact shared-component, strict velocity and Morck-control exclusions. [dependence-points-v1.csv](results/screen-v1/dependence-points-v1.csv) exposes raw input values for independent checks.

## Fixed analysis choices

The [configuration](config/dependence-config-v1.json) defines three primary blocks before calculating these scores: chemistry (charge change and Grantham distance), ADP geometry (ADP-to-rigor shift and ADP adenine proximity), and motor geometry (actin, relay landmark and essential-light-chain proximity). Each is evaluated against relative catalytic turnover, relative available heads and primary unloaded velocity, producing nine primary block rows. The all-ten-feature block contributes three sensitivity rows with reduced complete-case coverage. These groups are biological working hypotheses, not validated mechanistic modules; a block score cannot distinguish a single-feature effect from an interaction.

Complete cases are selected separately for each comparison, without imputation. Feature columns are centered and divided by their sample standard deviation (`ddof=1`), using features alone within that cohort; scaling is recomputed after every deletion. Constant dimensions are recorded and dropped. Outcomes retain raw mutant-to-matched-WT ratios. Rank sensitivity replaces each column with average ranks and rescales features; log sensitivity replaces only outcome ratios with natural logarithms. These sensitivities change the estimand and are never used to select the largest score.

Conventional distance correlation is nonnegative, uses double-centered Euclidean-distance matrices and is positively biased in small samples. The U-centered score uses off-diagonal U centering and the normalized inner product of the resulting matrices. Its squared covariance numerator is unbiased under the iid model; the estimated correlation ratio is **not** an unbiased estimator. It can be negative, is never square-rooted or clipped, and neither score indicates an increasing versus decreasing direction. Fewer than four complete cases, a constant side or zero U distance variance produce an explicit blank U score and status. A binary feature with one observation in one group has zero U distance variance; a relative floating-point tolerance prevents numerical residuals from producing an arbitrary ratio.

Study deletion removes rows belonging to the omitted study on either outcome; it is unavailable for feature-only pairs. C001 excludes seven exact reused mutant short-head ATPase components identified in the dependence ledger, with a separate union of confirmed mutant/WT reuse. The five Morck primary variants sharing slide/channel WT records are excluded from all rows involving primary velocity. Historical strict short-head MVEL20 subsets and the shared-component/uncertain-Na subset are preserved as sensitivities. These checks leave other study dependencies in place and are not independent replication. Expanded velocity and alternative metric analyses remain in the historical workflow and are not duplicated here.

## Reproduce

From the public repository root, activate the environment defined in `environment.yml` (named `MYH7 research`), then run:

```bash
PYTHONPATH=src python -m myh7_nonlinear_dependence.screen plan
PYTHONPATH=src python -m myh7_nonlinear_dependence.screen analyze
PYTHONPATH=src python -m pytest -q tests/myh7_nonlinear_dependence tests/myh7_correlation_screen
```

The committed [plan](../../data/public/myh7-nonlinear-dependence/dependence-plan-v1.json) uses the published correlation input pack and its historical results under `results/myh7-correlation-screen/results/screen-v1/`. No private repository, vault workbook, paper download or structure download is required. The estimator, orchestration and CLI are separate modules in [`src/myh7_nonlinear_dependence/`](../../src/myh7_nonlinear_dependence/).

For a fresh rerun using the committed plan:

```bash
PYTHONPATH=src python -m myh7_nonlinear_dependence.screen analyze \
  --plan data/public/myh7-nonlinear-dependence/dependence-plan-v1.json \
  --output-dir /tmp/myh7-nonlinear-rerun
```

Optional flags are `--source-dir`, `--config`, `--output-dir` and `--plan`. Default `plan` writes to the public input directory; default `analyze` writes to the frozen result directory. With `--output-dir`, both commands use that directory, and analysis defaults to its plan unless `--plan` is explicit. A custom `--source-dir` must contain all nine files listed in the [input guide](../../data/public/myh7-nonlinear-dependence/README.md); the default source uses the repository's split layout. The configuration is the fixed v1 contract, not a parameter-search interface.

`plan` validates historical identities, counts, paired variants, binary groups and effects before freezing nine input hashes, the implementation/configuration hashes, all 78 historical IDs and 12 block IDs. It estimates no new distance scores. A content digest detects plan tampering; `analyze` additionally reconstructs and compares the entire plan before calculating or writing any results. The saved plan is a commitment for reproducibility, not protection against someone intentionally rewriting both a plan and its digest. This extension is post hoc relative to previously viewed Spearman results, not a prospective preregistration.

For fixed inputs, code and package versions, all CSVs and manifests are deterministic, with no timestamps or absolute local paths. The [run manifest](results/screen-v1/dependence-run-manifest-v1.json) records input, code, config, plan and output hashes plus package versions. Changing the implementation requires a new plan; rerun `plan` explicitly and document why. A different numerical runtime may change final floating-point digits.

## Interpretation

The existing relative catalytic turnover/velocity lead has rho +0.6703, conventional distance correlation 0.7610 and U squared score 0.4001 in the same 14 variants. These are different statistics and should not be compared as interchangeable effect sizes. The fixed primary block U scores range from −0.1656 to +0.0690, while the all-feature sensitivity scores are negative in cohorts of 11, 11 and 9 variants. These results do not establish a new joint biological relationship; a near-zero or negative finite-sample score does not establish independence.

Population distance covariance characterizes general dependence under its mathematical assumptions, but these few heterogeneous literature variants do not supply independent, exchangeable measurements. Source controls, literal ATPase reuse, static WT geometry and missingness all limit interpretation. Deletion ranges are influence checks, not confidence intervals. Claims about significance, causality, mechanisms, interactions or prediction require a defensible independent measurement and validation design. The scientific next step remains the existing matched-assay validation protocol.

The [publication audit](../../docs/nonlinear-dependence-publication-audit-2026-10-02.md) records the public rerun and comparison against the independently verified development results. See the [current decision](../../docs/myh7-nonlinear-dependence-decision.md) for the scientific follow-up.
