# Nonlinear dependence publication audit

Date: 2026-10-02. This records the fourth completed public workflow, promoted after the independently verified development run.

## Public layout and provenance

- Estimators: `src/myh7_nonlinear_dependence/statistics.py`.
- Frozen comparison orchestration: `src/myh7_nonlinear_dependence/analysis.py`.
- Command entry point: `src/myh7_nonlinear_dependence/screen.py`.
- Tests: `tests/myh7_nonlinear_dependence/`.
- Configuration: `results/myh7-nonlinear-dependence/config/dependence-config-v1.json`.
- Public plan: `data/public/myh7-nonlinear-dependence/dependence-plan-v1.json`.
- Frozen outputs: `results/myh7-nonlinear-dependence/results/screen-v1/`.
- Method reviews and decision: `docs/myh7-nonlinear-dependence/` and the linked spec/decision documents.

The public plan was regenerated against seven published correlation input files and two historical result files. It hashes the scrubbed public outcomes and public historical plan, the unchanged fixed configuration, and all four executable package modules. Public provenance therefore intentionally differs from the development plan and manifest. The three numerical CSVs are unchanged.

No raw workbooks, paper PDFs, full-text extracts, email content, credentials or absolute machine paths were added. No new data license is declared. Earlier diagnostic and correlation result packs were preserved. Living repository indexes, runbooks, provenance, reproducibility and next-step documents were updated; the earlier publication audit remains a historical record.

## Numerical verification

The original development audit independently recomputed 90 comparisons and 1,761 sensitivities with separate loop formulas, including high-precision checks for mathematically degenerate binary U variance. Public packaging retains those formulas and scores. The public rerun reproduces the development CSVs byte for byte:

- `combined-correlations-v1.csv`: `95c42973e3ccd2a5316e80da8408bd651f767219ebf402a687ffed78f3b77d02`.
- `dependence-sensitivities-v1.csv`: `c7c76e272b5d9c9330ff1a0c3677b92423452f98ee3ea8c359e8704b29ed1032`.
- `dependence-points-v1.csv`: `da218860dd72baa5b3a528ac4c0a3b305afed65a59612cd051b5ba51069465f7`.

A separate temporary checkout containing only public package files, public input tables, historical public results and the new configuration/plan reproduced all five public artifacts byte for byte: three CSVs, the plan and the run manifest. It required no private repository, vault evidence or network download.

## Tests and runtime

The complete public repository suite passed: **630 tests**, with 37 existing Gaussian-process convergence warnings. The new package contributes 12 cases, including independent estimator fixtures, preserved negative scores, singleton binary degeneracy, curved and joint synthetic relationships, frozen original pair identities, exact shared-component exclusions, study-union deletions, invalid input rejection, plan tampering and deterministic CLI artifacts. A committed-plan test runs from outside the repository and reproduces all four result files.

The existing structure tests initially failed because DSSP could not find its executable and installed chemistry dictionaries in the relocated environment. With the scientific environment's `bin` on `PATH` and `LIBCIFPP_DATA_DIR` pointing to its `share/libcifpp`, all 32 affected-file tests passed, followed by the complete 630-test pass. No structure code or frozen features changed. This runtime lookup requirement is documented in [reproducibility](reproducibility.md).

Runtime: Python 3.12.13, NumPy 2.5.2 and SciPy 1.18.0, matching the existing public dependency pins. Different floating-point runtimes may change final digits; manifests record the numerical versions.

## Scientific decision

The [result decision](myh7-nonlinear-dependence-decision.md) remains descriptive. ATPase versus velocity retains rho +0.6703 in 14 variants, conventional distance correlation 0.7610 and normalized U-centered squared score 0.4001. The nine compact primary feature blocks provide no stable new lead under source-study deletion. Negative U estimates do not mean inverse biology or establish independence. No cohort p-values, confidence intervals, interaction attribution, predictive accuracy or mechanism claims are added. Independent matched-preparation ATPase/velocity measurements remain the next scientific step.
