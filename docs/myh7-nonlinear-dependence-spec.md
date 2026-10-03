# MYH7 nonlinear and multivariate dependence screen

Status: implemented and publicly promoted. Original spec and test boundaries auto-approved by Paavan on 2026-10-02.

## Problem Statement

The completed MYH7 screen describes monotonic pairwise relationships, but a curved relationship or a relationship involving several features together can escape Spearman rho. Paavan needs a reproducible way to look for these patterns in the existing small, heterogeneous measured-data cohorts and a compact CSV that can be discussed with his mentors.

## Solution

Add a descriptive distance-correlation screen alongside the frozen Spearman results. Report conventional distance correlation and a normalized U-centered, bias-corrected squared distance-correlation score for all existing comparisons and for fixed biological feature blocks against each canonical outcome. Retain negative U scores, explicit missingness, study and shared-measurement sensitivities, and simple row descriptions. These are exploratory dependence scores, not a trained ML predictor, inferential tests, or evidence of a mechanism.

## User Stories

1. As a researcher, I want the literature recommendation grounded in my actual cohort, so that I choose a method suitable for 15–17 variants.
2. As a researcher, I want nonlinear pairwise scores beside Spearman rho, so that curved patterns can be considered.
3. As a researcher, I want joint feature-block scores, so that multivariate relationships can be explored.
4. As a researcher, I want every previous comparison retained, so that weak results are visible.
5. As a researcher, I want binary comparisons labeled correctly, so that their median difference is never mistaken for rho.
6. As a researcher, I want the original sign of Spearman retained, so that association direction remains clear.
7. As a researcher, I want negative U-centered scores retained, so that small-sample noise is not concealed.
8. As a researcher, I want complete-case counts and variant identities, so that different coverage is visible.
9. As a researcher, I want biological feature blocks fixed before new estimates, so that outcome-driven selection is avoided.
10. As a researcher, I want features scaled without using labels, so that units do not dominate multivariate distances.
11. As a researcher, I want leave-one-variant checks, so that outlier influence is exposed.
12. As a researcher, I want leave-one-source-study checks, so that study-specific patterns are exposed.
13. As a researcher, I want exact shared-short-head ATPase exclusions, so that algebraic coupling is visible.
14. As a researcher, I want the shared Morck velocity-control exclusion checked, so that reused controls are acknowledged.
15. As a researcher, I want unestimable cases stated explicitly, so that missing scores are never treated as zero.
16. As a researcher, I want a frozen input manifest, so that a rerun detects changed evidence.
17. As a researcher, I want deterministic result files, so that I can compare independent reruns.
18. As a researcher, I want two or three plain-language sentences per CSV row, so that my mentors can interpret the results quickly.
19. As a researcher, I want a distinction between pairwise and block rows, so that joint scores are not mistaken for feature attribution.
20. As a researcher, I want the all-feature comparison marked as sensitivity, so that its reduced coverage is clear.
21. As a researcher, I want documented statistical limits, so that exploratory scores do not become claims of significance or prediction.
22. As a maintainer, I want independent execution and statistical verification, so that implementation errors are caught before review.
23. As a maintainer, I want a separate branch and PR, so that the work is reviewable and reversible.

## Implementation Decisions

- Use Euclidean-distance covariance with conventional double centering and U centering. The normalized U-centered squared score is bias corrected, but its ratio is not an unbiased estimator. It can be negative and has no direction of association.
- Reuse the existing frozen comparison inventory, analysis table, results, and dependence ledgers. Verify original pair IDs, complete-case variants and counts. Preserve all 78 original comparisons.
- Freeze three primary blocks by biological rationale: chemistry (charge change and Grantham distance), ADP geometry (ADP-to-rigor shift and adenine proximity), and motor geometry (actin, relay-landmark and ELC proximity). Evaluate each against relative catalytic turnover, relative available heads, and primary unloaded velocity. An all-ten-feature block is a separate dimensionality and coverage sensitivity. Do not choose blocks by new scores.
- Use complete cases without imputation. Standardize each varying feature column to unit sample standard deviation within the evaluated cohort. Constant dimensions are explicitly identified; a wholly constant side or fewer than four usable rows is unestimable for U scores. Recompute scaling for deletion sensitivities.
- Preserve raw ratio outcome scale. Include rank-distance sensitivity to distinguish scale/outlier sensitivity from Spearman's monotonic summary.
- Record individual variant deletion and source-study deletion scores. For outcome pairs, remove a study's rows if that study appears on either outcome. Feature-only pairs have no source-study deletion claim.
- Identify exact reused mutant and WT short-head ATPase components from the dependence ledger and rerun the catalytic-turnover/available-heads comparison after their removal. Rerun comparisons involving primary velocity after removing Morck shared-slide variants. These exclusions are sensitivity analyses, not independent replication.
- Do not calculate actual-cohort p-values, confidence intervals, multiplicity-adjusted discoveries, or model accuracy. Neither pooled nor within-study exchangeability has been established. Document how conditional blocked permutation inference would require a defensible design and a prespecified multiplicity family in a future matched dataset.
- Provide plan and analyze commands. The plan contains input hashes, method settings, feature groups and comparison IDs before new estimates. Analyze validates the plan before writing deterministic CSVs and a run manifest with input/output hashes and runtime versions.
- The combined CSV contains original method/effect, rho only where originally computed, new scores, sample coverage, sensitivity summaries, and a description field of two or three short sentences. Block rows have no invented Spearman coefficient.
- Preserve all frozen diagnostics, public sibling artifacts and unrelated existing research files. Developed first in the private research repository; public promotion was subsequently authorized by Paavan.
- GitHub Issues is the existing project tracker. Publish this spec with the requested ready-for-agent label. No unrelated agent-configuration rewrite is needed because the tracker is already identifiable from the repository's active issues.

## Testing Decisions

- Prefer the existing high-level command-line artifact seam: fixed input tables plus plan produce validated result CSVs. Tests assert observable outputs, rejection of invalid inputs, and reproducibility rather than internal helper structure.
- Use independent known estimator fixtures/reference formulas for conventional and U-centered distance scores, including symmetry, units, constant inputs and negative estimates.
- Verify a synthetic U-shaped relationship that Spearman misses and a synthetic joint interaction that weak marginal associations miss. These validate capability rather than real-data significance.
- Verify plan hash tampering, duplicated variants, nonfinite values, missingness, binary rows, block rows, invalid counts and unchanged historical results.
- Execute focused tests and a full real-data run; compare rerun CSV hashes. Independently verify score calculations and the shared-measurement exclusions. Reuse the correlation screen's command-line test conventions and installed scientific environment.

## Out of Scope

New evidence extraction, molecular dynamics, mutant structural acquisition, clinical claims, forward prediction, combined ensemble force, broad model tournaments, and discovery claims from p-values. The initial specification excluded public promotion, which Paavan subsequently authorized as a separate release step. Sending messages is outside the public workflow.

## Further Notes

This is an exploratory extension after the existing Spearman results were already seen, not a prospective biological preregistration. The fixed plan prevents further score-driven choices. Population distance covariance detects general dependence under its mathematical assumptions; a finite, correlated literature sample cannot establish independence from a near-zero score. A block association cannot identify an interaction, causal feature, blocked-head adoption or PPS stability. Stronger inference requires independent matched experimental preparations and validation across studies.

## Public packaging

Promoted to the public layout on 2026-10-02. Code is in `src/myh7_nonlinear_dependence/`, tests in `tests/myh7_nonlinear_dependence/`, configuration and frozen CSVs in `results/myh7-nonlinear-dependence/`, and the regenerated public plan in `data/public/myh7-nonlinear-dependence/`. The [runbook](../results/myh7-nonlinear-dependence/README.md) governs public execution paths; historical development paths above describe the original specification only.
