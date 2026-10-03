# MYH7 dependence beyond Spearman: method and result decision

Date: 2026-10-02. Scope: descriptive extension now published in the public research home; earlier screens and prediction diagnostics remain frozen historical evidence.

## Decision

Use conventional distance correlation and a normalized U-centered, bias-corrected squared distance-correlation score as the first nonlinear/multivariate extension. Retain all 78 historical pairs and add nine compact primary feature-block comparisons plus three all-feature sensitivities. This is the most useful next analysis for the present dataset because it requires few choices, uses the existing data, and addresses general dependence without fitting another predictor. It is not a universal ranking of statistical methods.

The [primary-source research synthesis](myh7-nonlinear-dependence/README.md) explains the alternatives and limitations. The [approved spec](myh7-nonlinear-dependence-spec.md) was approved before implementation. Run instructions and score definitions live in the [workflow](../results/myh7-nonlinear-dependence/README.md).

## Results

The [combined CSV](../results/myh7-nonlinear-dependence/results/screen-v1/combined-correlations-v1.csv) contains 90 comparisons, exact complete-case cohorts, the original statistics, new distance scores, scaling and deletion summaries, and three simple explanatory sentences per row. The [sensitivity CSV](../results/myh7-nonlinear-dependence/results/screen-v1/dependence-sensitivities-v1.csv) exposes the retained and omitted variants for each reduced-cohort analysis. These scores are descriptive; no actual-cohort p-values are calculated.

The original outcome results and new estimates are:

- **Relative catalytic turnover versus primary unloaded velocity, C002:** 14 variants; Spearman rho +0.670330, conventional distance correlation 0.760956, U-centered squared score 0.400149. Study-deletion U scores range from 0.132741 to 0.815689. Removing Morck shared-slide variants leaves nine variants with conventional score 0.771198 and U score 0.384834. The strict short-head MVEL20 sensitivity contains only four variants; its U score of 1 is highly sample-sensitive and does not establish a perfect biological relationship.
- **Relative catalytic turnover versus relative available heads, C001:** 14 variants; rho −0.548459, conventional score 0.626376, U score 0.190071. Removing seven variants with exactly reused mutant short-head ATPase components leaves seven variants, conventional score 0.605623 and U score 0.074521. Further risky-Na exclusions leave four variants with U score −0.5. These reduced cohorts cannot be read as independent replication or proof that the remaining relationship is absent.
- **Relative available heads versus velocity, C003:** 13 variants; rho −0.570250, conventional score 0.621401, U score 0.058192. Study-deletion U scores range from −0.036299 to 0.250729. Removing Morck leaves eight variants with U score 0.064588. This dependence magnitude is more fragile than its rank correlation alone suggests.

The nine primary chemistry, ADP-geometry and motor-geometry block U scores range from −0.165585 to +0.069016. All primary blocks with a positive full-cohort U score have study-deletion ranges spanning zero. The three all-ten-feature sensitivities have negative U scores: −0.057464 for catalytic turnover, −0.216063 for available heads and −0.186214 for velocity, with only 11, 11 and 9 complete variants respectively. This run therefore produces no stable new feature-block lead. This is not an inference that biological relationships are absent: the annotations, dimensionality, missingness and independent sample size may be insufficient to reveal them.

Do not compare U squared values with rho as if they shared a scale, interpret a negative U estimate as decreasing biology, or call a conventional score larger than rho evidence of nonlinearity. Even a block score cannot distinguish an interaction from a marginal feature effect. Structural feature–feature patterns may describe correlated aspects of site position rather than separate biological pathways.

## Interpretation and next action

The development audit independently recomputed all 90 comparisons and 1,761 sensitivities using separate loop formulas, with high-precision checks for degenerate binary cases. All five original artifacts matched its reruns, and 15 focused and historical development tests passed. The [public publication audit](nonlinear-dependence-publication-audit-2026-10-02.md) separately records packaging, public provenance regeneration and byte-for-byte comparison of the three numerical CSVs. Public plans and manifests intentionally differ because their inputs and implementation hashes refer to the public layout.

The measured ATPase/velocity association remains the most actionable lead for the [existing matched-preparation validation question](matched-atpase-velocity-protocol.md). Validate matched short-head ATPase and unloaded velocity with independent preparations and their own uncertainty and WT-control accounting before making broader mechanism claims.

The current labels combine studies, constructs and assays. WT reuse, exact short-head ATPase reuse and cross-study sources prevent automatic exchangeability assumptions. Permutation calibration, multiple-testing correction, cross-validation or a more flexible learner cannot automatically repair those measurement problems. A prospective matched panel can support a small preregistered dependence family and, if prediction becomes the objective, nested study-aware comparison with a training-fold mean baseline.

The new code runs cheaply on the current scientific environment and can be reused as evidence improves. It does not predict unseen variants, establish blocked-head adoption, measure PPS stability, or revive combined-force modeling.

## Public packaging

Promoted to the public layout on 2026-10-02. Code is in `src/myh7_nonlinear_dependence/`, tests in `tests/myh7_nonlinear_dependence/`, configuration and frozen CSVs in `results/myh7-nonlinear-dependence/`, and the regenerated public plan in `data/public/myh7-nonlinear-dependence/`. The [runbook](../results/myh7-nonlinear-dependence/README.md) governs public execution paths; historical development paths above describe the original specification only.
