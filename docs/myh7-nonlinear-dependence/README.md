# Which analysis should follow Spearman?

The recommended first addition is **descriptive distance correlation**, calculated both pairwise and between a small biological feature block and one measured outcome. It can expose curved or joint relationships that a monotonic scalar rank coefficient misses. This is a choice for the current MYH7 data, not a claim that one method is universally best. See the [original distance-correlation paper](https://arxiv.org/abs/0803.4101) and [estimator definitions](https://dcor.readthedocs.io/en/latest/theory.html).

## Read the evidence

- [Dependence-method review](dependence-methods.md): distance covariance, HSIC, MGC, information estimators, rank alternatives, permutation assumptions, multiplicity and dimensionality.
- [Predictive-ML review](ml-alternatives.md): ridge/elastic net, PLS, sparse CCA, GAM, kernel models and forests against the project's previous negative diagnostics.
- [Approved implementation spec](../myh7-nonlinear-dependence-spec.md).
- [Runnable workflow and combined CSV](../../results/myh7-nonlinear-dependence/README.md).
- [Actual-run findings and scientific next step](../myh7-nonlinear-dependence-decision.md).
- [Independent execution and numerical verification](../nonlinear-dependence-publication-audit-2026-10-02.md).

## Why this recommendation fits

The research homes have the same frozen analysis table, pair inventory and result tables. The nonlinear extension is now published alongside the historical screen; its three numerical CSVs preserve the verified development results byte for byte. There are 22 variants in the union, but only 17 relative catalytic-turnover, 16 relative-available-heads and 15 primary unloaded-velocity labels. The full ten-feature complete-case sets shrink to 11, 11 and 9 variants. Source-study counts are only 9, 6 and 9. A flexible predictor therefore introduces many more fitting choices than the available independent information can support.

The previous leakage-safe ridge and Gaussian-process diagnostics did not clear their prediction gates. A dependence statistic asks a different, narrower question: do variants that differ in these measurements or annotations also differ in the outcome? No predictor is trained. A score does not establish an interaction, identify which feature causes a joint relationship, or measure predictive accuracy.

The selected workflow adds two complementary distance scores:

- Conventional sample distance correlation, an unsigned descriptive magnitude between zero and one. Its small-sample upward bias makes it unsuitable for comparison with rho as if the two were interchangeable.
- The normalized U-centered, bias-corrected **squared** distance-correlation score. Negative estimates remain visible; their sign is estimator fluctuation, not decreasing biology. The covariance numerator is unbiased under the iid sampling model, but the normalized ratio is not an unbiased correlation.

These distinctions follow the [official estimator theory](https://dcor.readthedocs.io/en/latest/theory.html). Compact blocks reduce choices and avoid combining every partly redundant structural coordinate as primary evidence. Chemistry, ADP geometry and motor geometry are fixed by prior feature meaning, not by new scores. All ten features form a separate coverage/dimensionality sensitivity. The two independent reviews differ on whether the full-feature omnibus should be primary; this implementation favors compact primary groups because the full panel loses observations and adds correlated geometry. That is a documented project decision, not an empirical method-ranking result.

## Other methods worth knowing

**HSIC** is a strong kernel-based alternative when a biologically justified kernel is available; the kernel and bandwidth add decisions. **MGC** can adapt to local nonlinear geometry and has promising benchmark results, but those benchmarks do not prove superiority for this cohort. Both address dependence rather than generalizable prediction. See the [HSIC primary paper](https://proceedings.neurips.cc/paper/2007/file/d5cfead94f5350c12c322b5b664544c1-Paper.pdf) and [MGC primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6386524/).

Nearest-neighbor mutual information requires estimating local density from very few points. The common [scikit-learn routine](https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.mutual_info_regression.html) scores each feature separately; supplying a multicolumn matrix does not make it a joint interaction test. MIC/TIC and distance-rank tests are later robustness options, with their own estimation and calibration choices. Detailed primary-source discussion is in the dependence review.

For a future prediction question with substantially better independent measurements, one-component PLS or a strongly regularized, biologically restricted model is a reasonable first challenger. A single-feature curved hypothesis may favor a constrained GAM. PLS projections, transformations, selection and hyperparameters must all be trained within nested study-aware validation. More flexible kernels and trees can run now, but a successful fit would not overcome the dataset's limitations.

## What the data cannot answer

Seven ATPase/available-heads pairs reuse a short-head ATPase component. Available heads is constructed from a long-to-short turnover ratio, so an inverse relationship with short-head turnover may partly be algebraic. Five Morck velocity variants reuse WT slide/channel records. Mixed constructs, temperatures and assays add differences between studies. Static WT site geometry is not mutant dynamics or measured PPS stability.

The workflow reports variant deletion, source-study deletion, shared-component exclusions and scale sensitivity. These ranges are descriptions, not confidence intervals. No actual-cohort permutation p-values or discovery claims are calculated: neither pooled nor within-study exchangeability has been established. Even a [cluster-aware HSIC framework](https://www.nature.com/articles/s41598-022-26278-9) has specific sampling assumptions that heterogeneous unequal studies do not automatically satisfy. A valid multiplicity adjustment cannot repair invalid individual p-values.

The definite answer is **which exploratory analysis to run next**. The definite biological answer still requires matched independent experiments. The existing [matched ATPase/velocity validation question](../matched-atpase-velocity-protocol.md) remains the leading scientific follow-up; the new dependence screen can refine priorities without replacing that validation.
