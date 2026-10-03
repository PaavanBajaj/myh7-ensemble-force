# Nonlinear and multivariate dependence tests for the MYH7 screen

Date: 2026-10-02. Research question: which runnable methods can reveal relationships beyond the monotone patterns measured by Spearman rho, given the current MYH7 feature and outcome cohorts?

## Recommendation and meaning of “best”

**Implement U-centered distance covariance and its normalized, bias-corrected squared distance correlation first.** Use it for both the existing pairwise inventory and a small, frozen set of multivariate feature blocks against each outcome. Preserve Spearman rho alongside it to describe monotone direction. This is a concrete recommendation for this dataset, not a claim that one test has uniformly greatest power.

Distance covariance has an independence-characterizing population definition for random vectors of arbitrary dimensions with the required moment conditions. Consequently it can detect nonmonotone patterns and multivariate dependence that a univariate rank correlation can miss. These are population properties and asymptotic consistency statements, not a promise of reliable detection from 15 variants. [Székely, Rizzo and Bakirov, 2007](https://arxiv.org/abs/0803.4101).

My recommendation favors few choices, transparent scores, inexpensive computation, and easy comparison with the frozen screen. It requires no learned predictor, label-selected feature subset, Gaussian-kernel bandwidth, or neighbor-scale search. Distance correlation is therefore a suitable **general dependence test**, even though the requested output calls it a “new ML score.” Label that column with its actual statistic name, not accuracy or predictive performance.

The defensible definite answer is about **which analysis to run next**. No method can provide a definite biological relationship, validated mechanism, independence conclusion from a nonsignificant result, or generalizable predictor from these observations alone.

## Actual data constraints

Grounding: the repository [screen README](../../results/myh7-correlation-screen/README.md), [runbook](../../results/myh7-correlation-screen/README.md), and analysis table identify 17 canonical `kcat_rel`, 16 canonical `Na_rel`, and 15 primary `v_rel` variants; 10 features; and 78 frozen pairwise comparisons. Five geometry annotations describe WT structure, not measured mutant perturbations. Binary IHM membership has sparse groups. This is a modest number of dimensions with extremely few biological units, rather than a dataset with thousands of independent measurements.

The existing Spearman estimates are descriptive, without p-values. The 32 Morck slide/channel pairs represent five variants and must not be promoted to 32 independent variant rows. Study/construct/assay differences, shared WT controls, and missing features affect the effective information available.

Seven of 14 overlapping `kcat_rel`–`Na_rel` variants reuse exact mutant short-head ATPase components. Where matching WT components are also shared, `Na_rel = L_rel / S_rel`, so an association with short-head ATPase can arise by construction. Every dependence statistic inherits this problem. Excluding exact-S reuse and documenting uncertain denominators is necessary; changing the test does not repair measurement coupling.

These facts come from the checked-in ledgers, not from the statistical papers. Any implementation should freeze their hashes and retain one row per variant.

## Method comparison

### Distance covariance/correlation: recommended first

Use the 2014 U-centered estimator instead of relying only on the usual nonnegative, biased sample coefficient. Its squared distance-covariance numerator is unbiased under the sampling assumptions. Its normalized ratio is appropriately called **bias-corrected squared distance correlation**, not an unbiased estimator of correlation. The sample ratio can be negative; a negative estimate is finite-sample fluctuation, not an inverse biological relationship. Preserve it in output. [Székely and Rizzo, 2014](https://arxiv.org/abs/1310.2926), [dcor estimator documentation](https://dcor.readthedocs.io/en/latest/functions/dcor.u_distance_correlation_sqr.html).

The ordinary nonnegative distance correlation is useful as a familiar descriptive companion. It has no biological sign: strong increasing, decreasing, or curved relationships may all produce a high score. Comparisons between Spearman magnitude and distance correlation are qualitative because their scales and targets differ. Use the normalized U statistic for the chosen permutation test; its denominator remains fixed under row permutation.

### HSIC/kernel independence: strong second option

HSIC compares joint dependence through kernels and has quadratic sample-size computational cost. Suitable characteristic kernels permit general dependence detection; Gaussian kernels are a usual continuous-data choice. Its original paper develops independence testing and the authors' software offers shuffling as well as an approximation, with shuffling particularly relevant for small samples. [Gretton et al., 2007/2008](https://www.gatsby.ucl.ac.uk/~gretton/papers/GreFukTeoSonetal08.pdf), [author implementation](https://www.gatsby.ucl.ac.uk/~gretton/indepTestFiles/indep.htm).

For this dataset, freeze a Gaussian-kernel rule using each marginal's distances, keep it invariant through permutations, and calibrate by valid permutations instead of a large-sample approximation. Choosing bandwidths to maximize the observed feature-label association would need the complete bandwidth selection repeated under every permutation. I recommend postponing this additional tuning choice.

Distance covariance and kernel methods are closely related: distance-based statistics can be expressed through corresponding kernels, and the distance/kernel choice affects power. Agreement between these methods is useful robustness information, not independent replication. [Sejdinovic et al., 2013](https://arxiv.org/abs/1207.6076).

### MGC: useful nonlinear sensitivity, not a guaranteed winner here

MGC uses local distance correlations over neighbor scales and a smoothed scale-selection rule. Its theory establishes universal consistency and simulations show benefits over distance correlation for several nonlinear relationships. The applied paper benchmarks dimensions from 1 to 1000 and reports sample-efficiency benefits in its benchmark suite. Those benchmarks do not establish superiority for heterogeneous MYH7 studies with 15 variants. [Shen, Priebe and Vogelstein, JASA 2020](https://arxiv.org/abs/1710.09768), [Vogelstein et al., eLife 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6386524/).

SciPy implements `(n,p)` versus `(n,q)` MGC, returns a score, optimal scale, scale map, and permutation p-value, and documents a minimum of five samples. A software minimum is not an adequate sample-size justification. Its public interface exposes unrestricted permutations rather than study blocks, so the default p-value is inappropriate as this project's primary calibration. If implemented later, recompute the complete MGC scale-selection statistic under each allowed block permutation. [SciPy MGC documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.multiscale_graphcorr.html).

### Mutual information and MIC/TIC: defer

Population mutual information captures general dependence. Kraskov-type estimation uses nearest-neighbor distances with a bias/resolution tradeoff involving `k/N`. With roughly 15 rows, mixed discrete features, tied values, and a 10-dimensional joint feature vector, there is little data to estimate local densities. That practical limitation is my inference from the estimator design, not a formal impossibility theorem for n=15. [Kraskov, Stögbauer and Grassberger, 2004](https://doi.org/10.1103/PhysRevE.69.066138).

The original MIC equitability claims have been contested; authors distinguish relationship-strength ranking from testing independence. Newer MICe targets equitability, while TICe targets powerful independence testing and has consistency results. This distinction matters: a high exploratory dependence score is not automatically a valid test. Their good simulation performance does not answer this project's study-exchangeability question. [Kinney and Atwal, 2014](https://doi.org/10.1073/pnas.1309933111), [Reshef et al., 2016](https://www.jmlr.org/beta/papers/v17/15-308.html).

The common scikit-learn `mutual_info_regression` routine estimates MI separately between each feature and the outcome; it is not a joint feature-vector independence test. Its discrete-feature handling, neighbor count, and tie-breaking noise must be configured correctly. It should not be presented as an interaction detector across all features. [scikit-learn API](https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.mutual_info_regression.html).

### Rank alternatives

Kendall tau provides another monotone association description but does not solve the nonmonotone or joint-interaction question. Rank-transforming coordinates before distance correlation is a reasonable frozen sensitivity for marginal outliers and units, but yields a different statistic and may lose meaningful magnitudes. Ties need deterministic average ranks; a binary coordinate remains binary in information content.

HHG is a multivariate test based on ranks of distances, applicable in any dimension and consistent against dependent alternatives under its conditions. It is a credible later sensitivity. Do not confuse ranks of multidimensional distances with Spearman applied to a single coordinate. [Heller, Heller and Gorfine, 2013](https://arxiv.org/abs/1201.3522).

## A runnable, fixed first analysis

1. Retain all 78 frozen pairwise IDs and exactly their complete-case cohorts. Compute ordinary descriptive distance correlation and U-centered squared distance correlation beside the original rho or binary median difference. Do not silently replace a binary group difference with “Spearman rho.”
2. Freeze three omnibus comparisons: all 10 features versus `kcat_rel`, versus `Na_rel`, and versus primary `v_rel`. Complete-case rows must be explicit. If there are too few rows or no informative feature columns, emit a documented nonestimable result rather than filling structural values with zeros.
3. Optional explanatory blocks may split sequence descriptors and structural geometry. Specify exact memberships before estimates; the two sequence descriptors and seven continuous geometry/RSA columns are natural distinct blocks, with the binary IHM annotation separately identified. Every extra block increases the reporting/test family. Correlated geometry columns should not be interpreted as independent pathways.
4. Use Euclidean distances with exponent 1. For multivariate blocks, standardize each nonconstant feature to zero mean and unit standard deviation using only its observed marginal values in the frozen complete-case cohort. Record the standard-deviation convention. Keep that transform and feature membership fixed under permutations. This equalizes units but is a modeling decision; correlated geometry can still dominate distance contributions. Avoid covariance whitening at this sample size.
5. Primary outcome values can stay on their original ratios for continuity with the Spearman screen. Freeze natural-log outcomes as a secondary sensitivity. Scalar centering/rescaling does not change normalized distance correlation; nonlinear transformations do.
6. Keep variant and study deletion analyses as descriptive stability ranges. Refit the declared marginal scaling rule on each explicitly reduced cohort if that is the frozen deletion procedure, and record the changing cohort. These ranges are not confidence intervals.
7. Freeze no-Morck and no-exact-S-reuse sensitivities where relevant. Do not treat an overall nonlinear score with weak Spearman as evidence that a nonlinear mechanism exists; check plotted points, influential variants, source composition, and coupled measurements first.

These are project-specific design decisions derived from the constraints above, not prescriptions from a benchmark paper.

## Exact statistic and implementation checks

For `n > 3`, define Euclidean distances `a_ij = ||x_i-x_j||` and `b_ij = ||y_i-y_j||`, with zero diagonals. For off-diagonal entries:

`A_ij = a_ij - row_sum_i/(n-2) - column_sum_j/(n-2) + total_sum/((n-1)(n-2))`.

Set `A_ii = 0`; likewise define `B`. Then

`Omega_xy = sum_ij A_ij B_ij / [n(n-3)]`

and

`u_dcor_squared = Omega_xy / sqrt(Omega_xx Omega_yy)`.

If either self-product is zero, report a degenerate-variable status; do not manufacture an inferential result. These formulas and the possibility of negative U estimates are documented by the official implementation. [dcor theory](https://dcor.readthedocs.io/en/latest/theory.html).

Verify symmetry, affine scalar scale invariance, identical nondegenerate input giving one, exact U-centered zero diagonals/row sums, a known negative example, and agreement with independently computed reference fixtures. Test a deliberately symmetric U-shaped relationship where Spearman is near zero but distance dependence is positive. That example validates capability, not power for this data.

An ordinary NumPy implementation is tiny and makes strict within-study enumeration transparent; cross-check against `dcor.u_distance_correlation_sqr` or `hyppo.independence.Dcorr(bias=False)`. Hyppo offers block permutations and an `auto` approximation flag; explicitly disable approximation and inspect block semantics. Its hierarchical blocks can also exchange whole blocks unless marked fixed, so merely passing positive study IDs may not implement the intended restriction. [Hyppo Dcorr API](https://hyppo.neurodata.io/api/generated/hyppo.independence.dcorr).

## Permutations: what is and is not justified

**Pooled unrestricted permutations are not primary evidence here.** Under the null, the allowed reassignments must preserve the relevant sampling distribution. Independence across observations and exchangeability are different requirements; dependence does not automatically prohibit permutation if its entire joint distribution is invariant, but publication labels alone do not prove that invariance. Restricted exchangeability blocks help preserve dependence structure only where their assumptions hold. [Winkler et al., 2014](https://doi.org/10.1016/j.neuroimage.2014.01.060), [official PALM explanation](https://fsl.fmrib.ox.ac.uk/fsl/docs/statistics/palm/exchangeability_blocks.html).

For feature-outcome comparisons, a **within-outcome-study shuffle** is an assumption-dependent sensitivity: hold features fixed and permute the outcome assignment only among variants from the same study. A common WT denominator can be preserved as study context, but unequal shared-control reuse, variant-specific uncertainty, assay differences within a study, or selected variants may violate within-study exchangeability. Label the resulting number `within_study_permutation_p_assumption_dependent`, not “validated significance.”

For outcome-outcome pairs, intersect the two study assignments so a shuffle preserves both source contexts. This is conservative and may leave little or no movable data. Measurement components reused between the outcomes remain a separate violation of the scientific independence interpretation; run the exclusion sensitivity rather than explaining them away as blocks.

From the checked-in analysis table, the full outcome cohorts have these within-study index-permutation spaces:

- `kcat_rel`: 17 variants across 9 studies; `5! × (2!)^4 = 1,920`.
- `Na_rel`: 16 variants across 6 studies; `(5!)^2 × (2!)^2 = 57,600`.
- `v_rel`: 15 variants across 9 studies; `5! × (2!)^2 = 480`.

These are upper bounds before complete-case filtering or further exchangeability restrictions. Tied values reduce distinct outcomes. Singleton studies cannot move. The velocity calibration is consequently driven by only three movable studies, including the five-variant Morck cohort. A blocked test cannot establish absence of across-study dependence when its null reassignments only assess within-study pairing.

Enumerate all allowed index permutations when feasible, including the identity, and use `count(T_perm >= T_observed)/M` with a declared numerical tie tolerance. Counting equal-sized index permutations remains valid with ties; do not deduplicate statistic values and give each distinct statistic equal weight. For large spaces draw a fixed-seed uniform sample of permutations and use `(1 + count(T_perm >= T_observed))/(B + 1)`, with recorded `B` and Monte Carlo precision. Never emit zero p-values. [Phipson and Smyth, 2010](https://pubmed.ncbi.nlm.nih.gov/21044043/).

If there is only the identity, emit `not_testable_no_exchangeable_reassignments` and leave the p-value blank. A nominal p=1 obscures lack of testability. If the declared finer blocks cannot be established from provenance, keep the statistic descriptive and expose that limitation. Cluster bootstrap or residualizing study indicators does not automatically restore valid exchangeability; a new sampling/conditional-inference argument would be required.

## Multiplicity and power

The three frozen all-feature omnibus comparisons are the clearest small inferential family. Any formal claims must control that family, and the 78 exploratory pairwise tests must not be searched as unadjusted discoveries. Holm or conservative Bonferroni works under arbitrary dependence when individual p-values are valid. [Holm, 1979](https://www.jstor.org/stable/4615733). BY is a conservative FDR alternative under general test dependence; ordinary BH needs its appropriate dependency assumptions. No adjustment repairs an invalid permutation model. [Benjamini and Yekutieli, 2001](https://doi.org/10.1214/aos/1013699998).

Do not borrow high-dimensional asymptotic approximations as a solution to 15 heterogeneous observations. Theory shows standard joint distance/kernel statistics can reduce largely to componentwise linear behavior in certain regimes where dimension grows rapidly relative to n, motivating marginal aggregation. Other work demonstrates nonlinear capability in moderately high-dimensional regimes. These results are not contradictory guarantees at this project's finite sample size, and the present 10-feature setting does not directly reproduce either asymptotic regime. [Zhu et al., 2020](https://arxiv.org/abs/1902.03291), [Gao et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8491772/).

An omnibus vector test can catch interactions that marginal screens miss in principle, but adding irrelevant or redundant dimensions can reduce useful signal. Pairwise tests, a compact domain-defined block, and the full block should be interpreted together. No label-selected block search is allowed without incorporating that search in calibration. A null result means insufficient evidence under the chosen test and assumptions, not proof of independence.

## Acceptance criteria for downstream implementation/reporting

- Existing canonical labels, features, 78 IDs, public sibling artifacts, and original descriptive Spearman outputs stay frozen.
- Score column names distinguish ordinary distance correlation from bias-corrected **squared** distance correlation. Negative U values survive export.
- Every row includes cohort, variant count, feature membership, source/assay warnings, deletion stability, and explicit p-value status/assumptions.
- The CSV contains both binary group contrast and rho fields as applicable, and a plain two-to-three-sentence description explaining what the row compares and why its evidence is limited.
- Omnibus rows do not pretend that one multivariate score identifies a causal feature or predicts unseen variants.
- The final conclusion can prioritize a validation experiment but cannot promote a high score to a biological mechanism or successful predictor.

## Research handoff

**Completed:** primary-paper and official-implementation review; concrete first-method choice; exact estimator specification; actual full-cohort permutation counts; implementation and interpretation safeguards.

**Recommended first method:** U-centered distance covariance with normalized bias-corrected squared distance correlation, ordinary distance correlation as a companion, and frozen marginal/full-feature comparisons.

**Recommended calibration:** descriptive primary reporting plus explicitly assumption-dependent within-study permutation sensitivity, exact enumeration for feasible spaces, complete test-family correction where interpreted. More defensible prospective inference requires matched independent preparations and a prespecified test.

**Unresolved by existing data:** within-study exchangeability is not proven; unequal WT reuse covariances and reliable independent preparation counts are unavailable for portions of the cohort. No reviewed test removes those limitations.
