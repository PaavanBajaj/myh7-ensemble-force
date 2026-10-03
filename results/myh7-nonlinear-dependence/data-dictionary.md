# Dependence-screen data dictionary

All CSVs use UTF-8 with a header row. An empty score means undefined or not applicable, never zero. IDs separated by `|` are sorted. Fields ending in `_json` contain JSON objects inside correctly quoted CSV cells. Scores are formatted to 12 significant digits.

## Combined CSV

`comparison_id` preserves C001–C078 and adds B001–B012. `pair_type` distinguishes outcome/outcome, outcome/feature, feature/feature and block/outcome comparisons. `analysis_role` identifies a historical pair, a primary block or the all-feature sensitivity. `x`/`y` are exact variables or block IDs; `x_columns`/`y_columns` list the underlying feature names. `x_dimensions`/`y_dimensions` count requested columns; `x_varying_dimensions`/`y_varying_dimensions` count columns remaining after constants are removed.

`original_method`, `original_status`, `original_effect` and `original_effect_unit` are copied verbatim from the frozen screen. `spearman_rho` copies rho only for the 66 original Spearman rows. It is blank for 12 binary rows whose original effect is median(flag=1) minus median(flag=0), and for all 12 blocks. `spearman_blank_reason` explains those blanks. `group_0_n`/`group_1_n` preserve binary group sizes. `original_loo_variant_min/max` and `original_loo_study_min/max` preserve historical effect ranges, whose units follow `original_effect_unit`.

`n_paired`, `paired_variants`, `n_universe`, `n_missing_either` and `missing_variants` expose complete-case coverage against the 22-variant union. `missing_columns_by_variant_json` records exactly which required values are missing; `feature_missing_reasons_json` copies available structural missing reasons, otherwise stating that the value is missing in the frozen table. This fallback does not invent an unmeasured biological reason. `n_studies` and `study_ids` use the union of outcome sources. `x_study_composition`/`y_study_composition` count source membership on each outcome side and are blank for intrinsic features.

`distance_correlation` is conventional nonnegative sample distance correlation in [0,1], not its square. `u_centered_squared_distance_correlation` is the signed normalized U-centered inner product, a bias-corrected squared distance-correlation score in [−1,1]; its ratio is not unbiased. Neither gives an increasing/decreasing association direction. Negative U values are retained and usually indicate finite-sample fluctuation around a weak/null signal, not inverse correlation.

`status` applies to the conventional score: `estimated`, `n_lt2` or `constant_side`. `u_status` applies to the U score: `estimated`, `n_lt4`, `constant_side` or `zero_u_distance_variance`. Binary singleton groups can have mathematically zero U variance even while conventional distance correlation is defined. U zero variance is detected with a self-norm tolerance of 64 machine epsilons relative to the raw-distance norm; the cross-product is never thresholded.

`x_constant_columns`/`y_constant_columns` list removed constant dimensions. `x_scale_sample_sd_json`/`y_scale_sample_sd_json` record feature standard deviations used for the raw primary comparison. Outcomes are not standardized; scalar distance correlation is scale invariant. Scaling is descriptive within each retained cohort and does not use outcome labels.

`rank_distance_correlation`, `rank_u_centered_squared_distance_correlation` and `rank_status` describe average-rank sensitivity. `log_outcome_distance_correlation`, `log_outcome_u_centered_squared_distance_correlation` and `log_outcome_status` describe natural-log-outcome sensitivity; they are blank/not applicable for feature-only pairs. `rank_status`/`log_outcome_status` use the U-status vocabulary.

`loo_variant_distance_min/max`, `loo_variant_u_min/max`, `loo_study_distance_min/max` and `loo_study_u_min/max` summarize estimable deletion scores. The `loo_variant_estimable`/`loo_study_estimable` counts apply to U scores; conventional-score counts can differ if a deletion makes U variance zero. Feature-only study fields are blank with `loo_study_estimable=not_applicable_feature_pair`. Exact deletions and undefined scores remain visible in the sensitivity CSV. These ranges are not statistical confidence intervals or a cross-validation accuracy measure.

`interpretation_limit` states the shared scientific boundaries. `description` supplies three plain-language sentences about the comparison and its main limitations. A block row describes its features jointly and does not attribute the score to any one feature or prove an interaction.

## Sensitivity CSV

`scenario` identifies a variant/study deletion or a shared-control/assay exclusion. `omitted_unit` is the deleted variant or source study where applicable. `n_paired`, `paired_variants` and `excluded_variants` expose exactly what remained and what was removed relative to that comparison's complete-case primary cohort. `study_ids` is recalculated on retained rows. `original_spearman_rho` is included only when that precise selected sensitivity cohort had rho in the historical screen. `note` explains selection and limitations. Score, status, constant-dimension and scale fields have the same definitions as the combined CSV and are recomputed on each subset.

## Point CSV and manifests

The point CSV has one row per comparison, variant, side and underlying column, with `raw_value` and `source_study`. It can be joined to the combined CSV using `comparison_id` and makes multivariate inputs inspectable without treating repeated feature values across comparisons as extra independent observations.

The plan freezes method settings, all comparison definitions and variant IDs, nine source-file hashes, configuration hash and implementation hash. The run manifest records its own plan hash, the corresponding inputs/code/config hashes, output-file hashes and Python/NumPy/SciPy versions. Neither contains absolute local paths or changing timestamps. `association_estimates_computed` distinguishes planning from analysis; `inferential_p_values_computed=false` explicitly records the descriptive analysis boundary.
