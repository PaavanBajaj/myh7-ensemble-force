# Frozen `kcat_rel` diagnostic audit

| Gate | Result |
| --- | --- |
| Evidence audit | PASS (17 labels, 22 measurements, 9 studies) |
| ridge primary gate | FAIL: study_macro_mae, variant_mae, study_macro_rmse, variant_rmse, single_study_influence |
| gpr primary gate | FAIL: study_macro_mae, variant_mae, study_macro_rmse, variant_rmse, single_study_influence |

Selected model: none; negative diagnostic.

## Post-execution provenance repair (no diagnostic rerun)

This artifact-only repair records that the execution occurred at Git HEAD
`91affadb0f5bdb3fbf8f09b22209836f38ddf9c5`, before the Task 6 commit existed.
The CLI/reporting code used for that run was present as uncommitted working-tree
source. `run-manifest.json` now carries reconstructed content hashes for that
source snapshot and the relevant frozen model modules/config; it does not claim
that commit `7819ffa` existed during execution. The derived model table is
byte-identical to this run's `model-table.csv`; this repair did not fit, tune,
or otherwise alter any diagnostic result.

`inner-fold-audit.json` was added by deterministic post-execution
reconstruction from this run's canonical model table and outer tuning records.
It enumerates every model, outer view, outer fold, and source-group inner split,
including the scaler statistics recomputed only from that inner training set.
Its provenance field says `post_execution_reconstruction`; it is not presented
as an artifact captured during the original execution, and no diagnostic was
rerun.
