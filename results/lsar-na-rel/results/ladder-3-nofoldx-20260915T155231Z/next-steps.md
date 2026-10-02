# Next steps (Talk-v0 diagnostic ladder)

This run is an exploratory small-n private diagnostic. Risky allowed rows
R403Q, R663H, and G768R are included. Do not treat pooled OOF errors as a
clinical or mechanistic claim.

No FoldX-derived values are present in this release. The modeling family is
chemistry plus checksummed 8ACT features from the committed features.csv
snapshot. The IHM interface flag remains 2/16 positives
(D382Y, R403Q) at the frozen 4.5 Å cutoff; do not retune geometry from these
scores.

Re-derive metrics from metrics.json and oof-predictions.csv in this directory.
Ship the other pre-registered ladder run_ids on the same rows and folds.
Do not add a fourth feature combination or drop rows after seeing scores.
