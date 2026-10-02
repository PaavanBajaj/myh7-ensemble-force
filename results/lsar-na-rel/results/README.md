# Regenerated LSAR Na-rel results

These result packs were generated from this public checkout on 2026-09-15 with
the repository-local command:

```bash
PYTHONPATH=src python -m lsar_na_rel.cli talk-v0 --all-ladder
```

Every directory has the configured seven artifacts and a last-written
`SUCCESS` marker. Manifests record the exact input hashes, environment
versions, dirty-tree state, and a hash of runnable workflow content.

| Run ID | Result directory | GPR MAE | GPR RMSE | Baseline MAE | Baseline RMSE |
| --- | --- | ---: | ---: | ---: | ---: |
| `ladder-1-rsa` | `ladder-1-rsa-20260915T155230Z` | 0.13120166676711030 | 0.16785061947218140 | 0.13120166676711023 | 0.16785061947218136 |
| `ladder-2-rsa-ihm` | `ladder-2-rsa-ihm-20260915T155230Z` | 0.13120166676711037 | 0.16785061947218144 | 0.13120166676711023 | 0.16785061947218136 |
| `ladder-3-nofoldx` | `ladder-3-nofoldx-20260915T155231Z` | 0.13088618764022130 | 0.16730405908319063 | 0.13120166676711023 | 0.16785061947218136 |

The values above were independently recomputed from each
`oof-predictions.csv` and checked against `metrics.json` at absolute
tolerance below `1e-15`.

These are pooled out-of-fold descriptive errors for a fixed 16-row,
six-study-group cohort that includes three “Risky allowed” rows. The tiny
differences among runs do not establish clinical utility, mechanism,
statistical significance, or a winning model. No FoldX-derived values are
present. Talk-v0 is not a forward-prediction model and does not produce the
combined ensemble-force estimate.
