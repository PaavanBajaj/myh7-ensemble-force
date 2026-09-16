# LSAR Na-rel workflow

This is the public S3 workflow for the fixed 16-row continuous-LSAR cohort. It
generates outcome-free chemistry and 8ACT features, then evaluates three
registered Gaussian-process models with leave-one-`source_id`-out validation.
It is exploratory research software, not a clinical or mechanistic result.

## Supported environment

The exact lock currently supports `osx-arm64` only.

```bash
conda create --name lsar-na-rel-exact --file environment-osx-arm64.lock
conda activate lsar-na-rel-exact
export LIBCIFPP_DATA_DIR="$CONDA_PREFIX/share/libcifpp"
export PYTHONPATH=src
```

`environment.yml` contains the same direct dependency versions. The explicit
lock additionally freezes every transitive package URL, build, and SHA-256.

## Repository-local inputs

All default paths are derived from the checkout root:

- `data/public/lsar-na-rel/na-rel-primary.csv`
- `data/public/lsar-na-rel/labels.csv`
- `data/public/lsar-na-rel/features.csv`
- `data/public/lsar-na-rel/features-provenance.json`
- `tests/lsar_na_rel/fixtures/P12883.fasta`
- `workflows/lsar-na-rel/config/lsar-na-rel.json`
- `workflows/lsar-na-rel/cache/`
- `workflows/lsar-na-rel/results/`

No sibling repository or absolute user path is required.

## Generate or verify features

Generation downloads the full 8ACT biological assembly into the ignored cache,
verifies its pinned compressed and decompressed SHA-256 values, runs DSSP, and
atomically writes the feature snapshot and provenance.

```bash
PYTHONPATH=src python -m lsar_na_rel.features
```

Verify the committed snapshot without replacing it:

```bash
PYTHONPATH=src python -m lsar_na_rel.features --check
```

The generator reads only the `variant` column from the primary table. It does
not use `Na_rel`, `ln(Na_rel)`, labels, source groups, or model scores.

## Run the registered model ladder

```bash
PYTHONPATH=src python -m lsar_na_rel.cli talk-v0 --all-ladder
```

A single registered run can be selected with `--run-id ladder-1-rsa`,
`ladder-2-rsa-ihm`, or `ladder-3-nofoldx`. The config rejects added,
removed, reordered, or renamed registered runs.

Each successful timestamped result directory contains
`run-manifest.json`, `metrics.json`, `oof-predictions.csv`,
`observed-vs-predicted.png`, `run.log`, `next-steps.md`, and a
last-written `SUCCESS` marker.

## Frozen modeling contract

| Run ID | Feature columns |
| --- | --- |
| `ladder-1-rsa` | `wt_site_rsa_8act_mean_ab` |
| `ladder-2-rsa-ihm` | `wt_site_rsa_8act_mean_ab`, `ihm_interface_flag_8act` |
| `ladder-3-nofoldx` | `delta_charge`, `wt_site_rsa_8act_mean_ab`, `ihm_interface_flag_8act` |

All runs use training-fold-only standardization, an isotropic
`ConstantKernel * Matern(nu=1.5) + WhiteKernel` GPR, ten optimizer restarts,
`normalize_y=True`, and seed 20260811. The baseline is the training-fold mean.
The configured `talk-v0-foldx` profile is unavailable and is not a registered
run because no FoldX-derived values are present in this release. Exploratory
distance probes are intentionally excluded from the public runnable workflow.

## Test

```bash
PYTHONPATH=src pytest -q
```

The release suite contains 109 passing tests. Optimizer convergence warnings
from scikit-learn are expected under the frozen small-sample setup.
