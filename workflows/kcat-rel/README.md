# Frozen MYH7 `kcat_rel` diagnostic

This public package preserves the strict 17-variant primary cohort and the
completed diagnostic. It includes
the frozen evidence, approved no-MSA feature profile, leakage-safe validation,
secondary uncertainty analyses, and the one authorized diagnostic record. It
is a research result, not an approved forward-prediction model. The numerical
tables are unchanged from the private snapshot; publication edits removed
machine-specific paths from manifests and documentation.

## Double-check an exact label

Open [`data/canonical-labels.csv`](data/canonical-labels.csv). Each row places
the selected WT and mutant rates beside `kcat_rel` and `ln_kcat_rel`.
`canonical_measurement_ids` links the row to
[`data/measurement-evidence.csv`](data/measurement-evidence.csv), where assay
context, replicate identity, uncertainty, source location, DOI, and matched-WT
identity are recorded.

For the 12 summary-measurement labels, the stored fields follow:

```math
k_{\mathrm{cat,rel}}
= \frac{k_{\mathrm{cat,mut}}}{k_{\mathrm{cat,WT}}},
\qquad
y = \ln\!\left(k_{\mathrm{cat,rel}}\right).
```

Here `kcat_rel` stores $k_{\mathrm{cat,rel}}$ and `ln_kcat_rel` stores the
modeling target $y$. With independently reported WT and mutant errors, the
first-order propagation used by the executable audit is:

```math
\sigma_y
= \sqrt{
    \left(\frac{\sigma_{\mathrm{mut}}}{k_{\mathrm{cat,mut}}}\right)^2
    +
    \left(\frac{\sigma_{\mathrm{WT}}}{k_{\mathrm{cat,WT}}}\right)^2
  },
\qquad
\sigma_{k_{\mathrm{cat,rel}}}
= k_{\mathrm{cat,rel}}\,\sigma_y.
```

For the five Morck 2022 variants, the canonical ratio is the paper's published
average of same-day paired biological-replicate ratios. The printed replicate
ratios and rates are retained separately; the published average was calculated
from unrounded fits and can differ by 0.005 from the mean of the displayed
two-decimal ratios. `ln_kcat_rel` is the exact natural logarithm of that
published aggregate. Its delta-method uncertainty is
$\sigma_y = \sigma_{k_{\mathrm{cat,rel}}} / k_{\mathrm{cat,rel}}$.

## Frozen selection contract

- The primary cohort is exactly 17 variants. D382Y and P710R are not primary
  rows; they belong only to the later 19-row sensitivity run.
- Every selected row is pure-actin, steady-state, actin-activated catalytic
  turnover on an accepted short-motor construct.
- sS1/S1 is selected when available. Eligible 2-hep is used only when no
  eligible sS1/S1 measurement exists.
- WT and mutant values share source study, construct, assay, and experimental
  context. The Vera external-WT pairings are therefore absent.
- NADH-coupled and colorimetric phosphate detection remain explicit metadata.
- Each variant contributes one canonical label. Morck biological replicates
  remain separate evidence rows and are aggregated before modeling, preventing
  pseudoreplication.
- Reported uncertainty and error type are preserved. Primary diagnostic fits
  remain unweighted because error types are heterogeneous.

The raw primary-value audit and publisher source files remain outside this
repository. Public evidence rows cite source locations and DOIs.

## Audit

From this directory, using the already provisioned Talk-v0 environment:

```bash
conda run -n lsar-na-rel-talk-v0 python -m pytest
PYTHONPATH=src conda run -n lsar-na-rel-talk-v0 python -m kcat_rel.evidence \
  data/canonical-labels.csv data/measurement-evidence.csv
```

The executable audit recomputes all ratios, log labels, and propagated
uncertainties; enforces the cohort and construct policy; checks eligible assay
metadata and matched WT context; verifies source locations and DOIs; and checks
that all 22 measurement rows collapse to exactly 17 canonical labels from nine
source studies.

## Frozen feature commands and historical diagnostic record

Only the checksum-pinned UniProt P12883 FASTA and RCSB 8ACT archive may enter
feature generation.  The `features` command does not accept label or evidence
inputs; audit its generated pair before the target-bearing `diagnostic` command.
From this workflow directory, download the 8ACT archive into the ignored
`cache/` directory. The feature command verifies its configured SHA-256 before
using it.

```bash
mkdir -p cache
curl -L https://files.rcsb.org/download/8ACT.cif.gz -o cache/8ACT.cif.gz

PYTHONPATH=src conda run -n lsar-na-rel-talk-v0 python -m kcat_rel.cli features \
  --reference-fasta ../../tests/lsar_na_rel/fixtures/P12883.fasta \
  --structure cache/8ACT.cif.gz \
  --config config/kcat-rel-v0.json --output-dir /tmp/kcat-feature-check

PYTHONPATH=src conda run -n lsar-na-rel-talk-v0 python -m kcat_rel.cli audit-features \
  --features /tmp/kcat-feature-check/kcat-rel-v0-no-msa-features.csv \
  --feature-manifest /tmp/kcat-feature-check/kcat-rel-v0-no-msa-feature-manifest.json \
  --config config/kcat-rel-v0.json

```

### Historical diagnostic invocation — do not replay

The committed [`results/diagnostic-v0`](results/diagnostic-v0/) is the
historical, non-replayable record of the one authorized frozen diagnostic. The
following invocation documents that historical execution only; do not rerun it
or write to its destination:

```bash
# Historical record only — do not run.
PYTHONPATH=src conda run -n lsar-na-rel-talk-v0 python -m kcat_rel.cli diagnostic \
  --features data/derived/kcat-rel-v0-no-msa-features.csv \
  --feature-manifest data/derived/kcat-rel-v0-no-msa-feature-manifest.json \
  --canonical-labels data/canonical-labels.csv \
  --measurement-evidence data/measurement-evidence.csv \
  --config config/kcat-rel-v0.json --output results/diagnostic-v0
```

Any new execution must use a fresh, empty destination distinct from
`results/diagnostic-v0` (for example, `--output /tmp/kcat-new-run`). The CLI
rejects an occupied destination before opening evidence or starting fitting,
and rechecks it immediately before its atomic write; it never deletes or
alters the existing directory. `SUCCESS` is written last and means that the
diagnostic completed, not that its forward-prediction gate passed.

The required report set includes `inner-fold-audit.json`. For future runs this
is serialized from the `inner_folds` retained by each outer validation result.
For the already executed frozen run, it is explicitly marked as a
`post_execution_reconstruction` from the canonical model table, exact outer
folds, source-group inner splits, and inner-training-only scaler statistics.
This provenance repair did not rerun, refit, or retune the diagnostic.

The recorded frozen run is a negative diagnostic: neither Ridge nor GPR passed
the predeclared gate, so no forward-prediction implementation was created. See
[`results/diagnostic-v0/AUDIT.md`](results/diagnostic-v0/AUDIT.md) and its
hash-linked `run-manifest.json` for exact metrics and provenance.

### Read-only descriptive summary of frozen artifacts

The following command renders a descriptive, post-hoc presentation summary of
the frozen artifacts. It is **not** a diagnostic rerun: it does not fit,
retune, or alter the historical record. Its output directory must be outside
Git (and must be empty or absent); never use `results/diagnostic-v0` as its
destination.

```bash
MPLCONFIGDIR=/tmp/kcat-rel-mpl-cache \
PYTHONPATH=src conda run -n lsar-na-rel-talk-v0 \
python -m kcat_rel.modeling_summary \
  --features data/derived/kcat-rel-v0-no-msa-features.csv \
  --feature-manifest data/derived/kcat-rel-v0-no-msa-feature-manifest.json \
  --model-table data/derived/kcat-rel-v0-no-msa-model-table.csv \
  --study-oof results/diagnostic-v0/study-oof.csv \
  --variant-oof results/diagnostic-v0/variant-oof.csv \
  --metrics results/diagnostic-v0/metrics.json \
  --sensitivity results/diagnostic-v0/sensitivity.json \
  --tuning-records results/diagnostic-v0/tuning-records.csv \
  --run-manifest results/diagnostic-v0/run-manifest.json \
  --config config/kcat-rel-v0.json \
  --output-dir /absolute/path/outside-the-repository
```

## Known evidence limitations

- G768R is preprint evidence. Its displayed errors are retained, but the source
  does not label them as SEM versus SD for this result.
- For summary WT/mutant pairs, first-order independent-error propagation is
  used because their covariance is not reported.
- For Morck 2022, the paper's average paired-ratio uncertainty is retained;
  replicate fit errors are not treated as independent training observations.
