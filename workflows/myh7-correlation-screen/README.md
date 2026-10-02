# MYH7 correlation screen, 2026-10-02

This public workflow supports the [dated decision](../../docs/myh7-correlation-screen-decision.md) and [validation protocol](../../docs/matched-atpase-velocity-protocol.md). It is an exploratory comparison of measured `kcat_rel`, canonical LSAR-derived `Na_rel`, unloaded motility `v_rel`, and 10 distinct features. It does not train a model or revise any frozen diagnostic or published LSAR release.

## Read the results

Start with the [decision](../../docs/myh7-correlation-screen-decision.md), then the [data dictionary](data-dictionary.md). The [comparison inventory](data/comparison-inventory-v1.csv) was frozen before estimates; [results](data/comparison-results-v1.csv), [point-level rows](data/comparison-points-v1.csv), and 78 [labeled SVG plots](data/plots/C001.svg) through `C078.svg` expose every eligible pair. The [sensitivity table](data/sensitivity-results-v1.csv) reports shared-measurement, study, assay, and expanded-cohort checks. The full screen includes weak and null patterns rather than selecting only high correlations.

The primary outcome counts are 17 `kcat_rel`, 16 `Na_rel`, and 15 Table 1-aligned `v_rel`; 14, 14, and 13 variants overlap across the three outcome pairs in order. Three additional human HCM velocity leads have a separate sensitivity cohort. The [velocity cohort ledger](data/velocity-cohort-ledger-v1.csv) records 15 primary, 8 unresolved Table 1, 3 expanded, and one excluded LVNC control. The [kcat–Na dependence ledger](data/dependence-ledger-v1.csv) and [velocity dependence ledger](data/velocity-dependence-v1.csv) identify reused components and controls. The 32 selected Morck measurements are slide/channel pairs for five variants, not 32 biological variants.

## Feature selection and inputs

The [selection lock](feature-selection.md) predates all association estimates. The [primary-literature review](feature-literature.md) ranks candidates and defines exactly five new WT structural annotations for velocity: ADP-to-rigor site shift, ADP adenine proximity, rigor actin proximity, relay-landmark proximity, and ELC proximity. The [feature manifest](data/feature-manifest-v1.json) pins the P12883 sequence, 8ACT/8EFD/8EFE/8EFI coordinates, chain mapping, and output hashes. [Per-site mapping](data/site-mapping-v1.csv) and [missing reasons](data/feature-coverage-v1.csv) cover the 22-variant analysis union. Existing charge, Grantham, catalytic-motif, RSA, and IHM annotations are comparators; shared `delta_charge` was deduplicated.

[Source manifest](data/source-manifest-v1.json) lists public input hashes and structural provenance. Restricted source documents and the Morck workbook are excluded; the compact derived Morck pair table and its source locations are included. The published LSAR pack in this repository supplies the canonical `Na_rel` inputs.

## Reproduce

Run from this repository root. `screen.py features` needs Python with NumPy 2.5.2 and Biopython 1.88; the optional workbook extraction step needs openpyxl 3.1.5 and access to the source workbook. The checked-in Morck pair table lets readers reproduce the published screen without that workbook. Source files and downloaded coordinates must match the SHA-256 values in the manifests. Retrieve `8ACT.pdb1.gz`, `8EFD.cif`, `8EFE.cif`, and `8EFI.cif` from the exact RCSB URLs in the feature manifest; decompress the 8ACT archive into the ignored cache path below. The P12883 FASTA fixture is already in this repository.

```bash
mkdir -p workflows/lsar-na-rel/cache workflows/myh7-correlation-screen/cache
curl -fL https://files.rcsb.org/download/8ACT.pdb1.gz -o workflows/lsar-na-rel/cache/8ACT.pdb1.gz
gzip -dc workflows/lsar-na-rel/cache/8ACT.pdb1.gz > workflows/lsar-na-rel/cache/8ACT-assembly1.pdb
curl -fL https://files.rcsb.org/download/8EFD.cif -o workflows/myh7-correlation-screen/cache/8EFD.cif
curl -fL https://files.rcsb.org/download/8EFE.cif -o workflows/myh7-correlation-screen/cache/8EFE.cif
curl -fL https://files.rcsb.org/download/8EFI.cif -o workflows/myh7-correlation-screen/cache/8EFI.cif
```

The `features` command checks every downloaded structure against the pinned
hashes in the [feature manifest](data/feature-manifest-v1.json).

```bash
conda run -n lsar-na-rel-talk-v0 python workflows/myh7-correlation-screen/screen.py features \
  --variants workflows/myh7-correlation-screen/variant-registry-v1.csv \
  --fasta tests/lsar_na_rel/fixtures/P12883.fasta \
  --assembly-8act workflows/lsar-na-rel/cache/8ACT-assembly1.pdb \
  --cif-8efd workflows/myh7-correlation-screen/cache/8EFD.cif \
  --cif-8efe workflows/myh7-correlation-screen/cache/8EFE.cif \
  --cif-8efi workflows/myh7-correlation-screen/cache/8EFI.cif \
  --output-dir workflows/myh7-correlation-screen/data

python3 workflows/myh7-correlation-screen/screen.py outcomes \
  --kcat-labels workflows/kcat-rel/data/canonical-labels.csv \
  --kcat-evidence workflows/kcat-rel/data/measurement-evidence.csv \
  --na-primary data/public/lsar-na-rel/na-rel-primary.csv \
  --na-labels data/public/lsar-na-rel/labels.csv \
  --na-measurements data/public/lsar-na-rel/measurements.csv \
  --velocity-summary workflows/myh7-correlation-screen/data/velocity-summary-evidence-v1.csv \
  --morck-pairs workflows/myh7-correlation-screen/data/morck-velocity-pairs-v1.csv \
  --output-dir workflows/myh7-correlation-screen/data

python3 workflows/myh7-correlation-screen/screen.py plan \
  --outcomes workflows/myh7-correlation-screen/data/outcomes-v1.csv \
  --new-features workflows/myh7-correlation-screen/data/features-v1.csv \
  --kcat-features workflows/kcat-rel/data/derived/kcat-rel-v0-no-msa-features.csv \
  --na-features data/public/lsar-na-rel/features.csv \
  --config workflows/myh7-correlation-screen/comparison-config-v1.json \
  --output-dir workflows/myh7-correlation-screen/data

python3 workflows/myh7-correlation-screen/screen.py analyze \
  --table workflows/myh7-correlation-screen/data/analysis-table-v1.csv \
  --inventory workflows/myh7-correlation-screen/data/comparison-inventory-v1.csv \
  --plan-manifest workflows/myh7-correlation-screen/data/comparison-plan-manifest-v1.json \
  --outcomes workflows/myh7-correlation-screen/data/outcomes-v1.csv \
  --dependence-ledger workflows/myh7-correlation-screen/data/dependence-ledger-v1.csv \
  --velocity-alternatives workflows/myh7-correlation-screen/data/velocity-alternatives-v1.csv \
  --output-dir workflows/myh7-correlation-screen/data
```

The plan command records `association_estimates_computed=false`; the analyze command verifies the saved plan's table, inventory, and outcome hashes before calculating. Rerunning manifests updates their timestamps, while the CSV tables remain deterministic for fixed inputs. Run `conda run -n lsar-na-rel-talk-v0 pytest -q workflows/myh7-correlation-screen/tests` for the command-line output checks.

The optional publisher-workbook extraction is available as
`python workflows/myh7-correlation-screen/extract_morck_pairs.py --workbook <publisher-workbook.xlsx> --output <pairs.csv>`.
The committed pair table is the input for reproducing the published estimates;
the original workbook is not part of this repository.

## Interpretation boundary

Primary estimates are descriptive Spearman rho on raw ratios, or median group differences for binary annotations. Natural-log ratio Pearson is secondary. No p-values or multiplicity-adjusted discovery claims are made. Per-study and per-variant deletion ranges are sensitivity descriptions, not population confidence intervals. The static WT structures are annotation templates, not measured mutant effects. Assay heterogeneity, shared WT values, exact reuse of seven short-head ATPase components across `kcat_rel` and `Na_rel`, and sparse binary groups limit inference. The [new validation protocol](../../docs/matched-atpase-velocity-protocol.md) targets the most actionable measurement question.
