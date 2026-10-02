# Public workflow promotion audit, 2026-10-02

The public repository now contains all completed MYH7 workflows: the existing
LSAR release, the frozen 17-variant `kcat_rel` diagnostic, and the exploratory
correlation screen. No velocity predictor or combined-force workflow was built.

## Value preservation

- The published `kcat_rel` canonical labels, measurement evidence, feature CSV,
  model table, out-of-fold predictions, metrics, and sensitivity data are
  byte-identical to the frozen private workflow copies.
- The public correlation `outcomes` command reads the published `kcat_rel` and
  LSAR inputs plus the committed velocity source tables. Its `plan` and
  `analyze` commands reproduce all 78 comparison results, 1,227 plotted-point
  rows, and 125 sensitivity rows. The public and private result rows compare
  field for field with no numeric or categorical differences.
- The 78 SVG plots are byte-identical to the original screen. The ranking CSV
  excludes the binary IHM flag and retains the original 66 Spearman rows.

## Publication edits

Private absolute paths and sibling-checkout references were replaced with
repository-local paths. Vault paths and hashes of private vault-only files were
removed from the public source manifest. Restricted full text, PDFs, and the
original Morck workbook were not copied. The publicly downloadable Morck
workbook's checksum remains in the extraction script so its identity can be
verified. The compact source-linked Morck pair table remains a published input.

The frozen `kcat_rel` feature manifest and diagnostic run manifest were
path-scrubbed and rehashed. Their original execution revision and the
post-execution reconstruction account remain labeled as historical. The
diagnostic was not refit or rerun during publication. Rebuilt correlation
manifests use public input hashes and new packaging timestamps; the scientific
result tables remain unchanged.

Numerical source values were not re-audited or revised in this promotion.
The data dictionaries retain units, cohort definitions, uncertainty types,
missingness, and study provenance. The public release adds no separate data
license beyond the repository's existing license statement.
