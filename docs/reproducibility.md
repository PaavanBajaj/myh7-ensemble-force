# Reproducibility

Stage ladder defined for the original predictive program. “Reproducible” means
what a third party can do from the public tree (plus vault evidence they already
hold). S4 and S5 describe historical planned milestones, not current goals.

| Stage | Name | Required public content | “Reproducible” means |
| --- | --- | --- | --- |
| **S0** | Scaffold | README lineage, LICENSE, env stub, empty `docs/` stubs | Clone + read intent |
| **S1** | Provenance-ready | `provenance.md` + scrubbed schema/manifest + shareable derivatives | Vault holder can locate cited sources from public IDs/checksums |
| **S2** | Label workflow live | `src/<label>/`, `tests/<label>/`, `data/public/<label>/`, and `results/<label>/` with a runbook and locked env | Regenerate features/splits/baseline metrics from public recipe |
| **S3** | Label results shipped | `results/<label>/results/` pack | Third party regenerates the pack |
| **S4** | Ensemble-force shipped | `results/ensemble-force/` method + results | Same for the combined estimate |
| **S5** | v1 package | In-scope labels + ensemble + decisions + pinned env | Declared v1 path end-to-end |

## Current stage

**S3 — Label results shipped** for `lsar-na-rel` and the frozen `kcat_rel`
diagnostic. The public tree contains
repository-local inputs and feature provenance, tested feature/model commands,
direct dependency pins, an exact SHA-256 `osx-arm64` lock, and regenerated
result packs. The near-baseline LSAR prediction result does not establish a
useful predictor. The published `kcat_rel` diagnostic is likewise a negative
result. The correlation screen is a descriptive analysis with committed inputs,
all 78 results, and CLI commands to regenerate its tables. The Morck source
workbook and restricted literature remain outside Git; the compact derived
pair table is included. No ensemble-force workflow was implemented.

All completed workflows use the same public layout: executable modules in
`src/`, tests in `tests/`, shareable inputs and outputs in `data/public/`,
methods and evidence decisions in `docs/`, and runbooks, configurations, and
frozen run records in `results/`. Downloaded structures live in ignored
`.cache/<workflow>/` directories. The commands
in each workflow README run from the repository root. The historical `kcat_rel`
diagnostic is preserved as a frozen record; its feature and evidence audits
remain runnable, but its diagnostic should not be replayed into the frozen
result directory.

The explicit lock is currently released only for `osx-arm64`.
