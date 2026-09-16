# Reproducibility

Stage ladder for the public repository. “Reproducible” means what a third
party can do from the public tree (plus vault evidence they already hold).

| Stage | Name | Required public content | “Reproducible” means |
| --- | --- | --- | --- |
| **S0** | Scaffold | README lineage, LICENSE, env stub, empty `docs/` stubs | Clone + read intent |
| **S1** | Provenance-ready | `provenance.md` + scrubbed schema/manifest + shareable derivatives | Vault holder can locate cited sources from public IDs/checksums |
| **S2** | Label workflow live | Live `workflows/<label>/` with `method.md` + pipeline + locked env | Regenerate features/splits/baseline metrics from public recipe |
| **S3** | Label results shipped | `workflows/<label>/results/` pack | Third party regenerates the pack |
| **S4** | Ensemble-force shipped | `workflows/ensemble-force/` method + results | Same for the combined estimate |
| **S5** | v1 package | In-scope labels + ensemble + decisions + pinned env | Declared v1 path end-to-end |

## Current stage

**S3 — Label results shipped** for `lsar-na-rel`. The public tree contains
repository-local inputs and feature provenance, tested feature/model commands,
direct dependency pins, an exact SHA-256 `osx-arm64` lock, and regenerated
result packs. Other label and ensemble-force workflows have not started.

The explicit lock is currently released only for `osx-arm64`.
