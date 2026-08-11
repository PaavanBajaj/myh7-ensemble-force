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

**S1 — Provenance-ready.** Continuous LSAR public derivatives and source
manifest are published under `data/public/lsar-na-rel/` and
`data/public/manifests/lsar-na-rel-sources.json`. This is not S2: no pinned,
runnable regeneration pipeline is shipped in this release.
