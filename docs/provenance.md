# Provenance

How this public repository cites evidence. Vault evidence is
**cite, do not contain**: literature PDFs, full-text extracts, restricted
supplements, and raw working data stay in the encrypted private data vault.
Public materials point at scrubbed identifiers, checksums, and freeze
metadata only.

## Public packs

First public provenance pack for the continuous LSAR slice:

- Derivatives: `data/public/lsar-na-rel/`
- Source manifest: `data/public/manifests/lsar-na-rel-sources.json`
- Method account: `docs/lsar-na-rel-method.md`
- Feature-generation provenance: `data/public/lsar-na-rel/features-provenance.json`

The manifest cites primary sources by DOI and, where recorded, PMID/PMCID, plus
SHA-256 digests of the published derivative files themselves.

The [`kcat-rel`](../data/public/kcat-rel/) pack includes canonical labels,
measurement evidence, feature provenance, and the historical diagnostic
manifest under [`results/kcat-rel/results/diagnostic-v0/`](../results/kcat-rel/results/diagnostic-v0/).
The [correlation pack](../data/public/myh7-correlation-screen/) includes a
source manifest, feature and outcome manifests, the locked comparison plan,
and the analysis manifest alongside all result tables and plots.

## Rules

1. Cite sources by stable IDs and checksums published in public manifests.
2. Do not commit copyrighted full text, screenshots of papers, or vault paths
   that imply redistribution.
3. Scrubbed shareable derivatives live under `data/public/<workflow>/` after
   provenance review.
4. Do not assign a new data license unless one is explicitly supplied; document
   the current licensing limitation beside the pack.
5. Hash the public derivative files themselves; do not publish private vault
   file hashes as if they were the public pack.

The public `kcat_rel` and correlation packs contain source-linked numerical
derivatives and executable code. Private machine paths and vault file hashes
were removed during publication. The frozen numerical values and correlation
estimates were preserved.

## Nonlinear dependence publication

The [new input plan](../data/public/myh7-nonlinear-dependence/) reuses the scrubbed public outcome pack, correlation input tables and historical result tables. Its input hashes refer to these public files. The [frozen output pack](../results/myh7-nonlinear-dependence/results/screen-v1/) adds scores, deletion checks and raw comparison points without redistributing raw source evidence. See the [publication audit](nonlinear-dependence-publication-audit-2026-10-02.md).
