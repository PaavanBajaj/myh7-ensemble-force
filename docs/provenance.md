# Provenance

How this public repository cites evidence. Vault evidence is
**cite, do not contain**: literature PDFs, full-text extracts, restricted
supplements, and raw working data stay in the encrypted private data vault.
Public materials point at scrubbed identifiers, checksums, and freeze
metadata only.

## Status (S1)

First public provenance pack for the continuous LSAR slice:

- Derivatives: `data/public/lsar-na-rel/`
- Source manifest: `data/public/manifests/lsar-na-rel-sources.json`
- Method account: `workflows/lsar-na-rel/method.md`

The manifest cites primary sources by DOI and, where recorded, PMID/PMCID, plus
SHA-256 digests of the published derivative files themselves.

## Rules

1. Cite sources by stable IDs and checksums published in public manifests.
2. Do not commit copyrighted full text, screenshots of papers, or vault paths
   that imply redistribution.
3. Scrubbed shareable derivatives live under `data/public/` only after
   eligibility and copyright review.
4. Do not assign a new data license unless one is explicitly supplied; document
   the current licensing limitation beside the pack.
5. Hash the public derivative files themselves; do not publish private vault
   file hashes as if they were the public pack.
