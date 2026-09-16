# 8ACT sparse IHM-flag debugging report

**Date:** 2026-08-12
**Scope:** Diagnose why only 2/16 issue #28 variants have
`ihm_interface_flag_8act = 1`.
**Decision:** Retain the prespecified 4.5 Å primary cutoff. No existing source,
test, configuration, or feature file was modified.

## Bottom line

The sparse flag is **not an implementation error**.

- An independent coordinate scan exactly reproduced every committed BH–FH,
  BH–S2, and main 4.5 Å flag.
- The implementation uses the intended minimum non-hydrogen atom distance,
  not Cα, centroid, Cβ, or average distance.
- BH=A, FH=B, motor/lever=3–837, and proximal S2=838–chain end are applied
  correctly; chains C–F are excluded from the main flag.
- The two positives, `D382Y` and `R403Q`, are robust at both 4.0 and 4.5 Å and
  agree with the direct head–head interface described by Grinzato et al.
- Of the 14 negatives, 11 remain negative through 6.0 Å. Three become positive
  between 4.5 and 6.0 Å: `H251N`, `R249Q`, and `R719W`.
- Five variants make ≤4.5 Å contacts only to light chains
  (`D778V`, `L781P`, `S782N`, `A797T`, `F834L`). Counting “any different
  polymer” would inflate the main flag from 2 to 7 and would reintroduce the
  obsolete biological definition.

The correct interpretation is that this binary feature is a descriptor of
**direct WT-site heavy-chain IHM contact in one static 8ACT conformation**. It
is not a complete indicator of whether a mutation can alter IHM stability or
the accessible-head fraction through allostery.

## Reproduction and feedback loop

The exact symptom was reproduced from the committed feature snapshot:

```text
SPARSE_FLAG_REPRODUCED positives=['D382Y', 'R403Q'] count=2/16
```

The independent scan then computed all atom-pair distances directly from the
checksum-pinned `8ACT-assembly1.pdb` and compared classifications against:

1. `data/public/lsar-na-rel/features.csv`;
2. `ihm_partner_contacts()` at 4.0, 4.5, 5.0, 5.5, and 6.0 Å.

All comparisons matched:

```text
INDEPENDENT_DISTANCE_AUDIT_OK
```

No target, label, `Na_rel`, FoldX value, or model score was loaded or used.

## Structural mapping audited

- MYH7 heavy chains: A and B
- Blocked head (BH): A
- Free head (FH): B
- Essential light chains: C and D
- Regulatory light chains: E and F
- Motor/lever partner set: residues 3–837 on A/B
- Proximal-S2 partner set: residues 838–906 on A and 838–898 on B
- Main flag: OR of BH–FH and BH–S2 subtype contacts
- Primary contact rule: minimum heavy-atom distance ≤4.5 Å

The BH/FH assignment is independently supported by the asymmetric contact:

- A:R403 NH1 ↔ B:Y455 OH = 2.678 Å
- reverse B:R403 ↔ A:Y455 = 88.218 Å

## Exact distance calculation

For each WT focal residue on each heavy-chain copy:

```text
min_distance =
    min Euclidean distance(
        every non-H/non-D focal-residue atom,
        every non-H/non-D atom in the partner segment
    )
```

The per-variant value is the minimum across focal copies A and B for each
partner class. The main binary is the OR of the two subtype comparisons.

All 32 focal residues (16 positions × A/B) contained the expected standard
heavy-atom count, so no positive was lost because a focal side chain was
truncated in the deposited model. This does not eliminate coordinate
uncertainty: 8ACT is a 3.6 Å cryo-EM structure and the paper explicitly notes
labile “musical chairs” side chains and weaker local density.

## Sixteen-variant diagnostic summary

Distances are minimum heavy-atom distances in Å. `Light` is diagnostic only
and never enters the amended main flag.

| Variant | BH–FH | BH–S2 | IHM min | Light | RSA mean | 4.5 Å flag |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Y115H | 28.305 | 27.449 | 27.449 | 19.746 | 0.018018 | 0 |
| E497D | 10.695 | 18.675 | 10.695 | 20.342 | 0.087629 | 0 |
| R249Q | 5.994 | 5.097 | 5.097 | 17.183 | 0.205645 | 0 |
| H251N | 4.713 | 5.976 | 4.713 | 11.838 | 0.095109 | 0 |
| D382Y | 3.132 | 39.753 | 3.132 | 15.511 | 0.233129 | 1 |
| I457T | 10.315 | 12.305 | 10.315 | 10.999 | 0.000000 | 0 |
| R719W | 5.991 | 28.695 | 5.991 | 5.713 | 0.161290 | 0 |
| P710R | 15.667 | 22.028 | 15.667 | 17.217 | 0.110294 | 0 |
| D778V | 18.727 | 30.409 | 18.727 | 3.788 | 0.315951 | 0 |
| L781P | 18.122 | 31.053 | 18.122 | 3.241 | 0.152439 | 0 |
| S782N | 22.481 | 34.138 | 22.481 | 2.911 | 0.342308 | 0 |
| A797T | 29.229 | 31.045 | 29.229 | 3.136 | 0.028302 | 0 |
| F834L | 20.564 | 13.804 | 13.804 | 2.972 | 0.210660 | 0 |
| R403Q | 2.678 | 34.857 | 2.678 | 11.840 | 0.183468 | 1 |
| R663H | 8.165 | 6.372 | 6.372 | 20.458 | 0.608871 | 0 |
| G768R | 22.614 | 30.787 | 22.614 | 16.072 | 0.613095 | 0 |

Nearest-residue and atom-pair details for both structural copies are in
`ihm_sparse_flag_per_chain_diagnostics.csv`.

## Prespecified cutoff sensitivity

| Cutoff (Å) | Positive count | Positive variants |
| ---: | ---: | --- |
| 4.0 | 2 | D382Y; R403Q |
| 4.5 | 2 | D382Y; R403Q |
| 5.0 | 3 | H251N; D382Y; R403Q |
| 5.5 | 4 | R249Q; H251N; D382Y; R403Q |
| 6.0 | 5 | R249Q; H251N; D382Y; R719W; R403Q |

The increase is gradual, not a cliff at 4.5 Å:

- `H251N`: 4.713 Å BH–FH, only 0.213 Å beyond the primary cutoff.
- `R249Q`: 5.097 Å BH–S2, 0.597 Å beyond the primary cutoff.
- `R719W`: 5.991 Å BH–FH, 1.491 Å beyond the primary cutoff.

The first two are useful continuous-distance diagnostics. Their binary zeros
are still correct under the frozen 4.5 Å definition.

## Hypothesis results

1. **Genuine structural sparsity — supported.** Eleven of fourteen negatives
   remain beyond 6 Å from both approved heavy-chain partner classes.
2. **Near-threshold information loss — partly supported.** `H251N` and
   `R249Q` are close enough to merit continuous-distance provenance; the
   sensitivity curve does not justify changing the primary cutoff.
3. **Chain/segment mapping error — rejected.** A/B roles, Pro838 boundary, and
   chains C–F exclusion agree with the deposited assembly and public sources.
4. **Wrong distance calculation — rejected.** Independent minimum-heavy-atom
   distances reproduce production flags at every tested cutoff. Cα distances
   are substantially larger for direct contacts (for example, R403:
   2.678 Å heavy-atom versus 9.962 Å Cα).
5. **Missing focal side chains — rejected for these sites.** Every A/B focal
   residue has its expected heavy-atom count. Static-model and local-density
   uncertainty remain limitations, particularly for labile polar side chains.

## Biological cross-check

The 8ACT paper reports that interface interactions were generated with PDBsum,
manually checked, and extended to 4.5 Å for long-range weak interactions.
It also states that many cardiac-IHM side chains are dynamic and can occupy
alternative “musical chairs” conformations.

The coordinate results align with the paper’s mechanistic interpretation:

- R403 is a direct BH HCM-loop ↔ FH Y455 head–head lock.
- D382 is a direct BH PHHIS/head–head contact.
- H251 and R249 are functionally important Transducer/β-bulge mutations but
  are not equivalent to robust direct 4.5 Å locks in this static model.
- Pliant/lever mutations D778/L781/S782/A797/F834 lie near light-chain
  interfaces; those contacts are biologically interesting but deliberately
  excluded from the amended heavy-chain IHM flag.

One source-list correction to the supplied debugging brief: its PMC link and
Nature link both resolve to the same Grinzato et al. 2023 paper, not an
independent “later analysis.”

## Recommendation

1. Keep `ihm_interface_flag_8act` at the prespecified 4.5 Å cutoff.
2. Keep BH–FH and BH–S2 subtype bits as provenance.
3. Preserve the continuous per-chain and per-variant distances generated by
   this audit for interpretation and future prespecified analyses.
4. Do not add those distances to the current Talk-v0 model based on these
   outcomes; doing so requires a separate pre-model amendment.
5. Describe the binary as a direct-interface descriptor, not a complete
   mutation-mechanism classifier.
6. If uncertainty around H251/R249 matters, inspect density/alternative
   conformers or an independently specified structural ensemble rather than
   loosening the cutoff after observing sparsity.

## Generated artifacts

- `ihm_sparse_flag_debugging_report.md` — this report
- `ihm_sparse_flag_variant_summary.csv` — one row per variant
- `ihm_sparse_flag_per_chain_diagnostics.csv` — one row per A/B structural copy
- `ihm_sparse_flag_cutoff_sensitivity.csv` — prespecified cutoff counts

## Sources

1. Grinzato et al., *Nature Communications* 2023,
   https://doi.org/10.1038/s41467-023-38698-w
2. RCSB PDB 8ACT, https://doi.org/10.2210/pdb8act/pdb
3. Blankenfeldt et al., PDB 2FXM / Pro838 S2 boundary,
   https://doi.org/10.1073/pnas.0606741103
4. Adhikari et al., mutation/function context,
   https://pmc.ncbi.nlm.nih.gov/articles/PMC6582153/
