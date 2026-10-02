# Decisions and limitations

Living log for published methodology decisions and standing limitations.

## Project direction (2026-10-02)

The active direction is descriptive analysis of measured MYH7 biochemical
effects and investigation of possible variant mechanisms. The original
predictive ensemble-force plan is historical. No velocity model, forward
prediction model, or combined force calculation was implemented. The completed
ATPase diagnostic failed its prediction gates and is published alongside the
LSAR S3 release. The exploratory correlation screen is also published, with a
dated decision and protocol for matched-preparation measurements. These are
negative predictive and descriptive results, respectively.

The three completed workflows share one repository layout: code in `src/`,
tests in `tests/`, public data in `data/public/`, and method/configuration and
frozen run records in `workflows/`. See the [workflow index](workflows.md) for
the exact directories.

## Status (S3 for `lsar-na-rel`)

Promoted decisions for the continuous LSAR public pack
(`data/public/lsar-na-rel/`):

1. Modeled biochemical input is continuous LSAR
   (\(k_{\mathrm{cat,long}}/k_{\mathrm{cat,short}}\)). Canonical
   \(N_{a,\mathrm{rel}}=\mathrm{LSAR}_{mut}/\mathrm{LSAR}_{WT,\mathrm{study}}\);
   GPR target is \(y=\ln(N_{a,\mathrm{rel}})\). Primary pipeline table is
   `data/public/lsar-na-rel/na-rel-primary.csv`.
2. Spudich additional-\(N_a\), global \(\mathrm{LSAR}/0.57\), release
   \(a=(\mathrm{LSAR}_{mut}-\mathrm{LSAR}_{WT})/(1-\mathrm{LSAR}_{WT})\),
   and Table-implied reverse LSAR are stored as distinct derived/audit fields;
   release and reverse LSAR are not modeling targets.
3. One canonical label row per variant; replicates belong only in measurements.
4. Reported LSAR is preserved separately from recomputed long/short; neither is
   silently replaced to force Table 1 agreement.
5. Raw LSAR and within-study labels are not capped; display caps are recorded
   separately.
6. G256E remains approximate and excluded from exact-label modeling.
7. Q222K, E536D, and V606M remain unresolved unpublished.
8. D239N, R453C, and G741R are assay-absent blanks, not zeroes.
9. G768R’s later Pathak measurement is canonical and separate from Spudich’s
   unpublished 2024 epoch.
10. R403Q / R663H remain low-biological-\(n\) evidence (`review_low_n`).
11. Preprint/final versions and repeated presentations of the same summary are
    one evidence identity.
12. No new outside-Table-1 exact LSAR variant was recovered in the audited
    expansion screen summarized with this pack.
13. No separate public data license is declared for the compiled derivative
    beyond the repository MIT license’s ordinary scope (see pack README).
14. Feature generation reads variant identities, the P12883 sequence, the
    checksum-pinned 8ACT assembly, and `lsar-na-rel.json`; it never loads target
    or source-group columns as estimator inputs.
15. The registered ladder contains exactly three run IDs. Exploratory distance
    probes are not part of the public runnable workflow.
16. `ihm_interface_flag_8act` uses a frozen 4.5 Å minimum-heavy-atom cutoff and
    is positive for D382Y and R403Q only in this snapshot.
17. Validation is leave-one-`source_id`-out, with training-fold-only scaling and
    a training-fold mean baseline. The fixed seed is 20260811.
18. The GPR kernel, optimizer restarts, output artifact contract, run IDs, and
    feature order are all loaded from one validated configuration.
19. No FoldX-derived values are present. The registered public model uses
    `delta_charge` plus checksummed 8ACT features.
20. Exact environment locking is released for `osx-arm64` only.

## Limitations (standing)

- Research-use only; not a diagnostic test; not clinical advice.
- Vault evidence is cite-not-contain — literature PDFs and full-text extracts
  do not live in this repository.
- The LSAR result is a small, heterogeneous, 16-row descriptive analysis. Its
  best pooled OOF MAE is 0.130886 versus 0.131202 for the training-fold mean
  baseline. This near-baseline result does not establish useful prediction,
  clinical utility, mechanism, or causality.
- The 8ACT interface flag describes static WT structural proximity. It does not
  measure blocked/free-head occupancy or pre-power-stroke stability. Biochemical
  LSAR, structural contacts, and kinetic state readouts require separate
  interpretation.
- Scikit-learn may emit optimizer convergence warnings under the frozen setup;
  they are recorded as expected warnings and do not change the fitted contract.
- Construct family, assay chemistry, study WT baseline, and low-\(n\) rows
  still need pooling review before any fitted model.
- The public `kcat-rel` diagnostic failed its predeclared prediction gates.
  Its historical results are published without rerunning or refitting them.
- The exploratory correlation screen contains 78 comparisons across small,
  heterogeneous cohorts. Its rank associations are descriptive, not causal or
  confirmatory. The [decision](myh7-correlation-screen-decision.md) selects
  matched-preparation ATPase and velocity validation as the next question.
- No velocity predictor or ensemble-force workflow was implemented.
