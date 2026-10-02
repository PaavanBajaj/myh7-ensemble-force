# Matched MYH7 ATPase and unloaded velocity validation

**Date:** 2026-10-02  
**Status:** Proposed focused next experiment, selected after the exploratory screen. No measurements or production structural simulations are authorized by this document alone.  
**Decision basis:** [2026-10-02 screen decision](myh7-correlation-screen-decision.md).

## Question

Does the positive across-variant rank relationship between actin-activated `kcat_rel` and unloaded pure-actin `v_rel` persist when both outcomes are measured from matched human β-cardiac MYH7 short-head preparations, against concurrently measured WT, under one shared assay regime?

The existing 14-variant exploratory rho is +0.670, but combines different constructs, velocity metrics, and studies. Several individual variants disagree with a simple monotonic relationship. This protocol tests reproducibility of the **association and measurements**; it does not presume that one endpoint causes the other.

## Panel and measurements

Use human β-cardiac MYH7 sS1 (residues 1–808) with the same light-chain and attachment configuration for WT and **D239N, H251N, Y115H, E497D, R403Q, R663H, and G768R**. The seven variants span the observed velocity range and include E497D as a discordant case. Sequence-verify each construct. Retain G768R's preprint evidence status in reporting. Do not swap in a different variant after seeing new data without recording a new protocol version.

For each variant and a concurrently prepared WT reference, target at least three **independent protein preparations**, each with paired measurements of (1) actin-activated ATPase saturation `kcat` in s⁻¹ and (2) unloaded pure F-actin motility FAST `MVEL20` in nm/s. Use a common 23 °C, 25 mM KCl assay regime when both methods permit it; record the full ATPase and motility buffers, actin source, ATP concentration, immobilization density, light-chain occupancy, slide selection rule, and temperature. Prespecify the motor-density plateau and motility tracking/filtering rules before seeing variant outcomes. Record technical repeats and filament tracks beneath preparation IDs, not as independent biological samples. Randomize run order and mask variant identity during FAST analysis if feasible.

Compute per-preparation mutant/paired-WT ratios for both endpoints; aggregate at the preparation level, preserving absolute WT and mutant values and error types. When a batch shares a WT preparation across mutants, identify that common denominator explicitly. Do not divide by an external or global WT. Retain every exclusion and technical failure with its prespecified reason.

## Endpoints and interpretation

The primary descriptive endpoint is tie-aware Spearman rho between the seven variant-level median `kcat_rel` and `v_rel` values, with a leave-one-variant range. Report the individual variant ratios and paired preparation plots first. Also report the sign and range after leaving out a shared WT batch, and a preparation-level resampling interval that resamples whole preparations or WT batches rather than filaments. Because this is a seven-variant pilot, the interval describes measurement variability under the panel and does not establish a population-wide relationship. Do not present an unadjusted exploratory p-value as confirmation.

Treat the positive screen pattern as **replicated within this panel** only if at least six of seven variants have at least three usable paired preparations, the aggregate rho is positive, and every leave-one-variant estimate remains positive. A failure of these conditions leaves the question unresolved or contradicted; inspect measurement quality and heterogeneity without replacing the primary analysis. The discordant E497D result should be shown explicitly even if the overall sign is positive. Report absolute and relative endpoints, study/batch dependence, and all null or discordant cases.

Secondary checks may compare unfiltered `MVEL` with `MVEL20`, and assess whether apparent association depends on G768R or one shared WT batch. These are labeled sensitivity analyses. An ADP-release, optical-trap step-size, or PPS/IHM experiment needs a separate question and protocol informed by the paired measurements.

## Scope and stewardship

This is a proposed experimental validation design. It does not reopen the frozen `kcat_rel` or public `Na_rel` diagnostics, generate a velocity predictor, infer combined force, or assert a structural mechanism. Raw traces, fit outputs, and literature full text stay in the appropriate private vault; a later versioned evidence table can cite them with row IDs, context, and checksums. Future raw-data releases require a separate rights and provenance review.
