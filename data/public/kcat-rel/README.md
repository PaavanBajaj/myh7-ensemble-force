# Data dictionary

## `canonical-labels.csv`

One modeling label per frozen primary variant.

| Field | Meaning |
| --- | --- |
| `label_id` | Stable canonical-label identifier |
| `variant`, `gene` | MYH7 amino-acid substitution identity |
| `source_id` | Leakage group used by study-grouped validation |
| `canonical_measurement_ids` | `|`-delimited links to evidence rows |
| `construct_family`, `construct` | Selected accepted short-motor family and reported construct |
| `aggregation_method` | Ratio of summary means or published average of paired replicate ratios |
| `wt_kcat_s_1`, `mutant_kcat_s_1` | Source rates; `|` preserves paired replicate values |
| `wt_kcat_error`, `mutant_kcat_error`, `error_type` | Reported source uncertainty |
| `kcat_rel`, `kcat_rel_error` | Canonical relative turnover and propagated uncertainty |
| `ln_kcat_rel`, `ln_kcat_rel_error` | Natural-log modeling label and delta-method uncertainty |
| `uncertainty_method` | How uncertainty was carried to the canonical label |
| `selection_basis` | Cohort-consistent construct/assay choice |

## `measurement-evidence.csv`

The 22 source measurements underlying the 17 labels. Summary measurements stay
as summary rows; the ten Morck measurements stay as two biological replicates
for each of five variants.

| Field group | Fields |
| --- | --- |
| Identity | `measurement_id`, `variant`, `gene`, `source_id`, `biological_replicate` |
| Eligibility | `canonical_evidence`, `construct_family`, `construct`, `actin_system`, `assay_regime` |
| Assay stratum | `detection_chemistry`, `temperature_c` |
| Raw values | WT/mutant `kcat` values, errors, and `error_type` |
| Published paired ratios | `reported_kcat_rel`, `reported_kcat_rel_error` (Morck replicate rows only) |
| WT pairing | `wt_context_id`, `mutant_context_id`, `wt_context_match` |
| Provenance | `source_location`, `doi`, `evidence_epoch`, `notes` |

Empty `reported_kcat_rel` fields mean the canonical ratio is recomputed from
the reported WT and mutant summary means; they do not mean the label is missing.

## Mathematical definitions

The CSV column names use ASCII identifiers, while the corresponding quantities
are:

```math
k_{\mathrm{cat,rel}}
= \frac{k_{\mathrm{cat,mut}}}{k_{\mathrm{cat,WT}}},
\qquad
y
= \ln\!\left(k_{\mathrm{cat,rel}}\right).
```

The `ln_kcat_rel` column stores $y$.

For independent WT and mutant summary errors, `ln_kcat_rel_error` is

```math
\sigma_{\ln k_{\mathrm{cat,rel}}}
= \sqrt{
    \left(\frac{\sigma_{\mathrm{mut}}}{k_{\mathrm{cat,mut}}}\right)^2
    +
    \left(\frac{\sigma_{\mathrm{WT}}}{k_{\mathrm{cat,WT}}}\right)^2
  },
```

and `kcat_rel_error` is
$k_{\mathrm{cat,rel}}\,\sigma_{\ln k_{\mathrm{cat,rel}}}$. For the five
Morck variants, the published paired-ratio aggregate and its retained
uncertainty replace the independent-summary calculation.
