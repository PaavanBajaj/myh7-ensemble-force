"""Deterministic canonical S1 input loading, schema audit, and label join."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

EXPECTED_PRIMARY_HEADER = "variant,LSAR_mut,LSAR_WT,Na_rel,ln(Na_rel)"
EXPECTED_PRIMARY_COLUMNS = tuple(EXPECTED_PRIMARY_HEADER.split(","))

EXPECTED_SOURCE_IDS = frozenset(
    {
        "adhikari2019",
        "morck2022",
        "nandwani2025",
        "pathak2026",
        "sarkar2020",
        "vanderroest2021",
    }
)

RISK_ALLOWED_VARIANTS = frozenset({"R403Q", "R663H", "G768R"})

# Columns that must never become estimator inputs (group/clinical/evidence metadata).
FORBIDDEN_FEATURE_COLUMNS = frozenset(
    {
        "source_id",
        "gene",
        "evidence_class",
        "modeling_eligibility",
        "modeling_ready",
        "canonical_lsar",
        "canonical_lsar_basis",
        "reported_lsar",
        "reported_lsar_error",
        "recomputed_lsar",
        "table1_na_percent",
        "table1_cell_state",
        "table1_implied_lsar_audit_only",
        "spudich_additional_na_percent_uncapped",
        "spudich_additional_na_percent_capped_display",
        "na_rel_global_spudich057",
        "study_wt_lsar",
        "study_wt_lsar_error",
        "na_rel",
        "na_rel_error",
        "ln_na_rel",
        "ln_na_rel_error",
        "release_audit_only",
        "ln_release_audit_only",
        "study_normalization_percent",
        "short_construct",
        "long_construct",
        "rlc_phosphorylation_state",
        "assay_method",
        "low_biological_n",
        "evidence_epoch",
        "normalization_tier",
        "notes",
        "paper_id",
        "pmid",
        "doi",
        "clinvar",
        "pathogenicity",
        "phenotype",
        "y",
        "ln_na_rel_public",
        "ln(Na_rel)",
        "LSAR_mut",
        "LSAR_WT",
        "Na_rel",
        "variant",
    }
)

_VARIANT_RE = re.compile(r"^(?P<wt>[A-Z])(?P<position>\d+)(?P<mut>[A-Z])$")


@dataclass(frozen=True)
class ParsedVariant:
    label: str
    wt: str
    position: int
    mut: str


@dataclass(frozen=True)
class CohortBundle:
    """Canonical cohort with leakage-resistant separations.

    ``feature_column_names`` records the selected profile schema only.
    ``feature_values`` stays empty until later feature-generation issues.
    Target, groups, and eligibility never share that feature namespace.
    """

    primary: pd.DataFrame
    target: pd.Series
    groups: pd.Series
    eligibility: pd.DataFrame
    feature_column_names: tuple[str, ...]
    feature_values: pd.DataFrame | None = None


def parse_variant(label: str) -> ParsedVariant:
    match = _VARIANT_RE.fullmatch(str(label).strip())
    if match is None:
        raise ValueError(f"Invalid variant label: {label!r}")
    return ParsedVariant(
        label=match.group(0),
        wt=match.group("wt"),
        position=int(match.group("position")),
        mut=match.group("mut"),
    )


def _read_primary(primary_path: Path) -> pd.DataFrame:
    header = primary_path.read_text(encoding="utf-8").splitlines()[0]
    if header != EXPECTED_PRIMARY_HEADER:
        raise ValueError(
            f"Unexpected primary header in {primary_path}: {header!r}; "
            f"expected {EXPECTED_PRIMARY_HEADER!r}"
        )
    frame = pd.read_csv(primary_path)
    if tuple(frame.columns) != EXPECTED_PRIMARY_COLUMNS:
        raise ValueError(f"Unexpected primary columns: {tuple(frame.columns)!r}")
    return frame


def _audit_primary(primary: pd.DataFrame) -> pd.DataFrame:
    if len(primary) != 16:
        raise ValueError(f"Expected 16 primary rows, found {len(primary)}")
    variants = primary["variant"].astype(str)
    if variants.duplicated().any():
        dupes = variants[variants.duplicated()].tolist()
        raise ValueError(f"Duplicate variants in primary table: {dupes}")
    if "G256E" in set(variants):
        raise ValueError("G256E must be absent from the canonical primary cohort")

    for label in variants:
        parse_variant(label)

    na_rel = pd.to_numeric(primary["Na_rel"], errors="coerce")
    public_ln = pd.to_numeric(primary["ln(Na_rel)"], errors="coerce")
    if not np.isfinite(na_rel).all():
        raise ValueError("Na_rel must be finite numeric")
    if (na_rel <= 0).any():
        raise ValueError("Na_rel must be strictly positive")
    if not np.isfinite(public_ln).all():
        raise ValueError("ln(Na_rel) must be finite numeric")

    recomputed = np.log(na_rel.to_numpy(dtype=float))
    delta = np.max(np.abs(recomputed - public_ln.to_numpy(dtype=float)))
    if delta >= 1e-10:
        raise ValueError(
            f"Recomputed log(Na_rel) disagrees with public ln(Na_rel); max |delta|={delta}"
        )

    audited = primary.copy()
    audited["Na_rel"] = na_rel.astype(float)
    audited["ln(Na_rel)"] = public_ln.astype(float)
    audited["LSAR_mut"] = pd.to_numeric(audited["LSAR_mut"], errors="raise")
    audited["LSAR_WT"] = pd.to_numeric(audited["LSAR_WT"], errors="raise")
    return audited.reset_index(drop=True)


def _join_labels(primary: pd.DataFrame, labels_path: Path) -> pd.DataFrame:
    labels = pd.read_csv(labels_path)
    if "variant" not in labels.columns:
        raise ValueError("labels.csv missing variant column")
    if labels["variant"].duplicated().any():
        raise ValueError("labels.csv has duplicate variant keys")

    joined = primary.merge(labels, on="variant", how="left", validate="one_to_one")
    if joined["source_id"].isna().any():
        missing = joined.loc[joined["source_id"].isna(), "variant"].tolist()
        raise ValueError(f"Primary variants missing from labels.csv: {missing}")
    return joined


def _audit_join(joined: pd.DataFrame) -> None:
    groups = set(joined["source_id"].astype(str))
    if groups != EXPECTED_SOURCE_IDS:
        raise ValueError(
            f"Unexpected source_id set: {sorted(groups)}; expected {sorted(EXPECTED_SOURCE_IDS)}"
        )

    tiers = joined["normalization_tier"].astype(str)
    primary_n = int((tiers == "Primary").sum())
    risky_n = int((tiers == "Risky allowed").sum())
    if primary_n != 13 or risky_n != 3:
        raise ValueError(
            f"Expected 13 Primary + 3 Risky allowed rows; found {primary_n} + {risky_n}"
        )

    risky = set(joined.loc[tiers == "Risky allowed", "variant"].astype(str))
    if risky != RISK_ALLOWED_VARIANTS:
        raise ValueError(
            f"Risky allowed set mismatch: {sorted(risky)}; "
            f"expected {sorted(RISK_ALLOWED_VARIANTS)}"
        )


def _validate_feature_columns(feature_column_names: Sequence[str]) -> tuple[str, ...]:
    columns = tuple(str(name) for name in feature_column_names)
    if not columns:
        raise ValueError("feature_column_names must be non-empty (select a profile first)")
    if len(set(columns)) != len(columns):
        raise ValueError("feature_column_names must be unique")
    leaked = set(columns) & FORBIDDEN_FEATURE_COLUMNS
    if leaked:
        raise ValueError(f"Forbidden columns selected as features: {sorted(leaked)}")
    return columns


def load_canonical_cohort(
    *,
    primary_path: str | Path,
    labels_path: str | Path,
    feature_column_names: Sequence[str],
) -> CohortBundle:
    """Load and audit the public S1 cohort without generating feature values.

    Callers must select an approved profile (and therefore
    ``feature_column_names``) before loading the target via this function.
    """
    columns = _validate_feature_columns(feature_column_names)
    primary = _audit_primary(_read_primary(Path(primary_path)))
    joined = _join_labels(primary, Path(labels_path))
    _audit_join(joined)

    target = pd.Series(
        np.log(primary["Na_rel"].to_numpy(dtype=float)),
        index=primary.index,
        name="y",
    )
    groups = joined["source_id"].astype(str).rename("source_id")
    eligibility = joined[["normalization_tier"]].copy()

    return CohortBundle(
        primary=primary,
        target=target,
        groups=groups,
        eligibility=eligibility,
        feature_column_names=columns,
        feature_values=None,
    )
