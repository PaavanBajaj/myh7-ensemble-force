"""Official Talk-v0 labeled table: variant, model features, grouping, y.

``ihm_min_dist_A`` is the audited minimum heavy-atom distance (Å) to an IHM
partner (min of BH–FH and BH–S2). The binary ``ihm_interface_flag_8act`` is 1
iff that distance is ≤ 4.5 Å. Distances are joined from the existing sparse-flag
audit CSV; 8ACT geometry is not regenerated.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from lsar_na_rel.paths import DEFAULT_PATHS
from lsar_na_rel.structure import CONTACT_CUTOFF_ANGSTROM

IHM_MIN_DIST_COLUMN = "ihm_min_dist_A"

MODEL_TABLE_COLUMNS: tuple[str, ...] = (
    "variant",
    "y",
    "ln_na_rel_public",
    "source_id",
    "normalization_tier",
    "delta_charge",
    "wt_site_rsa_8act_mean_ab",
    "ihm_interface_flag_8act",
    IHM_MIN_DIST_COLUMN,
)


def build_talk_v0_model_table(
    *,
    primary_path: str | Path,
    labels_path: str | Path,
    features_path: str | Path,
    distance_path: str | Path,
) -> pd.DataFrame:
    primary = pd.read_csv(primary_path)
    labels = pd.read_csv(labels_path)
    features = pd.read_csv(features_path)
    distances = pd.read_csv(distance_path)

    if "variant" not in primary.columns:
        raise ValueError("primary table missing variant")
    for frame, name in (
        (labels, "labels"),
        (features, "features"),
        (distances, "distance audit"),
    ):
        if "variant" not in frame.columns:
            raise ValueError(f"{name} missing variant")
        if frame["variant"].duplicated().any():
            raise ValueError(f"{name} has duplicate variant keys")

    if IHM_MIN_DIST_COLUMN not in distances.columns:
        raise ValueError(f"{distance_path} missing {IHM_MIN_DIST_COLUMN}")

    na_rel = pd.to_numeric(primary["Na_rel"], errors="coerce").to_numpy(dtype=float)
    public_ln = pd.to_numeric(primary["ln(Na_rel)"], errors="coerce").to_numpy(dtype=float)
    y = np.log(na_rel)
    if np.max(np.abs(y - public_ln)) >= 1e-10:
        raise ValueError("log(Na_rel) disagrees with public ln(Na_rel)")

    variants = primary["variant"].astype(str)
    base = pd.DataFrame(
        {
            "variant": variants,
            "y": y,
            "ln_na_rel_public": public_ln,
        }
    )
    label_cols = labels.loc[:, ["variant", "source_id", "normalization_tier"]].copy()
    label_cols["variant"] = label_cols["variant"].astype(str)
    feat_cols = features.loc[
        :,
        [
            "variant",
            "delta_charge",
            "wt_site_rsa_8act_mean_ab",
            "ihm_interface_flag_8act",
        ],
    ].copy()
    feat_cols["variant"] = feat_cols["variant"].astype(str)
    dist_cols = distances.loc[:, ["variant", IHM_MIN_DIST_COLUMN]].copy()
    dist_cols["variant"] = dist_cols["variant"].astype(str)

    table = (
        base.merge(label_cols, on="variant", how="left", validate="one_to_one")
        .merge(feat_cols, on="variant", how="left", validate="one_to_one")
        .merge(dist_cols, on="variant", how="left", validate="one_to_one")
    )
    if table[list(MODEL_TABLE_COLUMNS)].isna().any().any():
        missing = table.loc[table.isna().any(axis=1), "variant"].tolist()
        raise ValueError(f"Missing values in Talk-v0 model table for: {missing}")

    dist = pd.to_numeric(table[IHM_MIN_DIST_COLUMN], errors="coerce").to_numpy(dtype=float)
    flag = pd.to_numeric(table["ihm_interface_flag_8act"], errors="coerce").to_numpy(
        dtype=int
    )
    expected_flag = (dist <= float(CONTACT_CUTOFF_ANGSTROM)).astype(int)
    if not np.array_equal(flag, expected_flag):
        raise ValueError(
            "ihm_interface_flag_8act does not match ihm_min_dist_A ≤ "
            f"{CONTACT_CUTOFF_ANGSTROM} Å"
        )

    return table.loc[:, list(MODEL_TABLE_COLUMNS)].reset_index(drop=True)


def write_talk_v0_model_table(path: str | Path, table: pd.DataFrame) -> None:
    ordered = table.loc[:, list(MODEL_TABLE_COLUMNS)]
    Path(path).write_text(ordered.to_csv(index=False), encoding="utf-8")


def main() -> None:
    table = build_talk_v0_model_table(
        primary_path=DEFAULT_PATHS.primary_path,
        labels_path=DEFAULT_PATHS.labels_path,
        features_path=DEFAULT_PATHS.features_path,
        distance_path=DEFAULT_PATHS.public_data / "ihm_sparse_flag_variant_summary.csv",
    )
    out = DEFAULT_PATHS.public_data / "talk-v0-model-table.csv"
    write_talk_v0_model_table(out, table)
    print(f"wrote {out} rows={len(table)}")


if __name__ == "__main__":
    main()
