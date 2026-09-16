"""Validated configuration for the public LSAR Na-rel workflow."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


APPROVED_LADDER_RUN_IDS: tuple[str, ...] = (
    "ladder-1-rsa",
    "ladder-2-rsa-ihm",
    "ladder-3-nofoldx",
)
SUPPORTED_KERNEL_TYPE = "ConstantKernel * Matern + WhiteKernel"


@dataclass(frozen=True)
class KernelSettings:
    type: str
    constant_value: float
    length_scale: float
    nu: float
    noise_level: float
    normalize_y: bool
    n_restarts_optimizer: int
    random_state: int


@dataclass(frozen=True)
class TalkV0Config:
    """Immutable, validated LSAR Na-rel workflow configuration."""

    raw: Mapping[str, Any]
    seed: int
    contact_cutoff_angstrom: float
    kernel: KernelSettings
    registered_run_ids: tuple[str, ...]
    required_artifacts: tuple[str, ...]

    @property
    def profiles(self) -> Mapping[str, list[str]]:
        return self.raw["profiles"]


def _require_mapping(payload: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = payload.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"LSAR Na-rel config {key!r} must be an object")
    return value


def _load_kernel(payload: Mapping[str, Any], seed: int) -> KernelSettings:
    raw = _require_mapping(payload, "kernel")
    kernel_type = str(raw.get("type", ""))
    if kernel_type != SUPPORTED_KERNEL_TYPE:
        raise ValueError(f"Unsupported kernel type: {kernel_type!r}")
    try:
        normalize_y = raw["normalize_y"]
        settings = KernelSettings(
            type=kernel_type,
            constant_value=float(raw["constant_value"]),
            length_scale=float(raw["length_scale"]),
            nu=float(raw["nu"]),
            noise_level=float(raw["noise_level"]),
            normalize_y=normalize_y,
            n_restarts_optimizer=int(raw["n_restarts_optimizer"]),
            random_state=int(raw["random_state"]),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Invalid kernel settings") from exc
    if not isinstance(settings.normalize_y, bool):
        raise ValueError("kernel.normalize_y must be boolean")
    if settings.nu != 1.5:
        raise ValueError("kernel.nu must be 1.5 for the registered workflow")
    if min(settings.constant_value, settings.length_scale, settings.noise_level) <= 0:
        raise ValueError("kernel scale values must be positive")
    if settings.n_restarts_optimizer < 0:
        raise ValueError("kernel.n_restarts_optimizer must be nonnegative")
    if settings.random_state != seed:
        raise ValueError("kernel.random_state must match top-level seed")
    return settings


def load_talk_v0_config(path: str | Path) -> TalkV0Config:
    config_path = Path(path)
    with config_path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"LSAR Na-rel config must be an object: {config_path}")
    if payload.get("target") != "y = log(Na_rel)":
        raise ValueError("target must be 'y = log(Na_rel)'")
    if payload.get("grouping_key") != "source_id":
        raise ValueError("grouping_key must be 'source_id'")
    try:
        seed = int(payload["seed"])
        cutoff = float(payload["contact_cutoff_angstrom"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("seed and contact_cutoff_angstrom must be numeric") from exc
    if cutoff <= 0:
        raise ValueError("contact_cutoff_angstrom must be positive")
    if not isinstance(payload.get("profiles"), dict):
        raise ValueError(f"LSAR Na-rel config missing profiles: {config_path}")
    if tuple(payload["profiles"]) != ("talk-v0-foldx", "talk-v0-no-foldx"):
        raise ValueError("profiles must contain talk-v0-foldx then talk-v0-no-foldx")
    if "exploratory_runs" in payload:
        raise ValueError("exploratory_runs are not part of the registered public workflow")
    ladder = _require_mapping(payload, "diagnostic_ladder")
    registered = tuple(ladder)
    if registered != APPROVED_LADDER_RUN_IDS:
        raise ValueError(
            f"registered run IDs must be {APPROVED_LADDER_RUN_IDS!r}; found {registered!r}"
        )
    output_policy = _require_mapping(payload, "output_policy")
    artifacts = output_policy.get("required_artifacts")
    if not isinstance(artifacts, list) or not artifacts or artifacts[-1] != "SUCCESS":
        raise ValueError("output_policy.required_artifacts must end with SUCCESS")
    required_artifacts = tuple(str(item) for item in artifacts)
    if len(set(required_artifacts)) != len(required_artifacts):
        raise ValueError("output_policy.required_artifacts must be unique")
    if output_policy.get("success_marker_last") is not True:
        raise ValueError("output_policy.success_marker_last must be true")
    kernel = _load_kernel(payload, seed)
    config = TalkV0Config(
        raw=payload,
        seed=seed,
        contact_cutoff_angstrom=cutoff,
        kernel=kernel,
        registered_run_ids=registered,
        required_artifacts=required_artifacts,
    )
    for run_id in APPROVED_LADDER_RUN_IDS:
        run_feature_columns(config, run_id)
    return config


def profile_feature_columns(config: TalkV0Config, profile_name: str) -> tuple[str, ...]:
    try:
        columns = config.profiles[profile_name]
    except KeyError as exc:
        known = ", ".join(sorted(config.profiles))
        raise KeyError(f"Unknown profile {profile_name!r}; known: {known}") from exc
    if not isinstance(columns, list) or not columns:
        raise ValueError(f"Profile {profile_name!r} must be a non-empty list")
    if len(set(columns)) != len(columns):
        raise ValueError(f"Profile {profile_name!r} has duplicate feature columns")
    return tuple(str(column) for column in columns)


def _resolve_named_columns(
    mapping: object, run_id: str, *, section: str
) -> tuple[str, ...] | None:
    if not isinstance(mapping, dict) or run_id not in mapping:
        return None
    columns = mapping[run_id]
    if not isinstance(columns, list) or not columns:
        raise ValueError(f"run_id {run_id!r} in {section} must be a non-empty list")
    if len(set(columns)) != len(columns):
        raise ValueError(f"run_id {run_id!r} in {section} has duplicate feature columns")
    return tuple(str(column) for column in columns)


def ladder_feature_columns(config: TalkV0Config, run_id: str) -> tuple[str, ...]:
    return run_feature_columns(config, run_id)


def run_feature_columns(config: TalkV0Config, run_id: str) -> tuple[str, ...]:
    resolved = _resolve_named_columns(
        config.raw.get("diagnostic_ladder"), run_id, section="diagnostic_ladder"
    )
    if resolved is None:
        known = ", ".join(APPROVED_LADDER_RUN_IDS)
        raise KeyError(f"Unknown run_id {run_id!r}; known: {known}")
    if run_id == "ladder-3-nofoldx":
        fallback = profile_feature_columns(config, "talk-v0-no-foldx")
        if resolved != fallback:
            raise ValueError(
                "ladder-3-nofoldx must match talk-v0-no-foldx column order"
            )
    return resolved
