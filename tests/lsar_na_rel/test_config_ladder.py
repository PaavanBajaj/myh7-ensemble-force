"""Ladder freeze: run_id → columns before any target is loaded."""

import json
from pathlib import Path

import pytest

from lsar_na_rel.config import (
    APPROVED_LADDER_RUN_IDS,
    load_talk_v0_config,
    ladder_feature_columns,
    profile_feature_columns,
    run_feature_columns,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = REPO_ROOT / "workflows" / "lsar-na-rel" / "config" / "lsar-na-rel.json"


@pytest.fixture(scope="module")
def config():
    return load_talk_v0_config(CONFIG_PATH)


def test_approved_ladder_run_ids_are_exactly_three_in_order():
    assert APPROVED_LADDER_RUN_IDS == (
        "ladder-1-rsa",
        "ladder-2-rsa-ihm",
        "ladder-3-nofoldx",
    )


def test_ladder_column_order_is_frozen(config):
    assert ladder_feature_columns(config, "ladder-1-rsa") == (
        "wt_site_rsa_8act_mean_ab",
    )
    assert ladder_feature_columns(config, "ladder-2-rsa-ihm") == (
        "wt_site_rsa_8act_mean_ab",
        "ihm_interface_flag_8act",
    )
    assert ladder_feature_columns(config, "ladder-3-nofoldx") == (
        "delta_charge",
        "wt_site_rsa_8act_mean_ab",
        "ihm_interface_flag_8act",
    )


def test_ladder_three_matches_talk_v0_no_foldx_profile(config):
    assert ladder_feature_columns(config, "ladder-3-nofoldx") == profile_feature_columns(
        config, "talk-v0-no-foldx"
    )


def test_unknown_run_id_raises(config):
    with pytest.raises(KeyError, match="ladder-4"):
        ladder_feature_columns(config, "ladder-4-extra")


def test_public_config_excludes_exploratory_probe_runs(config):
    assert "exploratory_runs" not in config.raw
    assert tuple(config.profiles) == ("talk-v0-foldx", "talk-v0-no-foldx")
    with pytest.raises(KeyError, match="probe-ihm-min-dist"):
        run_feature_columns(config, "probe-ihm-min-dist")


def test_runtime_settings_are_loaded_from_config(config):
    assert config.seed == 20260811
    assert config.contact_cutoff_angstrom == 4.5
    assert config.kernel.type == "ConstantKernel * Matern + WhiteKernel"
    assert config.kernel.nu == 1.5
    assert config.kernel.constant_value == 1.0
    assert config.kernel.length_scale == 1.0
    assert config.kernel.noise_level == 1.0
    assert config.kernel.normalize_y is True
    assert config.kernel.n_restarts_optimizer == 10
    assert config.kernel.random_state == config.seed
    assert config.registered_run_ids == APPROVED_LADDER_RUN_IDS
    assert config.required_artifacts[-1] == "SUCCESS"


def _write_changed_config(tmp_path: Path, change) -> Path:
    payload = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    change(payload)
    path = tmp_path / "config.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_rejects_mismatched_kernel_seed(tmp_path: Path):
    path = _write_changed_config(
        tmp_path, lambda payload: payload["kernel"].__setitem__("random_state", 7)
    )
    with pytest.raises(ValueError, match="random_state.*seed"):
        load_talk_v0_config(path)


def test_rejects_unsupported_kernel_type(tmp_path: Path):
    path = _write_changed_config(
        tmp_path, lambda payload: payload["kernel"].__setitem__("type", "RBF")
    )
    with pytest.raises(ValueError, match="Unsupported kernel"):
        load_talk_v0_config(path)


def test_rejects_nonpositive_contact_cutoff(tmp_path: Path):
    path = _write_changed_config(
        tmp_path, lambda payload: payload.__setitem__("contact_cutoff_angstrom", 0)
    )
    with pytest.raises(ValueError, match="contact_cutoff_angstrom"):
        load_talk_v0_config(path)


def test_rejects_registered_ladder_drift(tmp_path: Path):
    path = _write_changed_config(
        tmp_path,
        lambda payload: payload["diagnostic_ladder"].__setitem__("ladder-4", ["x"]),
    )
    with pytest.raises(ValueError, match="registered run IDs"):
        load_talk_v0_config(path)
