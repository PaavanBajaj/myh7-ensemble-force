"""Talk-v0 CLI: one run_id or all three ladder runs."""

from __future__ import annotations

import argparse
import os
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

from lsar_na_rel.config import APPROVED_LADDER_RUN_IDS, load_talk_v0_config
from lsar_na_rel.manifests import (
    FOLDX_FALLBACK_REASON,
    build_run_manifest,
    write_run_manifest,
)
from lsar_na_rel.matrix import build_talk_v0_matrix
from lsar_na_rel.paths import DEFAULT_PATHS
from lsar_na_rel.reporting import build_metrics_payload, write_reports
from lsar_na_rel.validation import run_loso_gpr_and_baseline

REPO_ROOT = DEFAULT_PATHS.repo_root
DEFAULT_PRIMARY = DEFAULT_PATHS.primary_path
DEFAULT_LABELS = DEFAULT_PATHS.labels_path
DEFAULT_FEATURES = DEFAULT_PATHS.features_path
DEFAULT_CONFIG = DEFAULT_PATHS.config_path
DEFAULT_OUTPUT_ROOT = DEFAULT_PATHS.results_dir

_BLAS_THREAD_KEYS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def set_blas_threads_one() -> None:
    for key in _BLAS_THREAD_KEYS:
        os.environ[key] = "1"


def select_run_ids(run_id: str | None, all_ladder: bool) -> tuple[str, ...]:
    if bool(run_id) == bool(all_ladder):
        raise ValueError("exactly one of --run-id or --all-ladder is required")
    if all_ladder:
        return APPROVED_LADDER_RUN_IDS
    assert run_id is not None
    if run_id not in APPROVED_LADDER_RUN_IDS:
        raise ValueError(
            f"--run-id must be a registered run_id: {', '.join(APPROVED_LADDER_RUN_IDS)}"
        )
    return (run_id,)


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def make_run_dir(output_root: Path, run_id: str, timestamp: str) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    run_dir = output_root / f"{run_id}-{timestamp}"
    if run_dir.exists():
        micros = datetime.now(timezone.utc).strftime("%f")
        run_dir = output_root / f"{run_id}-{timestamp}-{micros}"
    run_dir.mkdir(parents=True, exist_ok=False)
    return run_dir


def write_success_marker(run_dir: Path, required_artifacts: tuple[str, ...]) -> None:
    for name in required_artifacts:
        if name == "SUCCESS":
            continue
        path = run_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"required artifact missing or empty: {path}")
    (run_dir / "SUCCESS").write_text("SUCCESS\n", encoding="utf-8")


def _profile_for_run(run_id: str) -> str:
    if run_id == "ladder-3-nofoldx":
        return "talk-v0-no-foldx"
    return "diagnostic-ladder"


def _append_log(run_log: Path, message: str) -> None:
    with run_log.open("a", encoding="utf-8") as handle:
        handle.write(message)
        if not message.endswith("\n"):
            handle.write("\n")


def _portable_path(path: Path, repo_root: Path = REPO_ROOT) -> str:
    try:
        return path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return path.name


def _print_console_summary(
    *,
    run_id: str,
    oof,
    matrix,
    run_dir: Path,
    seed: int,
) -> None:
    payload = build_metrics_payload(run_id, matrix.feature_columns, seed, oof)
    n_rows = int(payload["n_rows"])
    n_folds = int(payload["n_folds"])
    n_groups = int(matrix.groups.nunique())
    lines = [
        f"profile={_profile_for_run(run_id)} run_id={run_id}",
        f"fallback_reason={FOLDX_FALLBACK_REASON}",
    ]
    if n_rows == 16 and n_groups == 6 and n_folds == 6:
        lines.append("16 rows / six study groups / six outer folds")
    else:
        lines.append(f"{n_rows} rows / {n_groups} study groups / {n_folds} outer folds")
    lines.extend(
        [
            (
                f"GPR MAE={payload['gpr_mae']:.6g} RMSE={payload['gpr_rmse']:.6g}; "
                f"baseline MAE={payload['baseline_mae']:.6g} "
                f"RMSE={payload['baseline_rmse']:.6g}"
            ),
            (
                f"delta_mae_gpr_minus_baseline={payload['delta_mae_gpr_minus_baseline']:.6g}; "
                f"delta_rmse_gpr_minus_baseline="
                f"{payload['delta_rmse_gpr_minus_baseline']:.6g}"
            ),
            f"artifacts={run_dir}",
            "exploratory small-n result; Risky rows included",
        ]
    )
    print("\n".join(lines))


def run_one(
    *,
    run_id: str,
    primary_path: Path,
    labels_path: Path,
    features_path: Path,
    output_root: Path,
    config_path: Path,
    command: str,
) -> Path:
    set_blas_threads_one()
    timestamp = utc_timestamp()
    run_dir = make_run_dir(output_root, run_id, timestamp)
    run_log = run_dir / "run.log"
    run_log.write_text(f"start run_id={run_id} timestamp={timestamp}\n", encoding="utf-8")
    try:
        config = load_talk_v0_config(config_path)
        _append_log(run_log, f"loaded config={_portable_path(config_path)}")
        matrix = build_talk_v0_matrix(
            run_id=run_id,
            config=config,
            primary_path=primary_path,
            labels_path=labels_path,
            features_path=features_path,
        )
        _append_log(
            run_log,
            f"matrix rows={len(matrix.y)} features={list(matrix.feature_columns)}",
        )
        oof = run_loso_gpr_and_baseline(
            matrix.X,
            matrix.y,
            matrix.groups,
            kernel_settings=config.kernel,
            random_state=config.seed,
        )
        _append_log(
            run_log,
            f"validation splitter={oof.splitter_name} n_folds={oof.n_folds}",
        )
        write_reports(
            run_dir,
            run_id=run_id,
            feature_columns=matrix.feature_columns,
            seed=config.seed,
            variants=matrix.variants,
            source_ids=matrix.groups,
            oof=oof,
        )
        _append_log(run_log, "wrote reports")
        manifest = build_run_manifest(
            run_id=run_id,
            feature_columns=matrix.feature_columns,
            seed=config.seed,
            splitter_name=oof.splitter_name,
            primary_path=primary_path,
            labels_path=labels_path,
            features_path=features_path,
            repo_root=REPO_ROOT,
            command=command,
            timestamp=timestamp,
            n_rows=len(matrix.y),
            n_folds=oof.n_folds,
        )
        write_run_manifest(run_dir / "run-manifest.json", manifest)
        _append_log(run_log, "wrote run-manifest.json")
        write_success_marker(run_dir, config.required_artifacts)
        _print_console_summary(
            run_id=run_id,
            oof=oof,
            matrix=matrix,
            run_dir=run_dir,
            seed=config.seed,
        )
        return run_dir
    except Exception:
        if run_log.exists():
            _append_log(run_log, traceback.format_exc())
        raise


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="lsar_na_rel.cli")
    sub = parser.add_subparsers(dest="command", required=True)
    talk = sub.add_parser("talk-v0", help="Run Talk-v0 diagnostic ladder")
    mode = talk.add_mutually_exclusive_group(required=True)
    mode.add_argument("--run-id", type=str, default=None)
    mode.add_argument("--all-ladder", action="store_true")
    talk.add_argument("--primary", type=Path, default=DEFAULT_PRIMARY)
    talk.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    talk.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    talk.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    talk.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    return parser


def main(argv: list[str] | None = None) -> int:
    argv_list = list(sys.argv[1:] if argv is None else argv)
    command = "python -m lsar_na_rel.cli " + " ".join(argv_list)
    parser = _build_parser()
    args = parser.parse_args(argv_list)
    if args.command != "talk-v0":
        parser.error(f"unknown command {args.command!r}")
    try:
        run_ids = select_run_ids(args.run_id, args.all_ladder)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    failed = False
    for run_id in run_ids:
        try:
            run_one(
                run_id=run_id,
                primary_path=Path(args.primary),
                labels_path=Path(args.labels),
                features_path=Path(args.features),
                output_root=Path(args.output_root),
                config_path=Path(args.config),
                command=command,
            )
        except Exception as exc:
            print(f"run_id={run_id} failed: {exc}", file=sys.stderr)
            failed = True
            break
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
