"""Frozen descriptive nonlinear dependence screen. No actual-cohort inference."""
import argparse
from pathlib import Path

from .analysis import CONFIG, OUTPUT, PLAN, SOURCE, analyze, plan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["plan", "analyze"])
    parser.add_argument("--source-dir", type=Path, default=SOURCE)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--plan", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "plan":
            plan(args.source_dir, args.config, args.output_dir or PLAN.parent)
        else:
            analyze(args.source_dir, args.config, args.plan or (args.output_dir / "dependence-plan-v1.json" if args.output_dir else PLAN), args.output_dir or OUTPUT)
    except (ValueError, KeyError, OSError, TypeError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
