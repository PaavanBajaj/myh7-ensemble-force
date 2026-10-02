"""Command-line entry point for the versioned MYH7 exploratory screen."""

from __future__ import annotations

import argparse
from pathlib import Path



def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    feature_parser = commands.add_parser("features", help="generate five label-blind structural features")
    for name in ("variants", "fasta", "assembly-8act", "cif-8efe", "cif-8efd", "cif-8efi", "output-dir"):
        feature_parser.add_argument(f"--{name}", required=True, type=Path)
    outcome_parser = commands.add_parser("outcomes", help="assemble canonical outcomes and source links")
    for name in ("kcat-labels", "kcat-evidence", "na-primary", "na-labels", "na-measurements", "velocity-summary", "morck-pairs", "output-dir"):
        outcome_parser.add_argument(f"--{name}", required=True, type=Path)
    plan_parser = commands.add_parser("plan", help="freeze pair inventory without effect estimates")
    for name in ("outcomes", "new-features", "kcat-features", "na-features", "config", "output-dir"):
        plan_parser.add_argument(f"--{name}", required=True, type=Path)
    analysis_parser = commands.add_parser("analyze", help="compute the frozen descriptive screen and sensitivities")
    for name in ("table", "inventory", "plan-manifest", "outcomes", "dependence-ledger", "velocity-alternatives", "output-dir"):
        analysis_parser.add_argument(f"--{name}", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "features":
        from .features import generate_features

        generate_features(
            variants_path=args.variants, fasta_path=args.fasta,
            assembly_path=args.assembly_8act, adp_path=args.cif_8efe,
            rigor_path=args.cif_8efd, actin_path=args.cif_8efi,
            output_dir=args.output_dir,
        )
    elif args.command == "outcomes":
        from .outcomes import assemble_outcomes

        assemble_outcomes(
            kcat_labels=args.kcat_labels, kcat_evidence=args.kcat_evidence,
            na_primary=args.na_primary, na_labels=args.na_labels,
            na_measurements=args.na_measurements,
            velocity_summary=args.velocity_summary, morck_pairs=args.morck_pairs,
            output_dir=args.output_dir,
        )
    elif args.command == "plan":
        from .planning import freeze_plan

        freeze_plan(
            outcomes_path=args.outcomes, new_features_path=args.new_features,
            kcat_features_path=args.kcat_features, na_features_path=args.na_features,
            config_path=args.config, output_dir=args.output_dir,
        )
    elif args.command == "analyze":
        from .analysis import analyze

        analyze(
            table_path=args.table, inventory_path=args.inventory,
            plan_manifest_path=args.plan_manifest, outcomes_path=args.outcomes,
            dependence_ledger_path=args.dependence_ledger,
            velocity_alternatives_path=args.velocity_alternatives,
            output_dir=args.output_dir,
        )


if __name__ == "__main__":
    main()
