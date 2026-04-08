#!/usr/bin/env python3
"""Example script for loading and simulating Human-GEM with cobrapy."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import cobra
from cobra.util.solver import linear_reaction_coefficients

from src.utils.path import get_project_root

REPO_ROOT = get_project_root()
DEFAULT_MODEL_PATH = REPO_ROOT / "model" / "Human-GEM.yml"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Load Human-GEM.yml, report a short model summary, run FBA, "
            "and show the highest-flux reactions."
        )
    )
    parser.add_argument(
        "--model",
        type=Path,
        default=DEFAULT_MODEL_PATH,
        help=f"Path to the Human-GEM YAML model (default: {DEFAULT_MODEL_PATH})",
    )
    parser.add_argument(
        "--objective",
        help="Optional reaction ID to use as the objective instead of the model default.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Number of top absolute-flux reactions to print (default: 10).",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        help="Optional path to save the top-flux table as CSV.",
    )
    return parser.parse_args()


def objective_to_string(model: cobra.Model) -> str:
    coefficients = linear_reaction_coefficients(model)
    if coefficients:
        return ", ".join(
            f"{reaction.id} ({coefficient:g})"
            for reaction, coefficient in coefficients.items()
        )
    return str(model.objective.expression)


def build_top_flux_rows(
    model: cobra.Model,
    solution: Any,
    top_n: int,
) -> list[dict[str, Any]]:
    non_zero_fluxes = solution.fluxes[solution.fluxes.abs() > 1e-9]
    sorted_fluxes = non_zero_fluxes.reindex(
        non_zero_fluxes.abs().sort_values(ascending=False).index
    ).head(top_n)

    rows: list[dict[str, Any]] = []
    for reaction_id, flux in sorted_fluxes.items():
        reaction = model.reactions.get_by_id(reaction_id)
        rows.append(
            {
                "reaction_id": reaction.id,
                "reaction_name": reaction.name,
                "flux": float(flux),
                "lower_bound": float(reaction.lower_bound),
                "upper_bound": float(reaction.upper_bound),
            }
        )
    return rows


def write_csv(rows: list[dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "reaction_id",
                "reaction_name",
                "flux",
                "lower_bound",
                "upper_bound",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()

    if args.top < 1:
        raise SystemExit("--top must be at least 1.")

    model_path = args.model.resolve()
    if not model_path.exists():
        raise SystemExit(f"Model file not found: {model_path}")

    model = cobra.io.load_yaml_model(str(model_path))

    if args.objective:
        try:
            model.objective = model.reactions.get_by_id(args.objective)
        except KeyError as exc:
            raise SystemExit(f"Objective reaction not found: {args.objective}") from exc

    print(f"Model file: {model_path}")
    print(f"Solver: {model.solver.interface.__name__}")
    print(
        "Model size: "
        f"{len(model.reactions)} reactions, "
        f"{len(model.metabolites)} metabolites, "
        f"{len(model.genes)} genes, "
        f"{len(model.exchanges)} exchange reactions"
    )
    print(f"Objective: {objective_to_string(model)}")
    print(f"Objective direction: {model.objective_direction}")

    solution = model.optimize()
    print(f"Optimization status: {solution.status}")
    print(f"Objective value: {float(solution.objective_value):.6f}")

    top_flux_rows = build_top_flux_rows(model, solution, args.top)
    print(f"Top {len(top_flux_rows)} reactions by absolute flux:")
    for row in top_flux_rows:
        print(
            f"  {row['reaction_id']:>8}  {row['flux']:>12.6f}  {row['reaction_name']}"
        )

    if args.output_csv:
        output_path = args.output_csv.resolve()
        write_csv(top_flux_rows, output_path)
        print(f"Saved CSV: {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
