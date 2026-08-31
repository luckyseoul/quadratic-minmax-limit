#!/usr/bin/env python3
"""Test an exact scalar-coupled recession ray against cusp-1/2 gaps."""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from r1_p11_scalar_coupled_exact_lp import scalar_affine_reduction
from r1_p11_scalar_coupled_lp import load_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ray-report", type=Path, required=True)
    parser.add_argument("--scalar-qrows", type=Path, required=True)
    parser.add_argument("--gap-rows", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--sense", choices=("minimum", "maximum"), required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    report = json.loads(args.ray_report.read_text())
    cases = [row for row in report["cases"] if row["case"] == args.case]
    if len(cases) != 1:
        raise ValueError("ray case is missing or duplicated")
    bound = cases[0]["bounds"][args.sense]
    ray = tuple(
        Fraction(value) for value in bound["unbounded_certificate"]["ray"]
    )

    scalar_rows = load_rows(args.scalar_qrows)
    _base, scalar_matrix, pivot_columns = scalar_affine_reduction(
        scalar_rows, 11, 28
    )
    scalar_dimension = len(scalar_matrix[0])
    reduced_direction = ray[-scalar_dimension:]
    if len(reduced_direction) != scalar_dimension:
        raise ValueError("ray is shorter than the scalar affine dimension")
    free_columns = tuple(
        index for index in range(len(scalar_rows[0])) if index not in pivot_columns
    )
    direction = [Fraction()] * len(scalar_rows[0])
    for column, value in zip(free_columns, reduced_direction):
        direction[column] = value

    gaps = load_rows(args.gap_rows)
    if any(len(row) != len(direction) for row in gaps):
        raise ValueError("gap and scalar basis widths differ")
    derivatives = [
        sum((coefficient * value for coefficient, value in zip(row, direction)), Fraction())
        for row in gaps
    ]
    nonzero = [index for index, value in enumerate(derivatives) if value]
    output = {
        "experiment": "r1_p11_scalar_half_gap_ray",
        "case": args.case,
        "sense": args.sense,
        "gap_rows": len(gaps),
        "gap_rank": 15,
        "nonzero_gap_rows": len(nonzero),
        "nonzero_gap_coefficients": sorted({index // 10 for index in nonzero}),
        "first_nonzero_row": nonzero[0] if nonzero else None,
        "first_nonzero_derivative": str(derivatives[nonzero[0]]) if nonzero else None,
        "ray_survives_half_gap": not nonzero,
    }
    payload = json.dumps(output, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
