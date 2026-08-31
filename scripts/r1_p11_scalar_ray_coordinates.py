#!/usr/bin/env python3
"""Recover full modular coordinates of a scalar LP recession direction."""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from r1_p11_scalar_coupled_exact_lp import scalar_affine_reduction
from r1_p11_scalar_coupled_lp import load_rows


def gp_vector(values: list[Fraction]) -> str:
    return "[" + ",".join(str(value) for value in values) + "]\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ray-report", type=Path, required=True)
    parser.add_argument("--scalar-qrows", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--sense", choices=("minimum", "maximum"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    report = json.loads(args.ray_report.read_text())
    cases = [row for row in report["cases"] if row["case"] == args.case]
    if len(cases) != 1:
        raise ValueError("ray case is missing or duplicated")
    case = cases[0]
    scalar_fixed_counts = {
        int(exponent): Fraction(count)
        for exponent, count in case.get("scalar_fixed_counts", {}).items()
    }
    rows = load_rows(args.scalar_qrows)
    _base, affine_rows, _pivots = scalar_affine_reduction(
        rows, 11, 28, scalar_fixed_counts
    )
    ray = [
        Fraction(value)
        for value in case["bounds"][args.sense]["unbounded_certificate"]["ray"]
    ]
    affine_dimension = len(affine_rows[0])
    scalar_ray = ray[-affine_dimension:]
    coefficient_direction = [
        sum((value * coordinate for value, coordinate in zip(row, scalar_ray)), Fraction())
        for row in affine_rows
    ]

    width = len(rows[0])
    identity_rows: list[int] = []
    for column in range(width):
        expected = [Fraction()] * width
        expected[column] = Fraction(1)
        matches = [index for index, row in enumerate(rows) if row == expected]
        if len(matches) != 1:
            raise ArithmeticError(f"scalar pivot column {column} has {len(matches)} identity rows")
        identity_rows.append(matches[0])
    coordinates = [coefficient_direction[index] for index in identity_rows]
    for exponent, row in enumerate(rows):
        reconstructed = sum(
            (value * coordinate for value, coordinate in zip(row, coordinates)),
            Fraction(),
        )
        if reconstructed != coefficient_direction[exponent]:
            raise ArithmeticError(f"modular-coordinate reconstruction fails at {exponent}")

    args.output.write_text(gp_vector(coordinates))
    print(
        json.dumps(
            {
                "experiment": "r1_p11_scalar_ray_coordinates",
                "sense": args.sense,
                "coordinate_count": len(coordinates),
                "nonzero_coordinates": sum(value != 0 for value in coordinates),
                "pivot_exponents": identity_rows,
                "scalar_fixed_counts": {
                    str(exponent): str(count)
                    for exponent, count in sorted(scalar_fixed_counts.items())
                },
                "first_nonzero_coefficient_exponents": [
                    exponent
                    for exponent, value in enumerate(coefficient_direction)
                    if value
                ][:20],
                "output": str(args.output),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
