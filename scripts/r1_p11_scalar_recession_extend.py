#!/usr/bin/env python3
"""Locate later exact shell constraints that cut a certified R1 LP ray."""
from __future__ import annotations

import argparse
import json
import re
from fractions import Fraction
from pathlib import Path

from r1_p11_scalar_coupled_exact_lp import build_scalar_trace_budget_model
from r1_p11_scalar_coupled_lp import load_rows


CHANNELS = ("circle-kernel", "circle-low", "circle-high")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ray-report", type=Path, required=True)
    parser.add_argument("--scalar-qrows", type=Path, required=True)
    parser.add_argument("--affine-prefix", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--sense", choices=("minimum", "maximum"), required=True)
    parser.add_argument("--parity-fourth-moment", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    source = json.loads(args.ray_report.read_text())
    case_rows = [row for row in source["cases"] if row["case"] == args.case]
    if len(case_rows) != 1:
        raise ValueError(f"case is missing or duplicated in {args.ray_report}")
    bound = case_rows[0]["bounds"][args.sense]
    if bound["status"] not in {"exact_unbounded", "exact_unbounded_ray"}:
        raise ValueError("source bound is not an exact unbounded certificate")
    ray = tuple(
        Fraction(value) for value in bound["unbounded_certificate"]["ray"]
    )

    harmonic = {
        channel: load_rows(
            Path(f"{args.affine_prefix}{channel}_qrows_exact_20260827.txt")
        )
        for channel in CHANNELS
    }
    model, _representatives = build_scalar_trace_budget_model(
        load_rows(args.scalar_qrows),
        harmonic,
        11,
        args.case,
        parity_fourth_moment=args.parity_fourth_moment,
    )
    if len(ray) != len(model.target):
        raise ValueError(f"ray width {len(ray)} != extended model width {len(model.target)}")

    violations: list[dict[str, str | int]] = []
    smallest: Fraction | None = None
    for constraint in model.constraints:
        coefficients = (
            constraint.coefficients
            if constraint.sense == ">="
            else tuple(-value for value in constraint.coefficients)
        )
        derivative = sum(
            (coefficient * value for coefficient, value in zip(coefficients, ray)),
            Fraction(),
        )
        if smallest is None or derivative < smallest:
            smallest = derivative
        if derivative < 0:
            exponent_match = re.search(r"(\d+)", constraint.name)
            violations.append(
                {
                    "constraint": constraint.name,
                    "exponent": int(exponent_match.group(1)) if exponent_match else -1,
                    "derivative": str(derivative),
                }
            )
    violations.sort(key=lambda row: (int(row["exponent"]), str(row["constraint"])))
    output = {
        "experiment": "r1_p11_scalar_recession_extend",
        "case": args.case,
        "sense": args.sense,
        "parity_fourth_moment": args.parity_fourth_moment,
        "extended_constraint_count": len(model.constraints),
        "violated_constraint_count": len(violations),
        "first_violated_exponent": violations[0]["exponent"] if violations else None,
        "first_violations": violations[:20],
        "minimum_recession_derivative": str(smallest),
    }
    payload = json.dumps(output, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
