#!/usr/bin/env python3
"""Reconstruct the p=11 ordinary theta form from its exact short prefix.

The profile/CRT engine supplies coefficients through exponent 87.  In the
W4-gap-reduced modular space those data leave one affine coordinate.  The
independently normalized cusp-1/2 coefficient fixes that coordinate, after
which every exported infinity coefficient is an exact nonnegative integer.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from r1_p11_scalar_coupled_exact_lp import scalar_affine_reduction
from r1_p11_scalar_coupled_lp import load_rows
from r1_p11_theta_endpoint_qsopt import load_exact_theta_prefix


PHASE_FACTORS = (
    Fraction(1),
    Fraction(),
    Fraction(),
    Fraction(),
    Fraction(1),
    Fraction(1),
    Fraction(1),
    Fraction(),
    Fraction(1),
    Fraction(1, 2),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--theta-report", type=Path, required=True)
    parser.add_argument("--scalar-qrows", type=Path, required=True)
    parser.add_argument("--scalar-half-target-rows", type=Path, required=True)
    parser.add_argument("--scalar-half-target-first", type=Fraction, required=True)
    parser.add_argument(
        "--fixed-through",
        type=int,
        help="use only profile coefficients through this exponent for reconstruction",
    )
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.p != 11:
        parser.error("the current exact exports are specific to p=11")

    available_prefix, available_fixed_counts = load_exact_theta_prefix(
        args.theta_report, args.p
    )
    fixed_through = (
        len(available_prefix) - 1
        if args.fixed_through is None
        else args.fixed_through
    )
    if not 28 <= fixed_through < len(available_prefix):
        parser.error(
            f"--fixed-through must lie in 28..{len(available_prefix) - 1}"
        )
    prefix = available_prefix[: fixed_through + 1]
    fixed_counts = {
        exponent: value
        for exponent, value in available_fixed_counts.items()
        if exponent <= fixed_through
    }
    coefficient_rows = load_rows(args.scalar_qrows)
    target_rows = load_rows(args.scalar_half_target_rows)
    if len(target_rows) != len(PHASE_FACTORS):
        raise ArithmeticError("unexpected scalar half-target row count")
    all_rows = coefficient_rows + target_rows
    bases, matrices, pivots = scalar_affine_reduction(
        all_rows,
        args.p,
        2 * (args.p + 3),
        fixed_counts,
    )
    coefficient_bases = bases[: len(coefficient_rows)]
    coefficient_matrices = matrices[: len(coefficient_rows)]
    target_bases = bases[len(coefficient_rows) :]
    target_matrices = matrices[len(coefficient_rows) :]
    free_dimension = len(coefficient_matrices[0])
    if free_dimension not in (0, 1):
        raise ArithmeticError(
            f"theta prefix leaves {free_dimension} scalar coordinates, expected at most one"
        )
    if any(
        base != factor * target_bases[0]
        or row != tuple(factor * value for value in target_matrices[0])
        for base, row, factor in zip(target_bases, target_matrices, PHASE_FACTORS)
    ):
        raise ArithmeticError("scalar half-target rows do not lie on the phase line")
    if free_dimension == 1:
        target_slope = target_matrices[0][0]
        if target_slope == 0:
            raise ArithmeticError("scalar half-target does not fix the remaining coordinate")
        coordinate: Fraction | None = (
            args.scalar_half_target_first - target_bases[0]
        ) / target_slope
        reconstructed = [
            base + row[0] * coordinate
            for base, row in zip(coefficient_bases, coefficient_matrices)
        ]
        half_target_rank = 1
        half_target_role = "fixes_last_affine_coordinate"
    else:
        coordinate = None
        if target_bases[0] != args.scalar_half_target_first:
            raise ArithmeticError(
                "prefix-determined modular form fails the independent half-target value"
            )
        reconstructed = coefficient_bases
        half_target_rank = 0
        half_target_role = "independent_prediction_matched"
    if any(value.denominator != 1 or value < 0 for value in reconstructed):
        bad = next(
            index
            for index, value in enumerate(reconstructed)
            if value.denominator != 1 or value < 0
        )
        raise ArithmeticError(
            f"reconstructed coefficient {bad} is not a nonnegative integer: "
            f"{reconstructed[bad]}"
        )
    integers = [int(value) for value in reconstructed]
    if integers[: len(prefix)] != prefix:
        mismatch = next(
            index
            for index, (actual, expected) in enumerate(zip(integers, prefix))
            if actual != expected
        )
        raise ArithmeticError(
            f"reconstructed theta prefix mismatch at {mismatch}: "
            f"{integers[mismatch]} != {prefix[mismatch]}"
        )
    available_match = integers[: len(available_prefix)] == available_prefix
    if not available_match:
        mismatch = next(
            index
            for index, (actual, expected) in enumerate(
                zip(integers, available_prefix)
            )
            if actual != expected
        )
        raise ArithmeticError(
            f"out-of-sample theta mismatch at {mismatch}: "
            f"{integers[mismatch]} != {available_prefix[mismatch]}"
        )

    nonzero = {str(index): value for index, value in enumerate(integers) if value}
    report = {
        "experiment": "r1_p11_scalar_theta_reconstruct",
        "status": "complete_exact_modular_reconstruction",
        "p": args.p,
        "profile_theta_available_through": len(available_prefix) - 1,
        "profile_theta_fixed_through": fixed_through,
        "scalar_space_dimension_before_prefix": len(coefficient_rows[0]),
        "prefix_affine_rank": len(coefficient_rows[0]) - free_dimension,
        "remaining_dimension_before_half_target": free_dimension,
        "theta_prefix_uniquely_determines_modular_form": free_dimension == 0,
        "half_target_rank": half_target_rank,
        "half_target_role": half_target_role,
        "half_target_value_matches": True,
        "remaining_dimension_after_half_target": 0,
        "pivot_columns": list(pivots),
        "remaining_coordinate": str(coordinate) if coordinate is not None else None,
        "coefficients_reconstructed_through": len(integers) - 1,
        "all_reconstructed_coefficients_nonnegative_integers": True,
        "profile_prefix_matches_reconstruction": True,
        "out_of_sample_profile_coefficients_checked": (
            len(available_prefix) - len(prefix)
        ),
        "out_of_sample_profile_coefficients_match": available_match,
        "nonzero_coefficients": nonzero,
        "inputs": {
            "theta_report": str(args.theta_report),
            "theta_report_sha256": sha256(args.theta_report),
            "scalar_qrows": str(args.scalar_qrows),
            "scalar_qrows_sha256": sha256(args.scalar_qrows),
            "scalar_half_target_rows": str(args.scalar_half_target_rows),
            "scalar_half_target_rows_sha256": sha256(args.scalar_half_target_rows),
            "scalar_half_target_first": str(args.scalar_half_target_first),
        },
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
