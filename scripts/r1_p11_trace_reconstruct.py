#!/usr/bin/env python3
"""Reconstruct the complete p=11 quartic trace theta series exactly.

The profile counter supplies exact ordinary counts and common-coordinate
fourth moments.  Coordinate transitivity converts those moments into the
trace of the positive raw shell operator.  The corresponding harmonic trace
lies in the common 32-dimensional affine modular-form space of the three
square-circle channels.  A short exact prefix therefore determines every
exported coefficient and the half-cusp target.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from r1_p11_scalar_coupled_lp import CHANNELS, component_cases, load_rows


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_fraction_vector(values: object, name: str) -> list[Fraction]:
    if not isinstance(values, list):
        raise ValueError(f"{name} is not a list")
    return [Fraction(value) for value in values]


def load_scalar_reconstruction(path: Path) -> list[int]:
    report = json.loads(path.read_text())
    if report.get("status") != "complete_exact_modular_reconstruction":
        raise ValueError("scalar report is not a complete exact reconstruction")
    if report.get("remaining_dimension_after_half_target") != 0:
        raise ArithmeticError("scalar reconstruction retains a free coordinate")
    if report.get("all_reconstructed_coefficients_nonnegative_integers") is not True:
        raise ArithmeticError("scalar reconstruction lacks integral positivity")
    limit = int(report["coefficients_reconstructed_through"])
    counts = [0] * (limit + 1)
    nonzero = report.get("nonzero_coefficients")
    if not isinstance(nonzero, dict):
        raise ValueError("scalar reconstruction lacks nonzero_coefficients")
    for exponent_text, value in nonzero.items():
        exponent = int(exponent_text)
        if not 0 <= exponent <= limit or not isinstance(value, int) or value <= 0:
            raise ValueError(f"invalid scalar coefficient {exponent_text}: {value}")
        counts[exponent] = value
    return counts


def solve_affine_prefix(
    matrix: list[tuple[Fraction, ...]],
    rhs: list[Fraction],
    fixed_through: int,
) -> tuple[list[Fraction], list[int], int | None]:
    """Incremental exact RREF, retaining the exponent of every new pivot."""
    if not matrix or len(matrix) != len(rhs):
        raise ValueError("empty or inconsistent affine system")
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("inconsistent affine matrix width")
    basis: dict[int, tuple[list[Fraction], Fraction, int]] = {}
    pivot_exponents: list[int] = []
    first_full_rank: int | None = None
    for exponent in range(min(fixed_through, len(matrix) - 1) + 1):
        row = list(matrix[exponent])
        value = rhs[exponent]
        for pivot in sorted(basis):
            factor = row[pivot]
            if factor:
                pivot_row, pivot_value, _pivot_exponent = basis[pivot]
                row = [left - factor * right for left, right in zip(row, pivot_row)]
                value -= factor * pivot_value
        pivot = next((index for index, entry in enumerate(row) if entry), None)
        if pivot is None:
            if value:
                raise ArithmeticError(
                    f"trace prefix is inconsistent with the modular space at exponent {exponent}"
                )
            continue
        scale = row[pivot]
        row = [entry / scale for entry in row]
        value /= scale
        for old_pivot, (old_row, old_value, old_exponent) in list(basis.items()):
            factor = old_row[pivot]
            if factor:
                basis[old_pivot] = (
                    [left - factor * right for left, right in zip(old_row, row)],
                    old_value - factor * value,
                    old_exponent,
                )
        basis[pivot] = (row, value, exponent)
        pivot_exponents.append(exponent)
        if len(basis) == width and first_full_rank is None:
            first_full_rank = exponent

    if len(basis) != width:
        raise ArithmeticError(
            f"trace prefix has rank {len(basis)} in a {width}-dimensional affine space"
        )
    coordinates = [Fraction()] * width
    for pivot, (_row, value, _exponent) in basis.items():
        coordinates[pivot] = value
    return coordinates, pivot_exponents, first_full_rank


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile-report", type=Path, required=True)
    parser.add_argument("--scalar-reconstruction", type=Path, required=True)
    parser.add_argument("--affine-prefix", type=Path, required=True)
    parser.add_argument(
        "--affine-suffix",
        default="_qrows_exact_e800_20260828.txt",
    )
    parser.add_argument("--fixed-through", type=int, default=92)
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.p != 11:
        parser.error("the profile and affine exports are specific to p=11")

    profile = json.loads(args.profile_report.read_text())
    if profile.get("status") != "complete_exact_theta_prefix":
        raise ValueError("profile report is not marked complete")
    if profile.get("p") != args.p:
        raise ValueError("profile report has the wrong prime")
    required_checks = (
        "crt_product_exceeds_every_bound",
        "known_coefficients_match",
        "known_raw_trace_coefficients_match",
    )
    if any(profile.get(check) is not True for check in required_checks):
        raise ArithmeticError("profile report failed an exact calibration")
    profile_harmonic = parse_fraction_vector(
        profile.get("harmonic_trace_coefficients"),
        "harmonic_trace_coefficients",
    )
    profile_raw = parse_fraction_vector(
        profile.get("raw_quartic_trace_coefficients"),
        "raw_quartic_trace_coefficients",
    )
    if len(profile_harmonic) != int(profile["max_exponent"]) + 1:
        raise ValueError("profile harmonic trace length is inconsistent")
    if not 0 <= args.fixed_through < len(profile_harmonic):
        parser.error("--fixed-through lies outside the exact profile prefix")

    paths = {
        channel: Path(f"{args.affine_prefix}{channel}{args.affine_suffix}")
        for channel in CHANNELS
    }
    rows = {channel: load_rows(path) for channel, path in paths.items()}
    lengths = {len(channel_rows) for channel_rows in rows.values()}
    widths = {
        len(row)
        for channel_rows in rows.values()
        for row in channel_rows
    }
    if len(lengths) != 1 or len(widths) != 1:
        raise ValueError("harmonic channel exports have inconsistent shapes")
    row_count = lengths.pop()
    width = widths.pop()
    if row_count < 2 or width < 2:
        raise ValueError("harmonic channel exports are too small")
    reference = rows[CHANNELS[0]]
    if any(
        [row[1:] for row in rows[channel]] != [row[1:] for row in reference]
        for channel in CHANNELS[1:]
    ):
        raise ArithmeticError("channel homogeneous q-row matrices differ")

    channel_dimensions = {channel: 0 for channel in CHANNELS}
    for case in component_cases(args.p):
        channel_dimensions[str(case["channel"])] += int(
            case["component_count"]
        ) * int(case["representation_dimension"])
    n = args.p * args.p + 1
    d = n // 2
    zdim = n * (n - 6) // 8
    if sum(channel_dimensions.values()) != zdim:
        raise ArithmeticError("channel dimensions do not sum to dim Z")
    weighted_base = [
        sum(
            (
                Fraction(channel_dimensions[channel]) * rows[channel][index][0]
                for channel in CHANNELS
            ),
            Fraction(),
        )
        for index in range(row_count)
    ]
    matrix = [tuple(row[1:]) for row in reference]
    coefficient_base = weighted_base[:-1]
    coefficient_matrix = matrix[:-1]
    target_base = weighted_base[-1]
    target_row = matrix[-1]
    if len(coefficient_base) <= args.fixed_through:
        raise ValueError("affine coefficient export is shorter than the fixed prefix")

    rhs = [
        profile_harmonic[exponent] - coefficient_base[exponent]
        for exponent in range(len(profile_harmonic))
    ]
    coordinates, pivot_exponents, first_full_rank = solve_affine_prefix(
        coefficient_matrix,
        rhs + [Fraction()] * (len(coefficient_matrix) - len(rhs)),
        args.fixed_through,
    )
    reconstructed_harmonic = [
        base
        + sum(
            (coefficient * coordinate for coefficient, coordinate in zip(row, coordinates)),
            Fraction(),
        )
        for base, row in zip(coefficient_base, coefficient_matrix)
    ]
    if reconstructed_harmonic[: len(profile_harmonic)] != profile_harmonic:
        mismatch = next(
            exponent
            for exponent, (actual, expected) in enumerate(
                zip(reconstructed_harmonic, profile_harmonic)
            )
            if actual != expected
        )
        raise ArithmeticError(
            f"out-of-sample harmonic trace mismatch at {mismatch}: "
            f"{reconstructed_harmonic[mismatch]} != {profile_harmonic[mismatch]}"
        )

    counts = load_scalar_reconstruction(args.scalar_reconstruction)
    limit = min(len(counts), len(reconstructed_harmonic))
    counts = counts[:limit]
    reconstructed_harmonic = reconstructed_harmonic[:limit]
    reconstructed_raw: list[Fraction] = []
    for exponent, (count, harmonic) in enumerate(zip(counts, reconstructed_harmonic)):
        radius_squared = Fraction(exponent, 2 * args.p)
        radial = Fraction(2 * count, d * (d + 2)) * radius_squared**2
        raw = harmonic + zdim * radial
        if raw < 0:
            raise ArithmeticError(f"negative reconstructed raw trace at exponent {exponent}: {raw}")
        reconstructed_raw.append(raw)
    if reconstructed_raw[: len(profile_raw)] != profile_raw:
        mismatch = next(
            exponent
            for exponent, (actual, expected) in enumerate(
                zip(reconstructed_raw, profile_raw)
            )
            if actual != expected
        )
        raise ArithmeticError(
            f"out-of-sample raw trace mismatch at {mismatch}: "
            f"{reconstructed_raw[mismatch]} != {profile_raw[mismatch]}"
        )

    trace_target = target_base + sum(
        (coefficient * coordinate for coefficient, coordinate in zip(target_row, coordinates)),
        Fraction(),
    )
    report = {
        "experiment": "r1_p11_trace_reconstruct",
        "status": "complete_exact_trace_modular_reconstruction",
        "p": args.p,
        "channel_dimensions": channel_dimensions,
        "trace_affine_dimension": len(coordinates),
        "profile_fixed_through": args.fixed_through,
        "prefix_affine_rank": len(pivot_exponents),
        "first_full_rank_exponent": first_full_rank,
        "pivot_exponents": pivot_exponents,
        "profile_trace_available_through": len(profile_harmonic) - 1,
        "out_of_sample_profile_trace_coefficients_checked": (
            len(profile_harmonic) - args.fixed_through - 1
        ),
        "out_of_sample_profile_trace_coefficients_match": True,
        "common_second_moment_tight_frame_checks": len(
            profile.get(
                "second_moment_tight_frame_identity_checked_at_nonempty_exponents",
                [],
            )
        ),
        "known_raw_trace_coefficients_match": True,
        "coefficients_reconstructed_through": limit - 1,
        "all_reconstructed_raw_trace_coefficients_nonnegative": True,
        "trace_half_cusp_target": str(trace_target),
        "ordinary_theta_coefficients": counts,
        "harmonic_trace_coefficients": [str(value) for value in reconstructed_harmonic],
        "raw_trace_coefficients": [str(value) for value in reconstructed_raw],
        "inputs": {
            "profile_report": str(args.profile_report),
            "profile_report_sha256": sha256(args.profile_report),
            "scalar_reconstruction": str(args.scalar_reconstruction),
            "scalar_reconstruction_sha256": sha256(args.scalar_reconstruction),
            "harmonic_qrows": {
                channel: {"path": str(path), "sha256": sha256(path)}
                for channel, path in paths.items()
            },
        },
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items() if not key.endswith("coefficients")}, indent=2))


if __name__ == "__main__":
    main()
