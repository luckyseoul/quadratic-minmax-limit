#!/usr/bin/env python3
"""Reconstruct the three p=11 broad-channel trace theta series exactly.

The marked profile counter supplies exact raw shell mass in the kernel, low,
and high eigenspaces of the square-circle operator through exponent 120.
Dividing by each broad dimension and removing the universal radial shift gives
a harmonic theta series in that channel's 32-dimensional affine modular
space.  Exact prefix elimination determines all exported coefficients and
the half-cusp target, with the remaining profile rows held out.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from r1_p11_scalar_coupled_lp import CHANNELS, load_rows
from r1_p11_trace_reconstruct import load_scalar_reconstruction, solve_affine_prefix


MASS_KEYS = {
    "circle-kernel": "kernel",
    "circle-low": "low",
    "circle-high": "high",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--channel-moments", type=Path, required=True)
    parser.add_argument("--scalar-reconstruction", type=Path, required=True)
    parser.add_argument("--trace-reconstruction", type=Path, required=True)
    parser.add_argument("--exact-row-directory", type=Path, required=True)
    parser.add_argument(
        "--row-suffix", default="_qrows_exact_e800_20260828.txt"
    )
    parser.add_argument("--fixed-through", type=int, default=92)
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.p != 11:
        parser.error("the exact profile and q-row data are specific to p=11")

    moments = json.loads(args.channel_moments.read_text())
    required = (
        moments.get("status") == "complete_exact_broad_channel_mass_prefix",
        moments.get("p") == args.p,
        moments.get("crt_product_exceeds_every_bound") is True,
        moments.get("ordinary_coefficients_match_prior_exact_profile_report") is True,
        moments.get("classified_shell_channel_masses_match") is True,
    )
    if not all(required):
        raise ArithmeticError("channel-moment report failed a required exact check")
    profile_rows = moments.get("rows")
    if not isinstance(profile_rows, list) or len(profile_rows) != int(
        moments["max_exponent"]
    ) + 1:
        raise ValueError("channel-moment rows have inconsistent length")
    if not 0 <= args.fixed_through < len(profile_rows):
        parser.error("--fixed-through lies outside the exact profile prefix")

    dimensions = {
        channel: int(moments["broad_channel_dimensions"][channel])
        for channel in CHANNELS
    }
    p = args.p
    n = p * p + 1
    d = n // 2
    zdim = n * (n - 6) // 8
    if sum(dimensions.values()) != zdim:
        raise ArithmeticError("broad dimensions do not sum to dim Z")

    counts = load_scalar_reconstruction(args.scalar_reconstruction)
    trace_report = json.loads(args.trace_reconstruction.read_text())
    if trace_report.get("status") != "complete_exact_trace_modular_reconstruction":
        raise ValueError("aggregate trace reconstruction is not exact/complete")
    aggregate_raw = [Fraction(value) for value in trace_report["raw_trace_coefficients"]]
    aggregate_target = Fraction(trace_report["trace_half_cusp_target"])

    channel_outputs: dict[str, object] = {}
    reconstructed_mass: dict[str, list[Fraction]] = {}
    target_weighted_sum = Fraction()
    common_pivots: list[int] | None = None
    first_full_rank: int | None = None
    qrow_inputs = {}

    for channel in CHANNELS:
        path = args.exact_row_directory / f"p11_affine_{channel}{args.row_suffix}"
        rows = load_rows(path)
        if len(rows) < 2:
            raise ValueError(f"q-row export is too short: {path}")
        width = len(rows[0])
        if width < 2 or any(len(row) != width for row in rows):
            raise ValueError(f"q-row export has inconsistent width: {path}")
        coefficient_rows = rows[:-1]
        target_row = rows[-1]
        base = [row[0] for row in coefficient_rows]
        matrix = [tuple(row[1:]) for row in coefficient_rows]
        target_base = target_row[0]
        target = tuple(target_row[1:])
        if len(base) <= args.fixed_through:
            raise ValueError(f"q-row export ends before fixed prefix: {path}")

        dimension = dimensions[channel]
        mass_key = MASS_KEYS[channel]
        profile_raw_average = [
            Fraction(row[mass_key]) / dimension for row in profile_rows
        ]
        profile_harmonic = []
        for exponent, (row, raw) in enumerate(zip(profile_rows, profile_raw_average)):
            shell_count = int(row["shell_count"])
            radius = Fraction(exponent, 2 * p)
            radial = Fraction(2 * shell_count, d * (d + 2)) * radius**2
            profile_harmonic.append(raw - radial)

        rhs = [
            profile_harmonic[exponent] - base[exponent]
            for exponent in range(len(profile_harmonic))
        ]
        coordinates, pivots, channel_full_rank = solve_affine_prefix(
            matrix,
            rhs + [Fraction()] * (len(matrix) - len(rhs)),
            args.fixed_through,
        )
        reconstructed_harmonic = [
            constant
            + sum(
                (coefficient * coordinate for coefficient, coordinate in zip(row, coordinates)),
                Fraction(),
            )
            for constant, row in zip(base, matrix)
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
                f"{channel} held-out profile mismatch at exponent {mismatch}"
            )

        limit = min(len(counts), len(reconstructed_harmonic))
        raw_average = []
        broad_mass = []
        for exponent in range(limit):
            radius = Fraction(exponent, 2 * p)
            radial = Fraction(2 * counts[exponent], d * (d + 2)) * radius**2
            raw = reconstructed_harmonic[exponent] + radial
            if raw < 0:
                raise ArithmeticError(
                    f"negative {channel} average raw mass at exponent {exponent}: {raw}"
                )
            raw_average.append(raw)
            broad_mass.append(dimension * raw)
        if broad_mass[: len(profile_rows)] != [
            Fraction(row[mass_key]) for row in profile_rows
        ]:
            raise ArithmeticError(f"{channel} reconstructed raw prefix changed")

        channel_target = target_base + sum(
            (coefficient * coordinate for coefficient, coordinate in zip(target, coordinates)),
            Fraction(),
        )
        target_weighted_sum += dimension * channel_target
        reconstructed_mass[channel] = broad_mass
        if common_pivots is None:
            common_pivots = pivots
            first_full_rank = channel_full_rank
        elif pivots != common_pivots or channel_full_rank != first_full_rank:
            raise ArithmeticError("broad channels have different prefix rank profiles")
        qrow_inputs[channel] = {"path": str(path), "sha256": sha256(path)}
        channel_outputs[channel] = {
            "dimension": dimension,
            "affine_dimension": len(coordinates),
            "first_full_rank_exponent": channel_full_rank,
            "pivot_exponents": pivots,
            "profile_fixed_through": args.fixed_through,
            "held_out_profile_coefficients_matched": (
                len(profile_harmonic) - args.fixed_through - 1
            ),
            "half_cusp_average_harmonic_target": str(channel_target),
            "harmonic_average_coefficients": [
                str(value) for value in reconstructed_harmonic[:limit]
            ],
            "raw_average_coefficients": [str(value) for value in raw_average],
            "broad_raw_mass_coefficients": [str(value) for value in broad_mass],
            "all_reconstructed_broad_masses_nonnegative": True,
        }

    limit = min(len(aggregate_raw), *(len(values) for values in reconstructed_mass.values()))
    for exponent in range(limit):
        total = sum(
            (reconstructed_mass[channel][exponent] for channel in CHANNELS),
            Fraction(),
        )
        if total != aggregate_raw[exponent]:
            raise ArithmeticError(
                f"broad masses fail aggregate trace at exponent {exponent}: "
                f"{total} != {aggregate_raw[exponent]}"
            )
    if target_weighted_sum != aggregate_target:
        raise ArithmeticError(
            "broad half-cusp targets fail aggregate trace: "
            f"{target_weighted_sum} != {aggregate_target}"
        )

    report = {
        "experiment": "r1_p11_broad_channel_reconstruct",
        "status": "complete_exact_broad_channel_modular_reconstruction",
        "p": p,
        "profile_available_through": len(profile_rows) - 1,
        "profile_fixed_through": args.fixed_through,
        "common_affine_dimension": len(common_pivots or []),
        "common_first_full_rank_exponent": first_full_rank,
        "common_pivot_exponents": common_pivots,
        "coefficients_reconstructed_through": limit - 1,
        "all_three_broad_mass_series_nonnegative": True,
        "broad_masses_sum_to_aggregate_trace_through_reconstruction": True,
        "broad_targets_sum_to_aggregate_trace_target": True,
        "aggregate_trace_target": str(aggregate_target),
        "channels": channel_outputs,
        "inputs": {
            "channel_moments": {
                "path": str(args.channel_moments),
                "sha256": sha256(args.channel_moments),
            },
            "scalar_reconstruction": {
                "path": str(args.scalar_reconstruction),
                "sha256": sha256(args.scalar_reconstruction),
            },
            "trace_reconstruction": {
                "path": str(args.trace_reconstruction),
                "sha256": sha256(args.trace_reconstruction),
            },
            "harmonic_qrows": qrow_inputs,
        },
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    compact = {key: value for key, value in report.items() if key != "channels"}
    compact["channel_targets"] = {
        channel: channel_outputs[channel]["half_cusp_average_harmonic_target"]
        for channel in CHANNELS
    }
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
