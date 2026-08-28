#!/usr/bin/env python3
"""Certify p=11 R1 endpoints with separate broad-channel conservation.

The marked profile census reconstructs the raw shell trace and transformed
target average on each eigenspace of the square-circle operator.  This script
imposes those three exact conservation laws on every shell, optimizes each
irreducible target with QSopt_ex over the rationals, verifies primal and dual
certificates independently, and converts the certified intervals into a
channel-mean-aware upper bound for the R1 spectral variance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from fractions import Fraction
from pathlib import Path

from r1_p11_scalar_coupled_exact_lp import condition_model
from r1_p11_scalar_coupled_lp import CHANNELS, component_cases
from r1_p11_trace_coupled_exact_lp import (
    build_broad_channel_conserved_model,
    load_qrows,
    maximum_weighted_variance_exact,
    run_qsopt,
    verify_certificate,
    write_lp,
)
from r1_p11_trace_reconstruct import load_scalar_reconstruction


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_broad_reconstruction(
    path: Path, p: int
) -> tuple[
    dict[str, int],
    dict[str, list[Fraction]],
    dict[str, Fraction],
    Fraction,
    dict[str, object],
]:
    report = json.loads(path.read_text())
    required_true = (
        "all_three_broad_mass_series_nonnegative",
        "broad_masses_sum_to_aggregate_trace_through_reconstruction",
        "broad_targets_sum_to_aggregate_trace_target",
    )
    if report.get("status") != "complete_exact_broad_channel_modular_reconstruction":
        raise ValueError("broad-channel report is not a complete exact reconstruction")
    if report.get("p") != p:
        raise ValueError("broad-channel report has the wrong prime")
    if any(report.get(key) is not True for key in required_true):
        raise ArithmeticError("broad-channel report failed an exact validation")
    rows = report.get("channels")
    if not isinstance(rows, dict) or set(rows) != set(CHANNELS):
        raise ValueError("broad-channel report has the wrong channel set")

    dimensions: dict[str, int] = {}
    masses: dict[str, list[Fraction]] = {}
    targets: dict[str, Fraction] = {}
    lengths: set[int] = set()
    for channel in CHANNELS:
        row = rows[channel]
        if not isinstance(row, dict):
            raise ValueError(f"invalid broad-channel row for {channel}")
        if row.get("all_reconstructed_broad_masses_nonnegative") is not True:
            raise ArithmeticError(f"{channel} reconstruction is not nonnegative")
        dimension = int(row["dimension"])
        if dimension <= 0:
            raise ValueError(f"invalid broad dimension for {channel}")
        values = row.get("broad_raw_mass_coefficients")
        if not isinstance(values, list):
            raise ValueError(f"{channel} lacks broad raw masses")
        parsed = [Fraction(value) for value in values]
        if any(value < 0 for value in parsed):
            raise ArithmeticError(f"{channel} contains a negative broad mass")
        dimensions[channel] = dimension
        masses[channel] = parsed
        targets[channel] = Fraction(row["half_cusp_average_harmonic_target"])
        lengths.add(len(parsed))
    if len(lengths) != 1 or not lengths or next(iter(lengths)) < 1:
        raise ValueError("broad-channel coefficient arrays have inconsistent lengths")
    available_through = next(iter(lengths)) - 1
    if available_through != int(report["coefficients_reconstructed_through"]):
        raise ValueError("broad-channel coefficient limit is inconsistent")

    aggregate_target = Fraction(report["aggregate_trace_target"])
    weighted_target = sum(
        (dimensions[channel] * targets[channel] for channel in CHANNELS),
        Fraction(),
    )
    if weighted_target != aggregate_target:
        raise ArithmeticError("broad transformed targets do not sum to the aggregate")
    return dimensions, masses, targets, aggregate_target, report


def json_ready(value: object) -> object:
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--broad-reconstruction", type=Path, required=True)
    parser.add_argument("--scalar-reconstruction", type=Path, required=True)
    parser.add_argument("--exact-row-directory", type=Path, required=True)
    parser.add_argument(
        "--row-suffix", default="_qrows_exact_e800_20260828.txt"
    )
    parser.add_argument("--coefficient-through", type=int)
    parser.add_argument("--case", action="append")
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--esolver", type=Path, default=Path("/usr/bin/esolver"))
    parser.add_argument(
        "--qsopt-library-directory",
        type=Path,
        default=Path("/usr/lib/x86_64-linux-gnu"),
    )
    parser.add_argument("--raw-unconditioned", action="store_true")
    args = parser.parse_args()
    if args.p != 11:
        parser.error("the exact reconstruction is specific to p=11")
    args.output_directory.mkdir(parents=True, exist_ok=True)

    dimensions, broad_masses, broad_targets, aggregate_target, broad_report = (
        load_broad_reconstruction(args.broad_reconstruction, args.p)
    )
    counts = load_scalar_reconstruction(args.scalar_reconstruction)
    recorded_scalar = broad_report.get("inputs", {}).get("scalar_reconstruction", {})
    if not isinstance(recorded_scalar, dict) or recorded_scalar.get("sha256") != sha256(
        args.scalar_reconstruction
    ):
        raise ArithmeticError(
            "the scalar reconstruction differs from the broad-channel input"
        )

    available_through = min(
        len(counts) - 1,
        *(len(broad_masses[channel]) - 1 for channel in CHANNELS),
    )
    coefficient_through = (
        available_through
        if args.coefficient_through is None
        else args.coefficient_through
    )
    if not 92 <= coefficient_through <= available_through:
        parser.error(f"--coefficient-through must lie in 92..{available_through}")
    counts = counts[: coefficient_through + 1]
    broad_masses = {
        channel: values[: coefficient_through + 1]
        for channel, values in broad_masses.items()
    }

    paths = {
        channel: args.exact_row_directory
        / f"p11_affine_{channel}{args.row_suffix}"
        for channel in CHANNELS
    }
    reductions = {}
    for channel, path in paths.items():
        base, matrix, target_base, target = load_qrows(path)
        if len(base) <= coefficient_through:
            raise ValueError(f"{channel} q-row export is too short")
        reductions[channel] = (
            base[: coefficient_through + 1],
            matrix[: coefficient_through + 1],
            target_base,
            target,
        )

    cases = component_cases(args.p)
    expected_dimensions = {
        channel: sum(
            int(case["component_count"])
            * int(case["representation_dimension"])
            for case in cases
            if str(case["channel"]) == channel
        )
        for channel in CHANNELS
    }
    if dimensions != expected_dimensions:
        raise ArithmeticError(
            f"reported broad dimensions {dimensions} differ from {expected_dimensions}"
        )
    if args.case:
        requested = set(args.case)
        unknown = requested - {str(case["name"]) for case in cases}
        if unknown:
            parser.error(f"unknown cases: {sorted(unknown)}")
        selected_indices = [
            index
            for index, case in enumerate(cases)
            if str(case["name"]) in requested
        ]
    else:
        selected_indices = list(range(len(cases)))

    report: dict[str, object] = {
        "experiment": "r1_p11_broad_endpoint_qsopt",
        "status": "running",
        "p": args.p,
        "method": (
            "exact channel-resolved shell moments; unique modular reconstruction; "
            "separate kernel/low/high raw-mass and transformed-target conservation; "
            "rational QSopt_ex primal/dual certificates"
        ),
        "coefficient_through": coefficient_through,
        "broad_reconstruction": {
            "path": args.broad_reconstruction,
            "sha256": sha256(args.broad_reconstruction),
            "profile_fixed_through": broad_report["profile_fixed_through"],
            "profile_available_through": broad_report["profile_available_through"],
            "held_out_coefficients_per_channel": (
                int(broad_report["profile_available_through"])
                - int(broad_report["profile_fixed_through"])
            ),
        },
        "scalar_reconstruction": {
            "path": args.scalar_reconstruction,
            "sha256": sha256(args.scalar_reconstruction),
        },
        "broad_dimensions": dimensions,
        "broad_targets": broad_targets,
        "aggregate_trace_target": aggregate_target,
        "harmonic_qrows": {
            channel: {"path": path, "sha256": sha256(path)}
            for channel, path in paths.items()
        },
        "exact_conditioning": not args.raw_unconditioned,
        "cases": [],
    }
    report_path = args.output_directory / "report.json"
    started_all = time.monotonic()
    intervals_by_case: dict[str, tuple[Fraction, Fraction]] = {}

    for distinguished_index in selected_indices:
        case = cases[distinguished_index]
        case_name = str(case["name"])
        print(f"{case_name}: building broad-conservation model", flush=True)
        started_case = time.monotonic()
        raw_model, representatives = build_broad_channel_conserved_model(
            reductions,
            cases,
            distinguished_index,
            counts,
            broad_masses,
            broad_targets,
            args.p,
        )
        model = raw_model if args.raw_unconditioned else condition_model(raw_model)
        case_row: dict[str, object] = {
            "case": case_name,
            "channel": str(case["channel"]),
            "component_count": int(case["component_count"]),
            "representation_dimension": int(case["representation_dimension"]),
            "variable_count": len(model.target),
            "constraint_count": len(model.constraints),
            "fixed_checks": list(model.fixed_checks),
            "symmetry_representatives": representatives,
            "bounds": {},
        }
        endpoints: dict[str, Fraction] = {}
        for sense in ("minimum", "maximum"):
            print(f"{case_name} {sense}: solving exactly", flush=True)
            lp_path = args.output_directory / f"{case_name}_{sense}.lp"
            solution_path = args.output_directory / f"{case_name}_{sense}.sol"
            log_path = args.output_directory / f"{case_name}_{sense}.log"
            objective = write_lp(lp_path, model, sense)
            started_solve = time.monotonic()
            log = run_qsopt(
                args.esolver,
                args.qsopt_library_directory,
                lp_path,
                solution_path,
            )
            log_path.write_text(log)
            certificate = verify_certificate(model, objective, solution_path)
            solver_value = Fraction(str(certificate["solver_objective"]))
            variable_value = solver_value if sense == "minimum" else -solver_value
            endpoint = model.target_base + variable_value
            endpoints[sense] = endpoint
            case_row["bounds"][sense] = {
                "status": "exact_qsopt_primal_dual_certified",
                "endpoint": endpoint,
                "endpoint_decimal": float(endpoint),
                "solver_seconds": time.monotonic() - started_solve,
                "certificate": certificate,
                "lp": lp_path,
                "lp_sha256": sha256(lp_path),
                "solution": solution_path,
                "solution_sha256": sha256(solution_path),
                "log": log_path,
                "log_sha256": sha256(log_path),
                "solver_log_tail": log.splitlines()[-12:],
            }
            report_path.write_text(json.dumps(json_ready(report), indent=2) + "\n")
        if endpoints["minimum"] > endpoints["maximum"]:
            raise ArithmeticError(f"reversed exact interval for {case_name}")
        intervals_by_case[case_name] = (
            endpoints["minimum"],
            endpoints["maximum"],
        )
        case_row["elapsed_seconds"] = time.monotonic() - started_case
        report["cases"].append(case_row)
        report_path.write_text(json.dumps(json_ready(report), indent=2) + "\n")

    if len(selected_indices) == len(cases):
        channel_variance: dict[str, dict[str, object]] = {}
        within_channel_variance = Fraction()
        all_intervals: list[tuple[Fraction, Fraction, int]] = []
        for channel in CHANNELS:
            intervals: list[tuple[Fraction, Fraction, int]] = []
            for case in cases:
                if str(case["channel"]) != channel:
                    continue
                lower, upper = intervals_by_case[str(case["name"])]
                weight = int(case["representation_dimension"])
                intervals.extend(
                    [(lower, upper, weight)] * int(case["component_count"])
                )
            total = dimensions[channel] * broad_targets[channel]
            variance = maximum_weighted_variance_exact(intervals, total)
            within_channel_variance += Fraction(variance["normalized_variance_max"])
            all_intervals.extend(intervals)
            channel_variance[channel] = {
                "dimension": dimensions[channel],
                "fixed_target_average": broad_targets[channel],
                "fixed_target_mass": total,
                **variance,
            }

        p = args.p
        n = p * p + 1
        zdim = n * (n - 6) // 8
        global_mean = aggregate_target / zdim
        between_channel_variance = sum(
            (
                dimensions[channel]
                * (broad_targets[channel] - global_mean) ** 2
                for channel in CHANNELS
            ),
            Fraction(),
        )
        target_variance = within_channel_variance + between_channel_variance
        lbar = Fraction(8 * (n - 2), n - 6)
        spherical = Fraction(8 * n, n + 4)
        scale = -Fraction(zdim) * (lbar - spherical) / aggregate_target
        phi_variance = scale**2 * target_variance
        phi_endpoint_values = [
            spherical - scale * endpoint
            for lower, upper, _weight in all_intervals
            for endpoint in (lower, upper)
        ]
        phi_lower = min(phi_endpoint_values)
        phi_upper = max(phi_endpoint_values)
        exact_threshold = Fraction(n, 2) * (lbar - 6) ** 2
        strong_threshold = Fraction(2 * n)
        report["variance_bound"] = {
            "method": (
                "exact interval-box relaxation with each broad-channel target "
                "mass fixed separately; within/between weighted variance identity"
            ),
            "global_target_mean": global_mean,
            "channel_bounds": channel_variance,
            "within_channel_target_variance_max": within_channel_variance,
            "between_channel_target_variance_exact": between_channel_variance,
            "total_target_variance_max": target_variance,
            "poisson_scale_from_trace": scale,
            "Phi_frobenius_variance_max": phi_variance,
            "Phi_eigenvalue_lower_bound": phi_lower,
            "Phi_eigenvalue_upper_bound": phi_upper,
            "principal_spectral_floor_margin": phi_lower - 6,
            "closes_principal_spectral_floor_at_p11": phi_lower >= 6,
            "R1_exact_threshold": exact_threshold,
            "R1_exact_margin": exact_threshold - phi_variance,
            "R1_strong_n_over_12_threshold": strong_threshold,
            "R1_strong_margin": strong_threshold - phi_variance,
            "closes_exact_R1_at_p11": phi_variance <= exact_threshold,
            "closes_strong_R1_at_p11": phi_variance <= strong_threshold,
        }

    report["status"] = "complete_exact_broad_qsopt_certified"
    report["elapsed_seconds"] = time.monotonic() - started_all
    report_path.write_text(json.dumps(json_ready(report), indent=2) + "\n")
    print(json.dumps(json_ready(report), indent=2), flush=True)


if __name__ == "__main__":
    main()
