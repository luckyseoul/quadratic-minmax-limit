#!/usr/bin/env python3
"""Certify p=11 R1 target bounds with exact shellwise trace conservation."""
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
    build_shellwise_conserved_model,
    load_qrows,
    run_qsopt,
    verify_certificate,
    write_lp,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_trace_reconstruction(
    path: Path, p: int
) -> tuple[list[int], list[Fraction], list[Fraction], dict[str, object]]:
    report = json.loads(path.read_text())
    if report.get("status") != "complete_exact_trace_modular_reconstruction":
        raise ValueError("trace report is not a complete exact reconstruction")
    if report.get("p") != p:
        raise ValueError("trace report has the wrong prime")
    required = (
        "out_of_sample_profile_trace_coefficients_match",
        "known_raw_trace_coefficients_match",
        "all_reconstructed_raw_trace_coefficients_nonnegative",
    )
    if any(report.get(item) is not True for item in required):
        raise ArithmeticError("trace report failed an exact validation")
    counts = report.get("ordinary_theta_coefficients")
    harmonic_values = report.get("harmonic_trace_coefficients")
    raw_values = report.get("raw_trace_coefficients")
    if not isinstance(counts, list) or any(
        not isinstance(value, int) or value < 0 for value in counts
    ):
        raise ValueError("trace report has invalid ordinary theta coefficients")
    if not isinstance(harmonic_values, list) or not isinstance(raw_values, list):
        raise ValueError("trace report lacks trace coefficient arrays")
    harmonic = [Fraction(value) for value in harmonic_values]
    raw = [Fraction(value) for value in raw_values]
    if not len(counts) == len(harmonic) == len(raw):
        raise ValueError("trace reconstruction arrays have inconsistent lengths")

    n = p * p + 1
    d = n // 2
    zdim = n * (n - 6) // 8
    for exponent, (count, hvalue, rvalue) in enumerate(
        zip(counts, harmonic, raw)
    ):
        radius_squared = Fraction(exponent, 2 * p)
        radial = Fraction(2 * count, d * (d + 2)) * radius_squared**2
        if hvalue + zdim * radial != rvalue or rvalue < 0:
            raise ArithmeticError(f"raw/harmonic trace identity fails at {exponent}")
    return counts, harmonic, raw, report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trace-reconstruction", type=Path, required=True)
    parser.add_argument("--affine-prefix", type=Path, required=True)
    parser.add_argument(
        "--affine-suffix",
        default="_qrows_exact_e800_20260828.txt",
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

    counts, harmonic_trace, _raw_trace, trace_report = load_trace_reconstruction(
        args.trace_reconstruction, args.p
    )
    available_through = len(counts) - 1
    coefficient_through = (
        available_through
        if args.coefficient_through is None
        else args.coefficient_through
    )
    if not 28 <= coefficient_through <= available_through:
        parser.error(f"--coefficient-through must lie in 28..{available_through}")
    counts = counts[: coefficient_through + 1]
    harmonic_trace = harmonic_trace[: coefficient_through + 1]

    paths = {
        channel: Path(f"{args.affine_prefix}{channel}{args.affine_suffix}")
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
        "experiment": "r1_p11_trace_endpoint_qsopt",
        "status": "running",
        "p": args.p,
        "method": (
            "exact profile fourth moments; uniquely reconstructed scalar and "
            "trace modular forms; shellwise positive raw-channel mass; exact "
            "mass conservation; rational QSopt_ex primal/dual certificates"
        ),
        "trace_reconstruction": str(args.trace_reconstruction),
        "trace_reconstruction_sha256": sha256(args.trace_reconstruction),
        "trace_profile_fixed_through": trace_report["profile_fixed_through"],
        "trace_first_full_rank_exponent": trace_report["first_full_rank_exponent"],
        "coefficient_through": coefficient_through,
        "harmonic_qrows": {
            channel: {"path": str(path), "sha256": sha256(path)}
            for channel, path in paths.items()
        },
        "exact_conditioning": not args.raw_unconditioned,
        "cases": [],
    }
    report_path = args.output_directory / "report.json"
    started_all = time.monotonic()

    for distinguished_index in selected_indices:
        case = cases[distinguished_index]
        case_name = str(case["name"])
        print(f"{case_name}: building exact conservation model", flush=True)
        started_case = time.monotonic()
        raw_model, representatives = build_shellwise_conserved_model(
            reductions,
            cases,
            distinguished_index,
            counts,
            harmonic_trace,
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
            case_row["bounds"][sense] = {
                "status": "exact_qsopt_primal_dual_certified",
                "endpoint": str(endpoint),
                "endpoint_decimal": float(endpoint),
                "solver_seconds": time.monotonic() - started_solve,
                "certificate": certificate,
                "lp": str(lp_path),
                "lp_sha256": sha256(lp_path),
                "solution": str(solution_path),
                "solution_sha256": sha256(solution_path),
                "log": str(log_path),
                "log_sha256": sha256(log_path),
                "solver_log_tail": log.splitlines()[-12:],
            }
            report_path.write_text(json.dumps(report, indent=2) + "\n")
        case_row["elapsed_seconds"] = time.monotonic() - started_case
        report["cases"].append(case_row)
        report_path.write_text(json.dumps(report, indent=2) + "\n")

    report["status"] = "complete_exact_qsopt_certified"
    report["elapsed_seconds"] = time.monotonic() - started_all
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
