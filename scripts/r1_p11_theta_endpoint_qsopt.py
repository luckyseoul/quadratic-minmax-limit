#!/usr/bin/env python3
"""Certify p=11 harmonic target bounds using exact ordinary-theta counts.

This joins the exact profile/CRT ordinary-theta prefix to the strongest
scalar/harmonic outer relaxation: all constituent channels share the scalar
theta series and a shellwise trace budget, with the coordinate-parity fourth
moment bound.  Every reported endpoint is accompanied by an independently
checked exact QSopt_ex primal/dual certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from fractions import Fraction
from pathlib import Path

from r1_p11_scalar_coupled_exact_lp import (
    build_scalar_trace_budget_model,
    condition_model,
)
from r1_p11_scalar_coupled_lp import CHANNELS, component_cases, load_rows
from r1_p11_trace_coupled_exact_lp import (
    run_qsopt,
    verify_certificate,
    write_lp,
)


PROVED_PREFIX_THROUGH = 28


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_exact_theta_prefix(path: Path, p: int) -> tuple[list[int], dict[int, Fraction]]:
    report = json.loads(path.read_text())
    if report.get("status") != "complete_exact_theta_prefix":
        raise ValueError("theta report is not marked complete_exact_theta_prefix")
    if report.get("p") != p:
        raise ValueError(f"theta report has p={report.get('p')}, expected {p}")
    if report.get("known_coefficients_match") is not True:
        raise ArithmeticError("theta report failed its known-coefficient calibration")
    if report.get("crt_product_exceeds_every_bound") is not True:
        raise ArithmeticError("theta report lacks a sufficient CRT uniqueness bound")
    crt_product = int(report["crt_modulus_product"])
    unrestricted_bound = int(report["maximum_unrestricted_bound"])
    if crt_product <= unrestricted_bound:
        raise ArithmeticError("CRT modulus product does not exceed the reconstruction bound")

    coefficients = report.get("theta_coefficients")
    max_exponent = int(report.get("max_exponent", -1))
    if not isinstance(coefficients, list) or len(coefficients) != max_exponent + 1:
        raise ValueError("theta coefficient list does not match max_exponent")
    if any(not isinstance(value, int) or value < 0 for value in coefficients):
        raise ValueError("theta coefficients must be nonnegative integers")
    for exponent_text, expected in report.get("known_coefficients", {}).items():
        exponent = int(exponent_text)
        if coefficients[exponent] != int(expected):
            raise ArithmeticError(f"known theta coefficient mismatch at exponent {exponent}")

    fixed = {
        exponent: Fraction(value)
        for exponent, value in enumerate(coefficients)
        if exponent > PROVED_PREFIX_THROUGH
    }
    return coefficients, fixed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--theta-report", type=Path, required=True)
    parser.add_argument("--scalar-qrows", type=Path, required=True)
    parser.add_argument("--affine-prefix", type=Path, required=True)
    parser.add_argument(
        "--affine-suffix",
        default="_qrows_exact_20260827.txt",
        help="suffix appended after each harmonic channel name",
    )
    parser.add_argument("--scalar-half-target-rows", type=Path, required=True)
    parser.add_argument("--scalar-half-target-first", type=Fraction, required=True)
    parser.add_argument("--case", action="append")
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--esolver", type=Path, default=Path("/usr/bin/esolver"))
    parser.add_argument(
        "--qsopt-library-directory",
        type=Path,
        default=Path("/usr/lib/x86_64-linux-gnu"),
    )
    parser.add_argument(
        "--raw-unconditioned",
        action="store_true",
        help="solve the equivalent unconditioned rational model",
    )
    args = parser.parse_args()
    if args.p != 11:
        parser.error("the current exact exports and profile engine are specific to p=11")

    args.output_directory.mkdir(parents=True, exist_ok=True)
    theta_coefficients, scalar_fixed_counts = load_exact_theta_prefix(
        args.theta_report, args.p
    )
    scalar_rows = load_rows(args.scalar_qrows)
    scalar_half_target_rows = load_rows(args.scalar_half_target_rows)
    harmonic_paths = {
        channel: Path(f"{args.affine_prefix}{channel}{args.affine_suffix}")
        for channel in CHANNELS
    }
    harmonic_rows = {
        channel: load_rows(path) for channel, path in harmonic_paths.items()
    }

    cases = component_cases(args.p)
    if args.case:
        requested = set(args.case)
        unknown = requested - {str(case["name"]) for case in cases}
        if unknown:
            parser.error(f"unknown cases: {sorted(unknown)}")
        cases = [case for case in cases if str(case["name"]) in requested]

    report: dict[str, object] = {
        "experiment": "r1_p11_theta_endpoint_qsopt",
        "status": "running",
        "p": args.p,
        "method": (
            "exact profile/CRT ordinary-theta prefix; W4-reduced scalar space; "
            "joint constituent trace budget; exact parity fourth-moment bound; "
            "conditioned rational QSopt_ex primal/dual certificates"
        ),
        "theta_report": str(args.theta_report),
        "theta_report_sha256": sha256(args.theta_report),
        "theta_fixed_through": len(theta_coefficients) - 1,
        "theta_fixed_count_count": len(scalar_fixed_counts),
        "scalar_qrows": str(args.scalar_qrows),
        "scalar_qrows_sha256": sha256(args.scalar_qrows),
        "scalar_half_target_rows": str(args.scalar_half_target_rows),
        "scalar_half_target_rows_sha256": sha256(args.scalar_half_target_rows),
        "scalar_half_target_first": str(args.scalar_half_target_first),
        "harmonic_qrows": {
            channel: {"path": str(path), "sha256": sha256(path)}
            for channel, path in harmonic_paths.items()
        },
        "exact_conditioning": not args.raw_unconditioned,
        "cases": [],
    }
    report_path = args.output_directory / "report.json"
    started_all = time.monotonic()

    for case in cases:
        case_name = str(case["name"])
        print(f"{case_name}: building exact joint model", flush=True)
        started_case = time.monotonic()
        raw_model, representatives = build_scalar_trace_budget_model(
            scalar_rows,
            harmonic_rows,
            args.p,
            case_name,
            parity_fourth_moment=True,
            scalar_half_target_rows=scalar_half_target_rows,
            scalar_half_target_first=args.scalar_half_target_first,
            scalar_fixed_counts=scalar_fixed_counts,
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
