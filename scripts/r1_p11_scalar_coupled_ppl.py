#!/usr/bin/env python3
"""Solve the p=11 scalar/harmonic LP with exact PPL certificates.

PPL's MIP simplex works directly over GMP rationals after each input row is
cleared to primitive integer coefficients.  For every endpoint this script
solves both the primal and its explicit dual, then checks primal feasibility,
dual feasibility, stationarity, and strong duality again with ``Fraction``.
No floating-point solver result is used as proof evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Iterable

try:
    import ppl
except ImportError as error:  # pragma: no cover - environment diagnostic
    raise SystemExit(
        "python-ppl is required; set PYTHONPATH to its dist-packages directory"
    ) from error

from r1_p11_scalar_coupled_exact_lp import (
    ExactModel,
    build_model,
    build_scalar_trace_budget_model,
)
from r1_p11_scalar_coupled_lp import component_cases, load_rows


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lcm(values: Iterable[int]) -> int:
    answer = 1
    for value in values:
        answer = math.lcm(answer, value)
    return answer


def primitive_integer_vector(
    values: Iterable[Fraction],
) -> tuple[tuple[int, ...], Fraction]:
    """Return primitive integers ``a`` and positive scale with a=scale*values."""
    values = tuple(values)
    denominator_lcm = lcm(value.denominator for value in values)
    integers = tuple(
        value.numerator * (denominator_lcm // value.denominator)
        for value in values
    )
    divisor = 0
    for value in integers:
        divisor = math.gcd(divisor, abs(value))
    if divisor == 0:
        return integers, Fraction(1)
    return (
        tuple(value // divisor for value in integers),
        Fraction(denominator_lcm, divisor),
    )


def canonical_rows(
    model: ExactModel,
) -> tuple[tuple[tuple[Fraction, ...], Fraction, str], ...]:
    """Convert constraints to ``a*x >= b``; preserve equality rows if any."""
    rows: list[tuple[tuple[Fraction, ...], Fraction, str]] = []
    for constraint in model.constraints:
        if constraint.sense == ">=":
            rows.append((constraint.coefficients, constraint.rhs, constraint.name))
        elif constraint.sense == "<=":
            rows.append(
                (
                    tuple(-value for value in constraint.coefficients),
                    -constraint.rhs,
                    constraint.name,
                )
            )
        elif constraint.sense == "=":
            rows.append((constraint.coefficients, constraint.rhs, constraint.name + "_ge"))
            rows.append(
                (
                    tuple(-value for value in constraint.coefficients),
                    -constraint.rhs,
                    constraint.name + "_le",
                )
            )
        else:
            raise ValueError(f"unknown constraint sense: {constraint.sense}")
    return tuple(rows)


def ppl_expression(values: Iterable[int], constant: int = 0) -> ppl.Linear_Expression:
    return ppl.Linear_Expression(tuple(int(value) for value in values), int(constant))


def add_rational_relation(
    problem: ppl.MIP_Problem,
    coefficients: tuple[Fraction, ...],
    rhs: Fraction,
    relation: str,
) -> None:
    integers, _scale = primitive_integer_vector((*coefficients, -rhs))
    expression = ppl_expression(integers[:-1], integers[-1])
    if relation == ">=":
        problem.add_constraint(expression >= 0)
    elif relation == "=":
        problem.add_constraint(expression == 0)
    else:  # pragma: no cover - internal API guard
        raise ValueError(relation)


def point_coordinates(point: ppl.Generator, dimension: int) -> tuple[Fraction, ...]:
    divisor = int(point.divisor())
    coefficients = tuple(int(value) for value in point.coefficients())
    return tuple(
        Fraction(coefficients[index] if index < len(coefficients) else 0, divisor)
        for index in range(dimension)
    )


def solve_primal(
    model: ExactModel,
    objective: tuple[Fraction, ...],
) -> tuple[str, Fraction | None, tuple[Fraction, ...] | None]:
    problem = ppl.MIP_Problem(len(objective))
    for coefficients, rhs, _name in canonical_rows(model):
        add_rational_relation(problem, coefficients, rhs, ">=")
    objective_integers, objective_scale = primitive_integer_vector(objective)
    problem.set_objective_function(ppl_expression(objective_integers))
    problem.set_optimization_mode("minimization")
    status = str(problem.solve()["status"])
    if status != "optimized":
        return status, None, None
    solver_value = Fraction(str(problem.optimal_value())) / objective_scale
    point = point_coordinates(problem.optimizing_point(), len(objective))
    return status, solver_value, point


def solve_dual(
    model: ExactModel,
    objective: tuple[Fraction, ...],
) -> tuple[str, Fraction | None, tuple[Fraction, ...] | None]:
    rows = canonical_rows(model)
    problem = ppl.MIP_Problem(len(rows))
    for index in range(len(rows)):
        coefficients = tuple(
            Fraction(1) if row_index == index else Fraction()
            for row_index in range(len(rows))
        )
        add_rational_relation(problem, coefficients, Fraction(), ">=")
    for column, expected in enumerate(objective):
        coefficients = tuple(row[0][column] for row in rows)
        add_rational_relation(problem, coefficients, expected, "=")
    dual_objective = tuple(row[1] for row in rows)
    objective_integers, objective_scale = primitive_integer_vector(dual_objective)
    problem.set_objective_function(ppl_expression(objective_integers))
    problem.set_optimization_mode("maximization")
    status = str(problem.solve()["status"])
    if status != "optimized":
        return status, None, None
    solver_value = Fraction(str(problem.optimal_value())) / objective_scale
    point = point_coordinates(problem.optimizing_point(), len(rows))
    return status, solver_value, point


def solve_unbounded_ray(
    model: ExactModel,
    objective: tuple[Fraction, ...],
) -> tuple[str, tuple[Fraction, ...] | None]:
    """Find an exact recession direction normalized by c*d <= -1."""
    problem = ppl.MIP_Problem(len(objective))
    for coefficients, _rhs, _name in canonical_rows(model):
        add_rational_relation(problem, coefficients, Fraction(), ">=")
    # Homogeneity lets every improving direction be normalized exactly to
    # c*d=-1.  The equality is both a sharper certificate and substantially
    # easier for the exact simplex than the unbounded half-space c*d<=-1.
    add_rational_relation(problem, objective, Fraction(-1), "=")
    problem.set_objective_function(ppl_expression((0,) * len(objective)))
    problem.set_optimization_mode("minimization")
    status = str(problem.solve()["status"])
    if status != "optimized":
        return status, None
    return status, point_coordinates(problem.optimizing_point(), len(objective))


def verify_unbounded_ray(
    model: ExactModel,
    objective: tuple[Fraction, ...],
    ray: tuple[Fraction, ...],
) -> dict[str, object]:
    rows = canonical_rows(model)
    for coefficients, _rhs, name in rows:
        derivative = sum(
            (coefficient * value for coefficient, value in zip(coefficients, ray)),
            Fraction(),
        )
        if derivative < 0:
            raise ArithmeticError(
                f"recession direction violates row {name}: {derivative} < 0"
            )
    objective_derivative = sum(
        (coefficient * value for coefficient, value in zip(objective, ray)),
        Fraction(),
    )
    if objective_derivative != -1:
        raise ArithmeticError(
            f"recession objective is not normalized to -1: {objective_derivative}"
        )
    return {
        "objective_derivative": str(objective_derivative),
        "recession_constraints_verified": len(rows),
        "nonzero_ray_variables": sum(value != 0 for value in ray),
        "ray": [str(value) for value in ray],
    }


def verify_certificate(
    model: ExactModel,
    objective: tuple[Fraction, ...],
    primal_value: Fraction,
    primal: tuple[Fraction, ...],
    dual_value: Fraction,
    dual: tuple[Fraction, ...],
) -> dict[str, object]:
    rows = canonical_rows(model)
    computed_primal = sum(
        (coefficient * value for coefficient, value in zip(objective, primal)),
        Fraction(),
    )
    if computed_primal != primal_value:
        raise ArithmeticError(
            f"primal objective mismatch: {computed_primal} != {primal_value}"
        )
    for coefficients, rhs, name in rows:
        lhs = sum(
            (coefficient * value for coefficient, value in zip(coefficients, primal)),
            Fraction(),
        )
        if lhs < rhs:
            raise ArithmeticError(f"primal row {name} fails: {lhs} < {rhs}")
    if any(value < 0 for value in dual):
        raise ArithmeticError("dual certificate has a negative multiplier")
    for column, expected in enumerate(objective):
        actual = sum(
            (dual[row_index] * rows[row_index][0][column] for row_index in range(len(rows))),
            Fraction(),
        )
        if actual != expected:
            raise ArithmeticError(
                f"dual stationarity fails at column {column}: {actual} != {expected}"
            )
    computed_dual = sum(
        (dual[index] * row[1] for index, row in enumerate(rows)), Fraction()
    )
    if computed_dual != dual_value:
        raise ArithmeticError(
            f"dual objective mismatch: {computed_dual} != {dual_value}"
        )
    if primal_value != dual_value:
        raise ArithmeticError(
            f"strong duality fails: {primal_value} != {dual_value}"
        )
    return {
        "objective": str(primal_value),
        "primal_constraints_verified": len(rows),
        "dual_stationarity_equations_verified": len(objective),
        "nonzero_primal_variables": sum(value != 0 for value in primal),
        "nonzero_dual_variables": sum(value != 0 for value in dual),
        "primal": [str(value) for value in primal],
        "dual": [str(value) for value in dual],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scalar-qrows", type=Path, required=True)
    parser.add_argument("--affine-prefix", type=Path, required=True)
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--case", action="append")
    parser.add_argument("--sense", choices=("minimum", "maximum"), action="append")
    parser.add_argument(
        "--ray-only",
        action="store_true",
        help="skip optimization and certify a normalized recession direction",
    )
    parser.add_argument(
        "--joint-trace-budget",
        action="store_true",
        help="couple all constituent modular forms through one shell trace budget",
    )
    parser.add_argument(
        "--parity-fourth-moment",
        action="store_true",
        help="tighten the trace budget by the exact scaled-coordinate parity DP",
    )
    parser.add_argument(
        "--scalar-half-target-rows",
        type=Path,
        help="exact ten-component scalar coefficient immediately after the cusp-1/2 gap",
    )
    parser.add_argument(
        "--scalar-half-target-first",
        type=Fraction,
        help="fixed first rational component of that cyclotomic coefficient",
    )
    parser.add_argument(
        "--scalar-fixed-count",
        metavar="EXPONENT=COUNT",
        action="append",
        default=[],
        help="fix an additional exact ordinary-theta coefficient (repeatable)",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if (args.scalar_half_target_rows is None) != (args.scalar_half_target_first is None):
        parser.error("--scalar-half-target-rows and --scalar-half-target-first are required together")
    scalar_fixed_counts: dict[int, Fraction] = {}
    for item in args.scalar_fixed_count:
        exponent_text, separator, count_text = item.partition("=")
        if not separator:
            parser.error(f"invalid --scalar-fixed-count {item!r}; expected EXPONENT=COUNT")
        try:
            exponent = int(exponent_text)
            count = Fraction(count_text)
        except (ValueError, ZeroDivisionError):
            parser.error(f"invalid --scalar-fixed-count {item!r}; expected integer values")
        if exponent < 0 or count.denominator != 1 or count < 0:
            parser.error(f"invalid --scalar-fixed-count {item!r}; values must be nonnegative integers")
        if exponent in scalar_fixed_counts and scalar_fixed_counts[exponent] != count:
            parser.error(f"conflicting fixed counts for scalar exponent {exponent}")
        scalar_fixed_counts[exponent] = count

    scalar_rows = load_rows(args.scalar_qrows)
    scalar_half_target_rows = (
        load_rows(args.scalar_half_target_rows)
        if args.scalar_half_target_rows is not None
        else None
    )
    cases = component_cases(args.p)
    if args.case:
        requested = set(args.case)
        unknown = requested - {str(case["name"]) for case in cases}
        if unknown:
            raise ValueError(f"unknown cases: {sorted(unknown)}")
        cases = [case for case in cases if str(case["name"]) in requested]
    senses = args.sense or ["minimum", "maximum"]
    report: dict[str, object] = {
        "experiment": "r1_p11_scalar_coupled_ppl",
        "status": "exact_ppl_certified",
        "p": args.p,
        "scalar_qrows": str(args.scalar_qrows),
        "scalar_qrows_sha256": sha256(args.scalar_qrows),
        "cases": [],
    }
    for case in cases:
        case_name = str(case["name"])
        channel = str(case["channel"])
        harmonic_path = Path(
            f"{args.affine_prefix}{channel}_qrows_exact_20260827.txt"
        )
        representatives: tuple[dict[str, int | str], ...] = ()
        if args.joint_trace_budget:
            harmonic_rows_by_channel = {
                name: load_rows(
                    Path(f"{args.affine_prefix}{name}_qrows_exact_20260827.txt")
                )
                for name in ("circle-kernel", "circle-low", "circle-high")
            }
            model, representatives = build_scalar_trace_budget_model(
                scalar_rows,
                harmonic_rows_by_channel,
                args.p,
                case_name,
                parity_fourth_moment=args.parity_fourth_moment,
                scalar_half_target_rows=scalar_half_target_rows,
                scalar_half_target_first=args.scalar_half_target_first,
                scalar_fixed_counts=scalar_fixed_counts,
            )
        else:
            model = build_model(
                scalar_rows,
                load_rows(harmonic_path),
                args.p,
                int(case["representation_dimension"]),
                scalar_fixed_counts=scalar_fixed_counts,
            )
        case_report: dict[str, object] = {
            "case": case_name,
            "channel": channel,
            "component_count": int(case["component_count"]),
            "representation_dimension": int(case["representation_dimension"]),
            "harmonic_qrows": str(harmonic_path),
            "harmonic_qrows_sha256": sha256(harmonic_path),
            "variable_count": len(model.target),
            "constraint_count": len(model.constraints),
            "fixed_checks": list(model.fixed_checks),
            "joint_trace_budget": args.joint_trace_budget,
            "parity_fourth_moment": args.parity_fourth_moment,
            "scalar_half_target_rows": (
                str(args.scalar_half_target_rows)
                if args.scalar_half_target_rows is not None
                else None
            ),
            "scalar_half_target_rows_sha256": (
                sha256(args.scalar_half_target_rows)
                if args.scalar_half_target_rows is not None
                else None
            ),
            "scalar_half_target_first": (
                str(args.scalar_half_target_first)
                if args.scalar_half_target_first is not None
                else None
            ),
            "scalar_fixed_counts": {
                str(exponent): str(count)
                for exponent, count in sorted(scalar_fixed_counts.items())
            },
            "symmetry_representatives": representatives,
            "bounds": {},
        }
        for sense in senses:
            objective = (
                model.target
                if sense == "minimum"
                else tuple(-value for value in model.target)
            )
            if args.ray_only:
                print(f"{case_name} {sense}: solving exact normalized ray", flush=True)
                started = time.monotonic()
                ray_status, ray = solve_unbounded_ray(model, objective)
                elapsed_ray = time.monotonic() - started
                if ray_status == "unfeasible":
                    case_report["bounds"][sense] = {
                        "status": "exact_no_improving_recession_ray",
                        "ray_status": ray_status,
                        "ray_seconds": elapsed_ray,
                        "recession_constraints_checked": len(canonical_rows(model)),
                    }
                    continue
                if ray_status != "optimized" or ray is None:
                    raise ArithmeticError(
                        f"{case_name} {sense} ray status is {ray_status}"
                    )
                case_report["bounds"][sense] = {
                    "status": "exact_unbounded_ray",
                    "ray_status": ray_status,
                    "ray_seconds": elapsed_ray,
                    "unbounded_certificate": verify_unbounded_ray(
                        model, objective, ray
                    ),
                }
                continue
            print(f"{case_name} {sense}: solving exact primal", flush=True)
            started = time.monotonic()
            primal_status, primal_value, primal = solve_primal(model, objective)
            elapsed_primal = time.monotonic() - started
            if primal_status == "unbounded":
                print(
                    f"{case_name} {sense}: exact primal is unbounded in "
                    f"{elapsed_primal:.2f}s; certifying a recession ray",
                    flush=True,
                )
                ray_status, ray = solve_unbounded_ray(model, objective)
                if ray_status != "optimized" or ray is None:
                    raise ArithmeticError(
                        f"{case_name} {sense} failed to certify unboundedness: "
                        f"ray status {ray_status}"
                    )
                case_report["bounds"][sense] = {
                    "status": "exact_unbounded",
                    "primal_status": primal_status,
                    "primal_seconds": elapsed_primal,
                    "unbounded_certificate": verify_unbounded_ray(
                        model, objective, ray
                    ),
                }
                continue
            if primal_status != "optimized" or primal_value is None or primal is None:
                raise ArithmeticError(
                    f"{case_name} {sense} primal status is {primal_status}"
                )
            print(
                f"{case_name} {sense}: primal {primal_value} in {elapsed_primal:.2f}s; "
                "solving exact dual",
                flush=True,
            )
            started = time.monotonic()
            dual_status, dual_value, dual = solve_dual(model, objective)
            elapsed_dual = time.monotonic() - started
            if dual_status != "optimized" or dual_value is None or dual is None:
                raise ArithmeticError(
                    f"{case_name} {sense} dual status is {dual_status}"
                )
            certificate = verify_certificate(
                model,
                objective,
                primal_value,
                primal,
                dual_value,
                dual,
            )
            endpoint = (
                model.target_base + primal_value
                if sense == "minimum"
                else model.target_base - primal_value
            )
            print(
                f"{case_name} {sense}: endpoint {endpoint}; exact certificate verified "
                f"in {elapsed_dual:.2f}s",
                flush=True,
            )
            case_report["bounds"][sense] = {
                "target_endpoint": str(endpoint),
                "primal_status": primal_status,
                "dual_status": dual_status,
                "primal_seconds": elapsed_primal,
                "dual_seconds": elapsed_dual,
                "certificate": certificate,
            }
        report["cases"].append(case_report)

    payload = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
