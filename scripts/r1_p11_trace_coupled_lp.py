#!/usr/bin/env python3
"""Couple p=11 R1 theta channels by conserved positive shell mass.

For each dual shell let ``A_s,c`` be its harmonic eigenvalue on a
multiplicity-free PSL constituent ``c`` and put

    q_s,c = A_s,c + 2 N_s r_s^2 / (d(d+2)).

The raw quartic operator is positive, so ``q_s,c >= 0``.  Its trace is
known from the trace-harmonic theta series, giving

    sum_c dim(c) q_s,c = tau_s.

Consequently ``q_s,c <= tau_s/dim(c)``.  This script combines those sharp
channel bounds with the exact modular coefficient matrix and the first four
proved shell operators.  Floating LP is used only for reconnaissance; any
closing bound must subsequently be certified over exact rationals.
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from r1_p11_shell_positivity_lp import (
    load_modular_export,
    load_theta_counts,
    parse_gp_vector,
)


CHANNELS = ("circle-kernel", "circle-low", "circle-high")


def parse_named_vector(path: Path, name: str) -> list[Fraction]:
    text = path.read_text()
    match = re.search(rf"{re.escape(name)}\s*=\s*(\[[^\n]*\])", text)
    if match is None:
        raise ValueError(f"{name} not found in {path}")
    payload = match.group(1)[1:-1]
    return [Fraction(token.strip()) for token in payload.split(",")]


def load_affine_qrows(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Load [constant, free q-coordinate row] lines plus a final target line."""
    rows = [
        parse_gp_vector(line)
        for line in path.read_text().splitlines()
        if line.strip()
    ]
    if len(rows) < 2:
        raise ValueError(f"affine q-row export is too short: {path}")
    width = len(rows[-1])
    if width < 2 or any(len(row) != width for row in rows):
        raise ValueError(f"inconsistent affine q-row widths: {path}")
    coefficient_rows = np.asarray(rows[:-1], dtype=np.float64)
    target_row = np.asarray(rows[-1], dtype=np.float64)
    return (
        coefficient_rows[:, 0],
        coefficient_rows[:, 1:],
        float(target_row[0]),
        target_row[1:],
    )


def component_cases(p: int) -> list[dict[str, int | str]]:
    kernel = (p - 1) * (p - 3) // 8
    cases: list[dict[str, int | str]] = [
        {
            "name": "circle-kernel-principal",
            "channel": "circle-kernel",
            "representation_dimension": p * p + 1,
            "component_count": kernel,
        }
    ]
    if p % 4 == 1:
        cases.extend(
            [
                {
                    "name": "circle-low-principal",
                    "channel": "circle-low",
                    "representation_dimension": p * p + 1,
                    "component_count": (p - 1) // 4,
                },
                {
                    "name": "circle-high-Weil",
                    "channel": "circle-high",
                    "representation_dimension": (p * p + 1) // 2,
                    "component_count": 1,
                },
                {
                    "name": "circle-high-principal",
                    "channel": "circle-high",
                    "representation_dimension": p * p + 1,
                    "component_count": (p - 5) // 4,
                },
            ]
        )
    else:
        cases.extend(
            [
                {
                    "name": "circle-low-Weil",
                    "channel": "circle-low",
                    "representation_dimension": (p * p + 1) // 2,
                    "component_count": 1,
                },
                {
                    "name": "circle-low-principal",
                    "channel": "circle-low",
                    "representation_dimension": p * p + 1,
                    "component_count": (p - 3) // 4,
                },
                {
                    "name": "circle-high-principal",
                    "channel": "circle-high",
                    "representation_dimension": p * p + 1,
                    "component_count": (p - 3) // 4,
                },
            ]
        )
    return [case for case in cases if int(case["component_count"]) > 0]


def solve_channel(
    matrix: np.ndarray,
    target: np.ndarray,
    coefficient_base: np.ndarray,
    target_base: float,
    counts: np.ndarray,
    trace_harmonic: list[Fraction],
    p: int,
    channel: str,
    representation_dimension: int,
) -> dict:
    d = (p * p + 1) // 2
    zdim = (p * p + 1) * (p * p - 5) // 8
    limit = min(len(counts), len(trace_harmonic), len(matrix))
    width = matrix.shape[1]
    fixed_through = 2 * (p + 3)

    aub: list[np.ndarray] = []
    bub: list[float] = []
    aeq: list[np.ndarray] = []
    beq: list[float] = []
    tau_min = float("inf")
    fixed_row_max = 0.0
    fixed_raw_violation = 0.0

    for exponent in range(limit):
        row = matrix[exponent]
        base = float(coefficient_base[exponent])
        count = int(counts[exponent])
        radius_sq_exact = Fraction(exponent, 2 * p)
        radial_exact = (
            Fraction(2 * count, d * (d + 2)) * radius_sq_exact**2
        )
        tau_exact = trace_harmonic[exponent] + zdim * radial_exact
        if tau_exact < 0:
            raise ArithmeticError(
                f"negative raw trace mass at exponent {exponent}: {tau_exact}"
            )
        radial = float(radial_exact)
        tau = float(tau_exact)
        if tau > 0:
            tau_min = min(tau_min, tau)

        # The exact affine reduction already imposes every cusp gap and every
        # infinity coefficient through exponent 28.  Checking, rather than
        # re-imposing, these rows avoids catastrophic redundant equalities.
        if exponent <= fixed_through:
            fixed_row_max = max(fixed_row_max, float(np.max(np.abs(row))))
            raw_upper = tau / representation_dimension
            raw_value = base + radial
            fixed_raw_violation = max(
                fixed_raw_violation,
                -raw_value,
                raw_value - raw_upper,
            )
            continue
        if count == 0 or tau_exact == 0:
            value = -radial if count else 0.0
            scale = max(1.0, float(np.max(np.abs(row))), abs(base), abs(value))
            aeq.append(row / scale)
            beq.append((value - base) / scale)
            continue

        raw_upper = float(tau_exact / representation_dimension)
        # 0 <= base + row*y + radial <= raw_upper.
        scale = raw_upper
        aub.append(row / scale)
        bub.append(1.0 - (base + radial) / scale)
        aub.append(-row / scale)
        bub.append((base + radial) / scale)

    A_ub = np.asarray(aub, dtype=np.float64).reshape((-1, width))
    b_ub = np.asarray(bub)
    A_eq = np.asarray(aeq, dtype=np.float64).reshape((-1, width))
    b_eq = np.asarray(beq)
    stacked = np.vstack([array for array in (A_ub, A_eq) if len(array)])
    column_max = np.max(np.abs(stacked), axis=0)
    variable_scale = np.ones(width)
    nonzero = column_max > 0
    variable_scale[nonzero] = 1.0 / column_max[nonzero]
    A_ub_scaled = A_ub * variable_scale
    A_eq_scaled = A_eq * variable_scale
    objective = target * variable_scale

    output: dict[str, object] = {
        "channel": channel,
        "representation_dimension": representation_dimension,
        "rows_used": limit,
        "minimum_raw_trace_mass": tau_min,
        "fixed_through_exponent": fixed_through,
        "fixed_free_row_max_abs": fixed_row_max,
        "fixed_raw_bound_violation": fixed_raw_violation,
    }
    for sense, c in (("minimum", objective), ("maximum", -objective)):
        result = linprog(
            c,
            A_ub=A_ub_scaled,
            b_ub=b_ub,
            A_eq=A_eq_scaled,
            b_eq=b_eq,
            bounds=[(None, None)] * width,
            method="highs",
            options={
                "dual_feasibility_tolerance": 1e-9,
                "primal_feasibility_tolerance": 1e-9,
            },
        )
        record: dict[str, object] = {
            "success": bool(result.success),
            "status": int(result.status),
            "message": result.message,
        }
        if result.success:
            coordinates = result.x * variable_scale
            value = float(target_base + target @ coordinates)
            record.update(
                value=value,
                normalized_eq_residual=float(
                    np.max(np.abs(A_eq @ coordinates - b_eq)) if len(A_eq) else 0.0
                ),
                normalized_ub_residual=float(
                    np.max(A_ub @ coordinates - b_ub) if len(A_ub) else 0.0
                ),
            )
        output[sense] = record
    return output


def trace_coordinates(
    matrix: np.ndarray, trace_harmonic: list[Fraction], limit: int
) -> tuple[np.ndarray, float]:
    solution = np.empty(matrix.shape[1])
    for column in range(matrix.shape[1]):
        candidates = np.flatnonzero(
            np.isclose(matrix[:limit, column], 1.0, atol=1e-12, rtol=0.0)
        )
        pivot = next(
            (
                int(row)
                for row in candidates
                if np.count_nonzero(np.abs(matrix[row]) > 1e-12) == 1
            ),
            None,
        )
        if pivot is None:
            raise ArithmeticError(f"no identity pivot for modular column {column}")
        solution[column] = float(trace_harmonic[pivot])
    residual = float(
        np.max(
            np.abs(
                matrix[:limit] @ solution
                - np.asarray([float(value) for value in trace_harmonic[:limit]])
            )
        )
    )
    return solution, residual


def maximum_weighted_variance(
    intervals: list[tuple[float, float, int]], total: float
) -> dict:
    """Maximize sum weights*(x-total/sum(weights))^2 over box plus sum."""
    expanded: list[tuple[float, float, int]] = []
    for lower, upper, weight_and_count in intervals:
        # The caller repeats components explicitly and passes their irrep dim.
        expanded.append((lower, upper, weight_and_count))
    weights = np.asarray([row[2] for row in expanded], dtype=np.float64)
    lower = np.asarray([row[0] for row in expanded], dtype=np.float64)
    upper = np.asarray([row[1] for row in expanded], dtype=np.float64)
    mean = total / float(weights.sum())
    best_value = -np.inf
    best_point: np.ndarray | None = None
    count = len(expanded)
    for pivot in range(count):
        other = [index for index in range(count) if index != pivot]
        for mask in range(1 << len(other)):
            point = np.empty(count)
            subtotal = 0.0
            for bit, index in enumerate(other):
                point[index] = upper[index] if mask & (1 << bit) else lower[index]
                subtotal += weights[index] * point[index]
            point[pivot] = (total - subtotal) / weights[pivot]
            tolerance = 1e-10 * max(1.0, abs(lower[pivot]), abs(upper[pivot]))
            if point[pivot] < lower[pivot] - tolerance or point[pivot] > upper[pivot] + tolerance:
                continue
            point[pivot] = min(upper[pivot], max(lower[pivot], point[pivot]))
            value = float(np.sum(weights * (point - mean) ** 2))
            if value > best_value:
                best_value = value
                best_point = point.copy()
    if best_point is None:
        raise ArithmeticError("target trace is outside the channel interval box")
    return {
        "normalized_mean": mean,
        "normalized_variance_max": best_value,
        "maximizer": best_point.tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--modular-export", type=Path, required=True)
    parser.add_argument("--theta-output", type=Path, required=True)
    parser.add_argument("--trace-output", type=Path, required=True)
    parser.add_argument(
        "--affine-directory",
        type=Path,
        required=True,
        help="directory containing p11_affine_<channel>_qrows_float_v2_20260827.txt",
    )
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    matrix, target = load_modular_export(args.modular_export)
    counts = np.asarray(load_theta_counts(args.theta_output), dtype=object)
    trace_harmonic = parse_named_vector(args.trace_output, "TRACE_HARMONIC_COEFS")
    limit = min(len(matrix), len(counts), len(trace_harmonic))
    trace_x, trace_residual = trace_coordinates(matrix, trace_harmonic, limit)
    trace_target = float(target @ trace_x)

    reductions = {}
    for channel in CHANNELS:
        path = (
            args.affine_directory
            / f"p11_affine_{channel}_qrows_float_v2_20260827.txt"
        )
        base, reduced_matrix, reduced_target_base, reduced_target = (
            load_affine_qrows(path)
        )
        if len(base) < limit:
            raise ValueError(f"not enough reduced coefficient rows in {path}")
        reductions[channel] = {
            "path": str(path),
            "base": base,
            "matrix": reduced_matrix,
            "target_base": reduced_target_base,
            "target": reduced_target,
        }

    cases = component_cases(args.p)
    solved = []
    expanded_intervals: list[tuple[float, float, int]] = []
    for case in cases:
        reduction = reductions[str(case["channel"])]
        result = solve_channel(
            reduction["matrix"],
            reduction["target"],
            reduction["base"],
            reduction["target_base"],
            counts,
            trace_harmonic,
            args.p,
            str(case["channel"]),
            int(case["representation_dimension"]),
        )
        merged = {**case, **result}
        solved.append(merged)
        if result["minimum"]["success"] and result["maximum"]["success"]:
            for _ in range(int(case["component_count"])):
                expanded_intervals.append(
                    (
                        float(result["minimum"]["value"]),
                        float(result["maximum"]["value"]),
                        int(case["representation_dimension"]),
                    )
                )

    p = args.p
    n = p * p + 1
    d = n // 2
    zdim = n * (n - 6) // 8
    lbar = 8 * (n - 2) / (n - 6)
    spherical = 8 * n / (n + 4)
    variance = maximum_weighted_variance(expanded_intervals, trace_target)
    scale = -zdim * (lbar - spherical) / trace_target
    variance["poisson_scale_from_trace"] = scale
    variance["Phi_frobenius_variance_max"] = (
        scale * scale * variance["normalized_variance_max"]
    )
    phi_endpoint_values = [
        spherical - scale * endpoint
        for lower, upper, _weight in expanded_intervals
        for endpoint in (lower, upper)
    ]
    variance["Phi_eigenvalue_lower_bound"] = min(phi_endpoint_values)
    variance["Phi_eigenvalue_upper_bound"] = max(phi_endpoint_values)
    variance["closes_principal_spectral_floor_numerically"] = (
        variance["Phi_eigenvalue_lower_bound"] >= 6.0
    )
    variance["R1_exact_threshold"] = n * (lbar - 6) ** 2 / 2
    variance["R1_strong_n_over_12_threshold"] = 2 * n
    variance["closes_exact_R1_numerically"] = (
        variance["Phi_frobenius_variance_max"]
        <= variance["R1_exact_threshold"]
    )
    variance["closes_strong_R1_numerically"] = (
        variance["Phi_frobenius_variance_max"]
        <= variance["R1_strong_n_over_12_threshold"]
    )

    output = {
        "experiment": "r1_p11_trace_coupled_lp",
        "status": "floating_reconnaissance_requires_exact_certificate",
        "p": p,
        "dimension": d,
        "Z_dimension": zdim,
        "component_count": len(expanded_intervals),
        "trace_modular_residual": trace_residual,
        "trace_target": trace_target,
        "affine_free_dimension": int(
            next(iter(reductions.values()))["matrix"].shape[1]
        ),
        "affine_reduction_files": {
            channel: reduction["path"]
            for channel, reduction in reductions.items()
        },
        "cases": solved,
        "variance_bound": variance,
    }
    payload = json.dumps(output, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
