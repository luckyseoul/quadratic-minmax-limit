#!/usr/bin/env python3
"""Best p=11 R1 bound from modularity plus shellwise quartic positivity.

This is a diagnostic, not an R1 proposition.  It combines the exact
coefficient matrix of the p=11 Kohnen space with the ordinary theta counts
of the Paley dual lattice.  For a unit admissible tensor W and a shell of N
vectors of squared norm r, its harmonic coefficient A_s(W) satisfies

    -2 N r^2 / (d(d+2)) <= A_s(W)
      <= N r^2 ((d-1)/d - 2/(d(d+2))).

The first four nonempty shells are inserted with their proved exact
operators.  A three-variable convex mixture covers the square-circle
kernel, low, and high channels.  The resulting LP therefore computes the
strongest bound available from *only* these inputs through the supplied
coefficient range; it does not assert that every LP point is a lattice
theta series.
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15631 import harmonic_min_shell_sum
from e1_gmin_m4_prop15634 import second_shadow_spectrum
from e1_gmin_m4_prop15635 import third_pair_harmonic_coefficient
from e1_gmin_m4_prop15640 import harmonic_spectrum


CHANNELS = ("circle-kernel", "circle-low", "circle-high")


def parse_gp_vector(text: str) -> list[float]:
    """Parse one PARI vector, including PARI's spaced scientific notation."""
    start = text.find("[")
    end = text.rfind("]")
    if start < 0 or end < start:
        raise ValueError("no GP vector found")
    payload = re.sub(r"(?<=\d)\s+E(?=[+-]\d)", "e", text[start : end + 1])
    values = ast.literal_eval(payload)
    if not isinstance(values, list):
        raise TypeError("expected a GP row vector")
    return [float(value) for value in values]


def load_modular_export(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = [parse_gp_vector(line) for line in path.read_text().splitlines() if line.strip()]
    if len(rows) < 2:
        raise ValueError("modular export needs coefficient rows plus target row")
    width = len(rows[-1])
    if any(len(row) != width for row in rows):
        raise ValueError("inconsistent modular-export row widths")
    return np.asarray(rows[:-1], dtype=np.float64), np.asarray(rows[-1], dtype=np.float64)


def load_theta_counts(path: Path) -> list[int]:
    text = path.read_text()
    match = re.search(r"COEFS\s*=\s*(\[[^\n]*\])", text)
    payload = match.group(1) if match else text[text.find("[") : text.rfind("]") + 1]
    values = ast.literal_eval(payload)
    if not isinstance(values, list):
        raise TypeError("expected an ordinary-theta coefficient vector")
    counts = [int(value) for value in values]
    if any(value < 0 for value in counts):
        raise ValueError("ordinary-theta counts must be nonnegative")
    return counts


def shell_interval(count: int, exponent: int, p: int, d: int) -> tuple[float, float]:
    radius_sq = exponent / (2 * p)
    radial = 2 * count * radius_sq**2 / (d * (d + 2))
    lower = -radial
    upper = count * radius_sq**2 * (d - 1) / d - radial
    return lower, upper


def exact_early_operators(p: int) -> dict[int, np.ndarray]:
    by_channel: dict[int, np.ndarray] = {}
    by_channel[p] = np.full(3, float(harmonic_min_shell_sum(p)))
    by_channel[2 * (p - 1)] = np.asarray(
        [
            float(next(row["eigenvalue"] for row in second_shadow_spectrum(p)
                       if row["channel"] == channel))
            for channel in CHANNELS
        ]
    )
    by_channel[2 * (p + 1)] = np.full(
        3, float(third_pair_harmonic_coefficient(p))
    )
    by_channel[3 * p - 6] = np.asarray(
        [
            float(next(row["eigenvalue"] for row in harmonic_spectrum(p)
                       if row["channel"] == channel))
            for channel in CHANNELS
        ]
    )
    return by_channel


def solve_bounds(
    matrix: np.ndarray,
    target: np.ndarray,
    counts: list[int],
    p: int,
) -> dict:
    d = (p * p + 1) // 2
    limit = min(len(counts), matrix.shape[0])
    matrix = matrix[:limit]
    counts = counts[:limit]
    early = exact_early_operators(p)

    # Variables are 66 modular coordinates followed by three nonnegative
    # channel weights.  Known shell coefficients equal the corresponding
    # convex mixture of the three proved channel eigenvalues.
    width = matrix.shape[1]
    nvars = width + 3
    aub: list[np.ndarray] = []
    bub: list[float] = []
    aeq: list[np.ndarray] = []
    beq: list[float] = []
    intervals: dict[int, tuple[float, float]] = {}

    for exponent, (row, count) in enumerate(zip(matrix, counts)):
        full = np.zeros(nvars, dtype=np.float64)
        full[:width] = row
        if exponent in early:
            full[width:] = -early[exponent]
            scale = max(1.0, float(np.max(np.abs(full))))
            aeq.append(full / scale)
            beq.append(0.0)
            continue
        lower, upper = shell_interval(count, exponent, p, d)
        intervals[exponent] = (lower, upper)
        if count == 0 or exponent == 0:
            scale = max(1.0, float(np.max(np.abs(full))))
            aeq.append(full / scale)
            beq.append(0.0)
            continue
        bound_scale = max(abs(lower), abs(upper))
        aub.append(full / bound_scale)
        bub.append(upper / bound_scale)
        aub.append(-full / bound_scale)
        bub.append(-lower / bound_scale)

    channel_sum = np.zeros(nvars, dtype=np.float64)
    channel_sum[width:] = 1.0
    aeq.append(channel_sum)
    beq.append(1.0)

    A_ub = np.asarray(aub)
    b_ub = np.asarray(bub)
    A_eq = np.asarray(aeq)
    b_eq = np.asarray(beq)

    # Scale each free modular coordinate so the largest normalized constraint
    # coefficient in that column is one.  This is only numerical conditioning;
    # all reported residuals are checked in the original coordinates below.
    stacked = np.vstack((A_ub, A_eq))
    column_max = np.max(np.abs(stacked), axis=0)
    variable_scale = np.ones(nvars)
    nz = column_max[:width] > 0
    variable_scale[:width][nz] = 1.0 / column_max[:width][nz]
    A_ub_scaled = A_ub * variable_scale
    A_eq_scaled = A_eq * variable_scale

    objective = np.zeros(nvars)
    objective[:width] = target
    scaled_objective = objective * variable_scale
    bounds = [(None, None)] * width + [(0.0, 1.0)] * 3

    results = {}
    for sense, c in (("minimum", scaled_objective), ("maximum", -scaled_objective)):
        result = linprog(
            c,
            A_ub=A_ub_scaled,
            b_ub=b_ub,
            A_eq=A_eq_scaled,
            b_eq=b_eq,
            bounds=bounds,
            method="highs",
            options={"dual_feasibility_tolerance": 1e-9,
                     "primal_feasibility_tolerance": 1e-9},
        )
        row: dict[str, object] = {
            "success": bool(result.success),
            "status": int(result.status),
            "message": result.message,
        }
        if result.success:
            original = result.x * variable_scale
            value = float(objective @ original)
            eq_residual = float(np.max(np.abs(A_eq @ original - b_eq)))
            ub_residual = float(np.max(A_ub @ original - b_ub))
            row.update(
                value=value,
                channel_weights={
                    channel: float(original[width + index])
                    for index, channel in enumerate(CHANNELS)
                },
                normalized_eq_residual=eq_residual,
                normalized_ub_residual=ub_residual,
            )
        results[sense] = row

    fixed_counts = {
        str(exponent): counts[exponent]
        for exponent in early
        if exponent < len(counts)
    }
    expected_counts = {
        str(p): 2 * (p * p + 1),
        str(2 * (p - 1)): p * (p + 1) * (p * p + 1),
        str(2 * (p + 1)): p * p * (p * p + 1),
        str(3 * p - 6): p * p * (p - 1) * (p + 7) * (p * p + 1) // 6,
    }
    return {
        "experiment": "r1_p11_shell_positivity_lp",
        "status": "diagnostic_not_R1_proof",
        "p": p,
        "dimension": d,
        "coefficient_rows_used": limit,
        "modular_coordinate_dimension": width,
        "ordinary_nonempty_shells": sum(count > 0 for count in counts),
        "fixed_early_counts": fixed_counts,
        "expected_early_counts": expected_counts,
        "early_counts_match": fixed_counts == expected_counts,
        "bounds": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--modular-export", type=Path, required=True)
    parser.add_argument("--theta-output", type=Path, required=True)
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    matrix, target = load_modular_export(args.modular_export)
    counts = load_theta_counts(args.theta_output)
    result = solve_bounds(matrix, target, counts, args.p)
    payload = json.dumps(result, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
