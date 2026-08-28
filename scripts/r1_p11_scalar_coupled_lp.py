#!/usr/bin/env python3
"""Reconnaissance LP coupling p=11 harmonic and ordinary theta series.

This avoids enumerating the ordinary theta coefficients.  The scalar theta
series is required to lie in a Fourier-pivot-normalized modular subspace, to
have the proved coefficients through exponent 28, and to have nonnegative
coefficients.  For every later shell, its harmonic eigenvalue ``A_s`` is
coupled to the scalar shell count ``N_s`` through the universal rank-one
operator bounds

    0 <= q_s = A_s + 2 N_s r_s^2 / (d(d+2))
       <= N_s r_s^2 (d-1)/d.

The scalar modular subspace may be an enlargement of the true theta space,
so a closing bound is sound after exact certification while a failed bound
is only diagnostic.  SciPy/HiGHS output from this script is reconnaissance,
not a proof certificate.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


CHANNELS = ("circle-kernel", "circle-low", "circle-high")


def component_cases(p: int) -> list[dict[str, int | str]]:
    if p % 4 != 3:
        raise NotImplementedError("the current component table is the p=11 branch")
    return [
        {
            "name": "circle-kernel-principal",
            "channel": "circle-kernel",
            "representation_dimension": p * p + 1,
            "component_count": (p - 1) * (p - 3) // 8,
        },
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


def parse_exact_vector(line: str) -> list[Fraction]:
    start = line.find("[")
    end = line.rfind("]")
    if start < 0 or end < start:
        raise ValueError("no exact vector found")
    payload = line[start + 1 : end].strip()
    if not payload:
        return []
    return [Fraction(token.strip()) for token in payload.split(",")]


def load_rows(path: Path) -> list[list[Fraction]]:
    rows = [
        parse_exact_vector(line)
        for line in path.read_text().splitlines()
        if line.strip()
    ]
    if not rows or len({len(row) for row in rows}) != 1:
        raise ValueError(f"empty or inconsistent exact row file: {path}")
    return rows


def known_scalar_counts(p: int) -> dict[int, int]:
    if p != 11:
        raise NotImplementedError("the current affine reductions are p=11")
    return {
        0: 1,
        p: 2 * (p * p + 1),
        2 * (p - 1): p * (p + 1) * (p * p + 1),
        2 * (p + 1): p * p * (p * p + 1),
        3 * p - 6: p * p * (p - 1) * (p + 7) * (p * p + 1) // 6,
    }


def solve_channel(
    scalar_rows_exact: list[list[Fraction]],
    harmonic_rows_exact: list[list[Fraction]],
    p: int,
    channel: str,
    representation_dimension: int,
    case_name: str,
) -> dict[str, object]:
    if len(harmonic_rows_exact) < 2:
        raise ValueError("harmonic export needs coefficient rows and a target row")
    scalar_full = np.asarray(scalar_rows_exact, dtype=np.float64)
    harmonic_all = np.asarray(harmonic_rows_exact, dtype=np.float64)
    harmonic = harmonic_all[:-1]
    target_all = harmonic_all[-1]
    base = harmonic[:, 0]
    hrows = harmonic[:, 1:]
    target_base = float(target_all[0])
    target = target_all[1:]

    limit = min(len(scalar_full), len(harmonic))
    scalar_full = scalar_full[:limit]
    harmonic = harmonic[:limit]
    base = base[:limit]
    hrows = hrows[:limit]
    hdim = hrows.shape[1]
    d = (p * p + 1) // 2
    fixed_through = 2 * (p + 3)
    known = known_scalar_counts(p)

    # The exported scalar basis is pivot-normalized.  Eliminate all scalar
    # coordinates fixed by exponents 0..28 exactly at the matrix level.  This
    # removes a badly conditioned rank-15 equality block before calling HiGHS.
    fixed_coordinates: dict[int, float] = {}
    for exponent in range(min(fixed_through + 1, limit)):
        exact_row = scalar_rows_exact[exponent]
        nonzero = [(index, value) for index, value in enumerate(exact_row) if value]
        expected = Fraction(known.get(exponent, 0))
        if not nonzero:
            if expected:
                raise ArithmeticError(f"scalar fixed row {exponent} is identically zero")
            continue
        if len(nonzero) != 1:
            raise ArithmeticError(
                f"scalar basis is not low-pivot normalized at exponent {exponent}"
            )
        coordinate, coefficient = nonzero[0]
        value = expected / coefficient
        prior = fixed_coordinates.get(coordinate)
        if prior is not None and prior != float(value):
            raise ArithmeticError(f"inconsistent fixed scalar coordinate {coordinate}")
        fixed_coordinates[coordinate] = float(value)
    fixed_columns = sorted(fixed_coordinates)
    free_columns = [
        index for index in range(scalar_full.shape[1]) if index not in fixed_coordinates
    ]
    fixed_values = np.asarray([fixed_coordinates[index] for index in fixed_columns])
    scalar_base = scalar_full[:, fixed_columns] @ fixed_values
    scalar = scalar_full[:, free_columns]
    sdim = scalar.shape[1]
    nvars = hdim + sdim

    aub: list[np.ndarray] = []
    bub: list[float] = []
    names_ub: list[str] = []
    fixed_checks: list[dict[str, object]] = []

    def add_ub(name: str, row: np.ndarray, rhs: float) -> None:
        if np.max(np.abs(row), initial=0.0) == 0.0:
            if 0.0 > rhs + 1e-12:
                raise ArithmeticError(f"constant infeasible inequality {name}: 0<={rhs}")
            return
        scale = max(1.0, abs(rhs), float(np.max(np.abs(row))))
        aub.append(row / scale)
        bub.append(rhs / scale)
        names_ub.append(name)

    # Exact scalar coefficients already proved geometrically through s=28.
    for exponent in range(min(fixed_through + 1, limit)):
        expected = float(known.get(exponent, 0))
        reconstructed = scalar_base[exponent]
        free_scalar_max = float(np.max(np.abs(scalar[exponent]), initial=0.0))
        radius_sq = exponent / (2 * p)
        radial_per_vector = 2 * radius_sq**2 / (d * (d + 2))
        raw = float(base[exponent] + radial_per_vector * expected)
        free_max = float(np.max(np.abs(hrows[exponent]), initial=0.0))
        upper = (
            expected
            * radius_sq**2
            * (d - 1)
            / (d * representation_dimension)
        )
        if (
            free_scalar_max > 1e-12
            or abs(reconstructed - expected) > 1e-9
            or free_max > 1e-12
            or raw < -1e-9
            or raw > upper + 1e-9
        ):
            raise ArithmeticError(
                f"fixed shell {exponent} fails: scalar={reconstructed}, "
                f"scalar_free={free_scalar_max}, harmonic_free={free_max}, "
                f"raw={raw}, upper={upper}"
            )
        if expected or abs(raw) > 0:
            fixed_checks.append(
                {
                    "exponent": exponent,
                    "count": int(expected),
                    "raw_value": raw,
                    "raw_upper": upper,
                }
            )

    # Later scalar coefficients are shell counts.  Raw quartic positivity and
    # the trace/rank-one bound couple each harmonic coefficient to that count.
    for exponent in range(fixed_through + 1, limit):
        nrow = scalar[exponent]
        nbase = float(scalar_base[exponent])
        hrow = hrows[exponent]
        radius_sq = exponent / (2 * p)
        radial_per_vector = 2 * radius_sq**2 / (d * (d + 2))
        # PSD and Schur scalarity give m_c*q_{s,c} <= tr(R_s).  The
        # pointwise trace bound is ||b_x||^2 <= r_s^2(d-1)/d.
        raw_upper_per_vector = (
            radius_sq**2 * (d - 1) / (d * representation_dimension)
        )

        row = np.zeros(nvars)
        row[hdim:] = -nrow
        add_ub(f"N{exponent}_nonnegative", row, nbase)

        # q_s >= 0: -A_s-rho_s <= 0.
        row = np.zeros(nvars)
        row[:hdim] = -hrow
        row[hdim:] = -radial_per_vector * nrow
        add_ub(
            f"q{exponent}_nonnegative",
            row,
            float(base[exponent] + radial_per_vector * nbase),
        )

        # q_s <= N_s r_s^2 (d-1)/d.
        row = np.zeros(nvars)
        row[:hdim] = hrow
        row[hdim:] = (radial_per_vector - raw_upper_per_vector) * nrow
        add_ub(
            f"q{exponent}_rank_one_upper",
            row,
            float(-base[exponent] + (raw_upper_per_vector - radial_per_vector) * nbase),
        )

    A_ub = np.asarray(aub, dtype=np.float64).reshape((-1, nvars))
    b_ub = np.asarray(bub, dtype=np.float64)
    stacked = A_ub
    column_max = np.max(np.abs(stacked), axis=0)
    variable_scale = np.ones(nvars)
    nonzero = column_max > 0
    variable_scale[nonzero] = 1.0 / column_max[nonzero]
    A_ub_scaled = A_ub * variable_scale
    objective = np.zeros(nvars)
    objective[:hdim] = target
    objective_scaled = objective * variable_scale

    results: dict[str, object] = {}
    for sense, c in (("minimum", objective_scaled), ("maximum", -objective_scaled)):
        result = linprog(
            c,
            A_ub=A_ub_scaled,
            b_ub=b_ub,
            bounds=[(None, None)] * nvars,
            method="highs",
            options={
                "dual_feasibility_tolerance": 1e-9,
                "primal_feasibility_tolerance": 1e-9,
            },
        )
        item: dict[str, object] = {
            "success": bool(result.success),
            "status": int(result.status),
            "message": result.message,
        }
        if result.success:
            original = result.x * variable_scale
            value = target_base + float(objective @ original)
            ub_slack = b_ub - A_ub @ original
            item.update(
                {
                    "value": value,
                    "minimum_normalized_inequality_slack": float(np.min(ub_slack)),
                    "active_inequalities": int(np.count_nonzero(ub_slack <= 1e-7)),
                }
            )
        results[sense] = item

    return {
        "case": case_name,
        "channel": channel,
        "representation_dimension": representation_dimension,
        "coefficient_rows_used": limit,
        "harmonic_coordinate_dimension": hdim,
        "scalar_coordinate_dimension": sdim,
        "fixed_scalar_rank": len(fixed_columns),
        "fixed_scalar_coordinates": fixed_columns,
        "equality_rows_after_elimination": 0,
        "inequality_rows": len(A_ub),
        "fixed_checks": fixed_checks,
        "bounds": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scalar-qrows", type=Path, required=True)
    parser.add_argument("--affine-prefix", type=Path, required=True)
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    scalar_rows = load_rows(args.scalar_qrows)
    output = {
        "experiment": "r1_p11_scalar_coupled_lp",
        "status": "floating_reconnaissance_not_proof",
        "p": args.p,
        "scalar_qrows": str(args.scalar_qrows),
        "channels": [],
    }
    for case in component_cases(args.p):
        channel = str(case["channel"])
        path = Path(f"{args.affine_prefix}{channel}_qrows_exact_20260827.txt")
        output["channels"].append(
            solve_channel(
                scalar_rows,
                load_rows(path),
                args.p,
                channel,
                int(case["representation_dimension"]),
                str(case["name"]),
            )
        )

    payload = json.dumps(output, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
