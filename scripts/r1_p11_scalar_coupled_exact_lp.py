#!/usr/bin/env python3
"""Exact LP for the p=11 scalar/harmonic theta coupling.

This is the rational counterpart of ``r1_p11_scalar_coupled_lp.py``.  It
eliminates the fifteen scalar coordinates fixed through exponent 28, writes
QSopt_ex LP files over Q, and independently verifies optimal primal/dual
certificates with ``Fraction`` arithmetic when ``--solve`` is requested.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from functools import lru_cache
from fractions import Fraction
from math import isqrt
from pathlib import Path
from typing import Mapping

import numpy as np

from r1_p11_scalar_coupled_lp import (
    CHANNELS,
    component_cases,
    known_scalar_counts,
    load_rows,
)
from r1_p11_trace_coupled_exact_lp import (
    Constraint,
    ExactModel,
    run_qsopt,
    verify_certificate,
    write_lp,
)


def scalar_affine_reduction(
    rows: list[list[Fraction]],
    p: int,
    fixed_through: int,
    scalar_fixed_counts: Mapping[int, Fraction] | None = None,
) -> tuple[list[Fraction], list[tuple[Fraction, ...]], tuple[int, ...]]:
    """Impose the proved scalar coefficients by exact affine elimination.

    The original scalar export was Fourier-pivot normalized, so each proved
    coefficient fixed one coordinate directly.  Additional cusp conditions
    change the basis and destroy that convenient shape.  RREF over
    ``Fraction`` keeps the reduction basis-independent.
    """
    if not rows or not rows[0]:
        raise ValueError("empty scalar coefficient matrix")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("inconsistent scalar coefficient width")
    proved_prefix = known_scalar_counts(p)
    known = {
        exponent: Fraction(proved_prefix.get(exponent, 0))
        for exponent in range(min(fixed_through + 1, len(rows)))
    }
    for exponent, count in (scalar_fixed_counts or {}).items():
        if not 0 <= exponent < len(rows):
            raise ValueError(
                f"fixed scalar exponent {exponent} is outside the exported rows"
            )
        value = Fraction(count)
        if value.denominator != 1 or value < 0:
            raise ValueError(
                f"fixed scalar count at exponent {exponent} is not a nonnegative integer"
            )
        if exponent in known and known[exponent] != value:
            raise ValueError(
                f"fixed scalar count at exponent {exponent} conflicts with the proved prefix"
            )
        known[exponent] = value
    augmented: list[list[Fraction]] = []
    for exponent in sorted(known):
        if exponent < len(rows):
            augmented.append(list(rows[exponent]) + [known[exponent]])

    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(width):
        selected = next(
            (
                row_index
                for row_index in range(pivot_row, len(augmented))
                if augmented[row_index][column]
            ),
            None,
        )
        if selected is None:
            continue
        augmented[pivot_row], augmented[selected] = (
            augmented[selected],
            augmented[pivot_row],
        )
        scale = augmented[pivot_row][column]
        augmented[pivot_row] = [value / scale for value in augmented[pivot_row]]
        for row_index, row in enumerate(augmented):
            if row_index == pivot_row or not row[column]:
                continue
            multiple = row[column]
            augmented[row_index] = [
                value - multiple * pivot_value
                for value, pivot_value in zip(row, augmented[pivot_row])
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(augmented):
            break

    for row in augmented:
        if not any(row[:-1]) and row[-1]:
            raise ArithmeticError("proved scalar coefficients are inconsistent")

    fixed_columns = tuple(pivot_columns)
    free_columns = tuple(index for index in range(width) if index not in fixed_columns)
    particular = [Fraction()] * width
    for row_index, column in enumerate(fixed_columns):
        particular[column] = augmented[row_index][-1]
    kernel: list[tuple[Fraction, ...]] = []
    for free_column in free_columns:
        vector = [Fraction()] * width
        vector[free_column] = Fraction(1)
        for row_index, pivot_column in enumerate(fixed_columns):
            vector[pivot_column] = -augmented[row_index][free_column]
        kernel.append(tuple(vector))

    bases: list[Fraction] = []
    matrix: list[tuple[Fraction, ...]] = []
    for row in rows:
        bases.append(sum((value * coordinate for value, coordinate in zip(row, particular)), Fraction()))
        matrix.append(
            tuple(
                sum((value * coordinate for value, coordinate in zip(row, direction)), Fraction())
                for direction in kernel
            )
        )
    return bases, matrix, fixed_columns


def build_model(
    scalar_rows: list[list[Fraction]],
    harmonic_rows: list[list[Fraction]],
    p: int,
    representation_dimension: int,
    scalar_fixed_counts: Mapping[int, Fraction] | None = None,
) -> ExactModel:
    if len(harmonic_rows) < 2:
        raise ValueError("harmonic export needs coefficient rows and one target row")
    harmonic_coefficients = harmonic_rows[:-1]
    target_row = harmonic_rows[-1]
    hdim = len(target_row) - 1
    if hdim < 1 or any(len(row) != hdim + 1 for row in harmonic_rows):
        raise ValueError("inconsistent harmonic affine row width")

    limit = min(len(scalar_rows), len(harmonic_coefficients))
    fixed_through = 2 * (p + 3)
    scalar_base, scalar_matrix, _fixed = scalar_affine_reduction(
        scalar_rows[:limit], p, fixed_through, scalar_fixed_counts
    )
    sdim = len(scalar_matrix[0])
    nvars = hdim + sdim
    n = p * p + 1
    d = n // 2
    known = known_scalar_counts(p)
    constraints: list[Constraint] = []
    fixed_checks: list[dict[str, str | int]] = []

    def row_join(
        harmonic: tuple[Fraction, ...] | None = None,
        scalar: tuple[Fraction, ...] | None = None,
        scalar_factor: Fraction = Fraction(1),
    ) -> tuple[Fraction, ...]:
        return tuple(harmonic or (Fraction(),) * hdim) + tuple(
            scalar_factor * value
            for value in (scalar or (Fraction(),) * sdim)
        )

    def add(name: str, row: tuple[Fraction, ...], sense: str, rhs: Fraction) -> None:
        if len(row) != nvars:
            raise ValueError(f"wrong constraint width for {name}")
        if all(value == 0 for value in row):
            valid = (
                (sense == "=" and rhs == 0)
                or (sense == ">=" and 0 >= rhs)
                or (sense == "<=" and 0 <= rhs)
            )
            if not valid:
                raise ArithmeticError(f"constant infeasible row {name}: 0 {sense} {rhs}")
            return
        constraints.append(Constraint(name, row, sense, rhs))

    for exponent in range(limit):
        hbase = harmonic_coefficients[exponent][0]
        hrow = tuple(harmonic_coefficients[exponent][1:])
        nbase = scalar_base[exponent]
        nrow = scalar_matrix[exponent]
        radius_sq = Fraction(exponent, 2 * p)
        gamma = Fraction(2, d * (d + 2)) * radius_sq**2
        beta = (
            Fraction(d - 1, d * representation_dimension) * radius_sq**2
        )

        if exponent <= fixed_through:
            expected = Fraction(known.get(exponent, 0))
            if nbase != expected or any(nrow) or any(hrow):
                raise ArithmeticError(f"fixed affine row fails at exponent {exponent}")
            raw = hbase + gamma * expected
            upper = beta * expected
            if not 0 <= raw <= upper:
                raise ArithmeticError(
                    f"fixed raw bound fails at exponent {exponent}: {raw} not in [0,{upper}]"
                )
            if expected or raw:
                fixed_checks.append(
                    {
                        "exponent": exponent,
                        "count": int(expected),
                        "raw_value": str(raw),
                        "raw_upper": str(upper),
                    }
                )
            continue

        # N_s = nbase + nrow*z >= 0.
        add(f"N{exponent}_lo", row_join(scalar=nrow), ">=", -nbase)

        # q_s = hbase+hrow*y+gamma*N_s >= 0.
        qrow = tuple(hrow) + tuple(gamma * value for value in nrow)
        add(f"q{exponent}_lo", qrow, ">=", -(hbase + gamma * nbase))

        # q_s <= beta*N_s.
        upper_row = tuple(hrow) + tuple((gamma - beta) * value for value in nrow)
        add(
            f"q{exponent}_hi",
            upper_row,
            "<=",
            -hbase + (beta - gamma) * nbase,
        )

    return ExactModel(
        target_base=target_row[0],
        target=tuple(target_row[1:]) + (Fraction(),) * sdim,
        constraints=tuple(constraints),
        fixed_checks=tuple(fixed_checks),
    )


@lru_cache(maxsize=None)
def parity_fourth_sum_table(
    p: int, max_exponent: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Exact minima of sum r_i^4 at fixed sum r_i^2 and common parity.

    For ``r=2p*x=(pI+C)z``, every coordinate of ``r`` has parity equal to
    the scaled norm exponent ``s`` and ``sum r_i^2=2ps``.  The dynamic
    program below is the exhaustive separable recurrence over the ``p^2+1``
    coordinates.  NumPy is used only for exact int64 min-plus updates; all
    attainable costs are at most ``(2*p*max_exponent)^2``.
    """
    if p <= 0 or p % 2 == 0 or max_exponent < 0:
        raise ValueError("parity fourth-sum table needs odd positive p")
    n = p * p + 1
    maximum_square_sum = 2 * p * max_exponent
    if maximum_square_sum**2 >= 2**62:
        raise OverflowError("parity fourth-sum DP exceeds the int64 proof range")
    infinity = np.int64(2**62)
    tables: list[tuple[int, ...]] = []
    for parity in (0, 1):
        levels = [
            (magnitude * magnitude, magnitude**4)
            for magnitude in range(parity, isqrt(maximum_square_sum) + 1, 2)
        ]
        current = np.full(maximum_square_sum + 1, infinity, dtype=np.int64)
        current[0] = 0
        for _coordinate in range(n):
            following = np.full_like(current, infinity)
            for square, fourth in levels:
                if square > maximum_square_sum:
                    continue
                np.minimum(
                    following[square:],
                    current[: maximum_square_sum + 1 - square] + fourth,
                    out=following[square:],
                )
            current = following
        tables.append(tuple(int(value) for value in current))
    return tables[0], tables[1]


def projected_rank_one_trace_upper(
    p: int,
    exponent: int,
    parity_tables: tuple[tuple[int, ...], tuple[int, ...]] | None = None,
) -> Fraction:
    """Pointwise upper bound for ``||projection_Z(xx^T)||^2`` on a shell."""
    n = p * p + 1
    radius_sq = Fraction(exponent, 2 * p)
    cauchy_upper = Fraction(n - 2, n) * radius_sq**2
    if parity_tables is None:
        return cauchy_upper
    scaled_square_sum = 2 * p * exponent
    minimum_fourth_sum = parity_tables[exponent % 2][scaled_square_sum]
    if minimum_fourth_sum >= 2**62:
        # No integer vector with the required common parity exists, so the
        # corresponding lattice shell is empty and its trace contribution is
        # identically zero.
        return Fraction()
    exact_upper = (
        Fraction(n, n - 2) * radius_sq**2
        - Fraction(4 * (n - 1), n - 2)
        * Fraction(minimum_fourth_sum, 16 * p**4)
    )
    if not 0 <= exact_upper <= cauchy_upper:
        raise ArithmeticError(
            f"invalid parity trace upper at exponent {exponent}: {exact_upper}"
        )
    return exact_upper


def build_scalar_trace_budget_model(
    scalar_rows: list[list[Fraction]],
    harmonic_rows_by_channel: dict[str, list[list[Fraction]]],
    p: int,
    distinguished_case_name: str,
    parity_fourth_moment: bool = False,
    scalar_half_target_rows: list[list[Fraction]] | None = None,
    scalar_half_target_first: Fraction | None = None,
    scalar_fixed_counts: Mapping[int, Fraction] | None = None,
) -> tuple[ExactModel, tuple[dict[str, int | str], ...]]:
    """Couple every constituent through one dimension-weighted trace budget.

    Constituents in the same representation case are interchangeable.  To
    optimize one distinguished target it is therefore exact for this convex
    relaxation to keep that constituent and replace every other same-case
    constituent by their average.  Each resulting representative retains an
    independent harmonic modular form; only the scalar theta series is shared.
    """
    cases = component_cases(p)
    matching = [
        index
        for index, case in enumerate(cases)
        if str(case["name"]) == distinguished_case_name
    ]
    if len(matching) != 1:
        raise ValueError(f"unknown distinguished case: {distinguished_case_name}")
    distinguished_index = matching[0]

    representatives: list[dict[str, int | str]] = []
    distinguished_representative = -1
    for case_index, case in enumerate(cases):
        count = int(case["component_count"])
        if case_index == distinguished_index:
            distinguished_representative = len(representatives)
            representatives.append(
                {**case, "role": "distinguished", "multiplicity": 1}
            )
            if count > 1:
                representatives.append(
                    {**case, "role": "average", "multiplicity": count - 1}
                )
        else:
            representatives.append(
                {**case, "role": "average", "multiplicity": count}
            )
    if distinguished_representative < 0:  # pragma: no cover - guarded above
        raise ArithmeticError("distinguished representative was not constructed")

    harmonic_coefficients: dict[str, list[list[Fraction]]] = {}
    harmonic_targets: dict[str, list[Fraction]] = {}
    harmonic_dimensions: set[int] = set()
    limits = [len(scalar_rows)]
    for channel in CHANNELS:
        rows = harmonic_rows_by_channel[channel]
        if len(rows) < 2:
            raise ValueError(f"harmonic export is too short for {channel}")
        harmonic_coefficients[channel] = rows[:-1]
        harmonic_targets[channel] = rows[-1]
        harmonic_dimensions.add(len(rows[-1]) - 1)
        limits.append(len(rows) - 1)
    if len(harmonic_dimensions) != 1:
        raise ValueError("harmonic channels have inconsistent affine widths")
    hdim = harmonic_dimensions.pop()
    if hdim < 1:
        raise ValueError("harmonic affine dimension is zero")
    limit = min(limits)

    fixed_through = 2 * (p + 3)
    if (scalar_half_target_rows is None) != (scalar_half_target_first is None):
        raise ValueError("scalar half-cusp target rows and value must be supplied together")
    if scalar_half_target_rows is not None and (
        len(scalar_half_target_rows) != 10
        or any(len(row) != len(scalar_rows[0]) for row in scalar_half_target_rows)
    ):
        raise ValueError("the p=11 scalar half-cusp target must be a 10 by scalar-width block")
    scalar_reduction_rows = scalar_rows[:limit] + (scalar_half_target_rows or [])
    scalar_all_base, scalar_all_matrix, _fixed = scalar_affine_reduction(
        scalar_reduction_rows, p, fixed_through, scalar_fixed_counts
    )
    scalar_base = scalar_all_base[:limit]
    scalar_matrix = scalar_all_matrix[:limit]
    scalar_target_base: Fraction | None = None
    scalar_target_matrix: tuple[Fraction, ...] | None = None
    if scalar_half_target_rows is not None:
        target_bases = scalar_all_base[limit:]
        target_matrices = scalar_all_matrix[limit:]
        phase_factors = (
            Fraction(1), Fraction(), Fraction(), Fraction(), Fraction(1),
            Fraction(1), Fraction(1), Fraction(), Fraction(1), Fraction(1, 2),
        )
        if any(
            base != factor * target_bases[0]
            or row != tuple(factor * value for value in target_matrices[0])
            for base, row, factor in zip(target_bases, target_matrices, phase_factors)
        ):
            raise ArithmeticError("unexpected scalar half-cusp phase line")
        scalar_target_base = target_bases[0]
        scalar_target_matrix = target_matrices[0]
    sdim = len(scalar_matrix[0])
    harmonic_variable_count = len(representatives) * hdim
    nvars = harmonic_variable_count + sdim
    n = p * p + 1
    d = n // 2
    zdim = n * (n - 6) // 8
    known = known_scalar_counts(p)
    constraints: list[Constraint] = []
    fixed_checks: list[dict[str, object]] = []
    parity_tables = (
        parity_fourth_sum_table(p, limit - 1) if parity_fourth_moment else None
    )

    represented_dimension = sum(
        int(row["multiplicity"]) * int(row["representation_dimension"])
        for row in representatives
    )
    if represented_dimension != zdim:
        raise ArithmeticError(
            f"representative dimensions sum to {represented_dimension}, not {zdim}"
        )

    def add(name: str, row: tuple[Fraction, ...], sense: str, rhs: Fraction) -> None:
        if len(row) != nvars:
            raise ValueError(f"wrong constraint width for {name}")
        if all(value == 0 for value in row):
            valid = (
                (sense == "=" and rhs == 0)
                or (sense == ">=" and 0 >= rhs)
                or (sense == "<=" and 0 <= rhs)
            )
            if not valid:
                raise ArithmeticError(f"constant infeasible row {name}: 0 {sense} {rhs}")
            return
        constraints.append(Constraint(name, row, sense, rhs))

    def scalar_block(values: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
        return (Fraction(),) * harmonic_variable_count + values

    if scalar_target_matrix is not None and scalar_target_base is not None:
        add(
            "scalar_half_target_eq",
            scalar_block(scalar_target_matrix),
            "=",
            scalar_half_target_first - scalar_target_base,
        )

    for exponent in range(limit):
        nbase = scalar_base[exponent]
        nrow = scalar_matrix[exponent]
        radius_sq = Fraction(exponent, 2 * p)
        gamma = Fraction(2, d * (d + 2)) * radius_sq**2
        trace_per_vector_upper = projected_rank_one_trace_upper(
            p, exponent, parity_tables
        )

        q_bases: list[Fraction] = []
        q_rows: list[tuple[Fraction, ...]] = []
        for representative_index, representative in enumerate(representatives):
            channel = str(representative["channel"])
            harmonic = harmonic_coefficients[channel][exponent]
            if len(harmonic) != hdim + 1:
                raise ValueError(f"inconsistent harmonic row at {channel}:{exponent}")
            q_bases.append(harmonic[0] + gamma * nbase)
            row = [Fraction()] * nvars
            offset = representative_index * hdim
            row[offset : offset + hdim] = harmonic[1:]
            row[harmonic_variable_count:] = [gamma * value for value in nrow]
            q_rows.append(tuple(row))

        trace_base = sum(
            (
                Fraction(
                    int(representative["multiplicity"])
                    * int(representative["representation_dimension"])
                )
                * q_base
                for representative, q_base in zip(representatives, q_bases)
            ),
            Fraction(),
        )
        trace_row = [Fraction()] * nvars
        for representative, qrow in zip(representatives, q_rows):
            weight = Fraction(
                int(representative["multiplicity"])
                * int(representative["representation_dimension"])
            )
            for index, value in enumerate(qrow):
                trace_row[index] += weight * value

        if exponent <= fixed_through:
            expected = Fraction(known.get(exponent, 0))
            if nbase != expected or any(nrow) or any(any(row) for row in q_rows):
                raise ArithmeticError(f"fixed affine row fails at exponent {exponent}")
            trace_upper = trace_per_vector_upper * expected
            if any(value < 0 for value in q_bases) or not 0 <= trace_base <= trace_upper:
                raise ArithmeticError(
                    f"fixed trace budget fails at exponent {exponent}: "
                    f"q={q_bases}, trace={trace_base}, upper={trace_upper}"
                )
            if expected or trace_base:
                fixed_checks.append(
                    {
                        "exponent": exponent,
                        "count": int(expected),
                        "raw_values": [str(value) for value in q_bases],
                        "raw_trace": str(trace_base),
                        "raw_trace_upper": str(trace_upper),
                    }
                )
            continue

        add(f"N{exponent}_lo", scalar_block(nrow), ">=", -nbase)
        for representative_index, (q_base, qrow) in enumerate(zip(q_bases, q_rows)):
            add(f"q{exponent}_r{representative_index}_lo", qrow, ">=", -q_base)

        # sum_c m_c q_{s,c} <= B_s N_s.  The scalar block already
        # contains zdim*gamma*nrow through the weighted q rows.
        for index, value in enumerate(nrow):
            trace_row[harmonic_variable_count + index] -= (
                trace_per_vector_upper * value
            )
        add(
            f"trace{exponent}_hi",
            tuple(trace_row),
            "<=",
            trace_per_vector_upper * nbase - trace_base,
        )

    distinguished = representatives[distinguished_representative]
    target_row = harmonic_targets[str(distinguished["channel"])]
    target = [Fraction()] * nvars
    offset = distinguished_representative * hdim
    target[offset : offset + hdim] = target_row[1:]
    return (
        ExactModel(
            target_base=target_row[0],
            target=tuple(target),
            constraints=tuple(constraints),
            fixed_checks=tuple(fixed_checks),
        ),
        tuple(representatives),
    )


def condition_model(model: ExactModel) -> ExactModel:
    """Exactly equilibrate rows and free variables without changing the LP.

    The q-row coordinates span many orders of magnitude.  In particular, a
    unit floor in the row norm leaves entire late-shell rows tiny and makes a
    floating-point bootstrap singular before an exact solver can take over.
    Alternating exact infinity-norm row and column scalings keeps every pass
    equivalent over Q while exposing a numerically usable matrix.
    """
    constraints = model.constraints
    target = model.target
    width = len(target)

    for _pass in range(4):
        row_scaled: list[Constraint] = []
        for constraint in constraints:
            magnitude = max(
                [abs(constraint.rhs)]
                + [abs(value) for value in constraint.coefficients]
            )
            if magnitude == 0:
                raise ArithmeticError(f"zero constraint during conditioning: {constraint.name}")
            row_scaled.append(
                Constraint(
                    constraint.name,
                    tuple(value / magnitude for value in constraint.coefficients),
                    constraint.sense,
                    constraint.rhs / magnitude,
                )
            )

        column_max = [
            max(
                (abs(row.coefficients[index]) for row in row_scaled),
                default=Fraction(),
            )
            for index in range(width)
        ]
        variable_scale = [
            Fraction(1) / value if value else Fraction(1) for value in column_max
        ]
        constraints = tuple(
            Constraint(
                row.name,
                tuple(
                    coefficient * scale
                    for coefficient, scale in zip(row.coefficients, variable_scale)
                ),
                row.sense,
                row.rhs,
            )
            for row in row_scaled
        )
        target = tuple(
            coefficient * scale
            for coefficient, scale in zip(target, variable_scale)
        )

    # Leave row norms at one after the final column pass.
    final_constraints: list[Constraint] = []
    for constraint in constraints:
        magnitude = max(
            [abs(constraint.rhs)]
            + [abs(value) for value in constraint.coefficients]
        )
        final_constraints.append(
            Constraint(
                constraint.name,
                tuple(value / magnitude for value in constraint.coefficients),
                constraint.sense,
                constraint.rhs / magnitude,
            )
        )

    return ExactModel(
        target_base=model.target_base,
        target=target,
        constraints=tuple(final_constraints),
        fixed_checks=model.fixed_checks,
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scalar-qrows", type=Path, required=True)
    parser.add_argument("--affine-prefix", type=Path, required=True)
    parser.add_argument("--channel", choices=CHANNELS, action="append")
    parser.add_argument("--case", action="append")
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--solve", action="store_true")
    parser.add_argument(
        "--raw-unconditioned",
        action="store_true",
        help="write the equivalent unscaled exact LP (mainly for regression)",
    )
    parser.add_argument(
        "--scalar-fixed-count",
        metavar="EXPONENT=COUNT",
        action="append",
        default=[],
        help="fix an additional exact ordinary-theta coefficient (repeatable)",
    )
    parser.add_argument(
        "--esolver", type=Path, default=Path("/home/nick/.local/qsopt-ex/usr/bin/esolver")
    )
    parser.add_argument(
        "--qsopt-library-directory",
        type=Path,
        default=Path("/home/nick/.local/qsopt-ex/usr/lib/x86_64-linux-gnu"),
    )
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)

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
    cases = component_cases(args.p)
    if args.channel:
        cases = [case for case in cases if str(case["channel"]) in args.channel]
    if args.case:
        unknown = set(args.case) - {str(case["name"]) for case in cases}
        if unknown:
            raise ValueError(f"unknown or channel-filtered cases: {sorted(unknown)}")
        cases = [case for case in cases if str(case["name"]) in args.case]
    report: dict[str, object] = {
        "experiment": "r1_p11_scalar_coupled_exact_lp",
        "status": "exact_qsopt_certified" if args.solve else "exact_lps_generated",
        "p": args.p,
        "scalar_fixed_counts": {
            str(exponent): str(count)
            for exponent, count in sorted(scalar_fixed_counts.items())
        },
        "channels": [],
    }
    for case in cases:
        channel = str(case["channel"])
        case_name = str(case["name"])
        representation_dimension = int(case["representation_dimension"])
        harmonic_path = Path(
            f"{args.affine_prefix}{channel}_qrows_exact_20260827.txt"
        )
        raw_model = build_model(
            scalar_rows,
            load_rows(harmonic_path),
            args.p,
            representation_dimension,
            scalar_fixed_counts=scalar_fixed_counts,
        )
        model = raw_model if args.raw_unconditioned else condition_model(raw_model)
        channel_row: dict[str, object] = {
            "case": case_name,
            "channel": channel,
            "representation_dimension": representation_dimension,
            "variable_count": len(model.target),
            "constraint_count": len(model.constraints),
            "exact_conditioning": not args.raw_unconditioned,
            "fixed_checks": list(model.fixed_checks),
            "bounds": {},
        }
        for sense in ("minimum", "maximum"):
            lp_path = args.output_directory / f"{case_name}_{sense}.lp"
            solution_path = args.output_directory / f"{case_name}_{sense}.sol"
            objective = write_lp(lp_path, model, sense)
            item: dict[str, object] = {
                "lp": str(lp_path),
                "lp_sha256": sha256(lp_path),
            }
            if args.solve:
                log = run_qsopt(
                    args.esolver,
                    args.qsopt_library_directory,
                    lp_path,
                    solution_path,
                )
                certificate = verify_certificate(model, objective, solution_path)
                signed = Fraction(str(certificate["solver_objective"]))
                variable_value = signed if sense == "minimum" else -signed
                item.update(
                    {
                        "value": str(model.target_base + variable_value),
                        "certificate": certificate,
                        "solution": str(solution_path),
                        "solution_sha256": sha256(solution_path),
                        "solver_log_tail": log.splitlines()[-8:],
                    }
                )
            channel_row["bounds"][sense] = item
        report["channels"].append(channel_row)

    report_path = args.output_directory / "report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
