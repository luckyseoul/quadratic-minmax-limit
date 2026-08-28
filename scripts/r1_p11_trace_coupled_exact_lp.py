#!/usr/bin/env python3
"""Build and verify exact QSopt_ex certificates for the p=11 R1 channel LP.

The input affine rows are exact rational q-coordinate reductions.  For every
dual shell this imposes positivity of the raw quartic channel eigenvalue and
the conserved-trace upper bound.  QSopt_ex minimizes over Q; its primal and
dual certificates are then checked independently with ``Fraction`` arithmetic.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import itertools
import json
import os
import re
import subprocess
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable

from r1_p11_scalar_coupled_lp import CHANNELS, component_cases


@dataclass(frozen=True)
class Constraint:
    name: str
    coefficients: tuple[Fraction, ...]
    sense: str
    rhs: Fraction


@dataclass(frozen=True)
class ExactModel:
    target_base: Fraction
    target: tuple[Fraction, ...]
    constraints: tuple[Constraint, ...]
    fixed_checks: tuple[dict[str, str | int], ...]


def parse_exact_vector(line: str) -> list[Fraction]:
    start = line.find("[")
    end = line.rfind("]")
    if start < 0 or end < start:
        raise ValueError("no exact GP vector found")
    payload = line[start + 1 : end].strip()
    if not payload:
        return []
    return [Fraction(token.strip()) for token in payload.split(",")]


def parse_named_vector(path: Path, name: str) -> list[Fraction]:
    match = re.search(
        rf"{re.escape(name)}\s*=\s*(\[[^\n]*\])", path.read_text()
    )
    if match is None:
        raise ValueError(f"{name} not found in {path}")
    return parse_exact_vector(match.group(1))


def load_counts(path: Path) -> list[int]:
    text = path.read_text()
    match = re.search(r"COEFS\s*=\s*(\[[^\n]*\])", text)
    if match is None:
        raise ValueError(f"COEFS not found in {path}")
    values = ast.literal_eval(match.group(1))
    if not isinstance(values, list) or any(int(x) < 0 for x in values):
        raise ValueError("invalid ordinary-theta counts")
    return [int(x) for x in values]


def load_qrows(
    path: Path,
) -> tuple[list[Fraction], list[tuple[Fraction, ...]], Fraction, tuple[Fraction, ...]]:
    rows = [parse_exact_vector(line) for line in path.read_text().splitlines() if line]
    if len(rows) < 2:
        raise ValueError(f"q-row file is too short: {path}")
    width = len(rows[-1])
    if width < 2 or any(len(row) != width for row in rows):
        raise ValueError(f"inconsistent q-row width: {path}")
    coefficient_rows = rows[:-1]
    target_row = rows[-1]
    return (
        [row[0] for row in coefficient_rows],
        [tuple(row[1:]) for row in coefficient_rows],
        target_row[0],
        tuple(target_row[1:]),
    )


def load_trace_target(path: Path, trace: list[Fraction]) -> Fraction:
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    if len(lines) != 2:
        raise ValueError("trace-target export must contain exactly two vectors")
    pivots = [int(x) for x in parse_exact_vector(lines[0])]
    target = parse_exact_vector(lines[1])
    if len(pivots) != len(target):
        raise ValueError("trace-target pivot/functional width mismatch")
    if max(pivots) >= len(trace):
        raise ValueError("trace series does not reach every target pivot")
    return sum((coefficient * trace[pivot] for pivot, coefficient in zip(pivots, target)), Fraction())


def build_model(
    base: list[Fraction],
    matrix: list[tuple[Fraction, ...]],
    target_base: Fraction,
    target: tuple[Fraction, ...],
    counts: list[int],
    trace: list[Fraction],
    p: int,
    representation_dimension: int,
) -> ExactModel:
    n = p * p + 1
    d = n // 2
    zdim = n * (n - 6) // 8
    fixed_through = 2 * (p + 3)
    limit = min(len(base), len(matrix), len(counts), len(trace))
    constraints: list[Constraint] = []
    fixed_checks: list[dict[str, str | int]] = []

    def add_constraint(
        name: str, row: tuple[Fraction, ...], sense: str, rhs: Fraction
    ) -> None:
        if all(value == 0 for value in row):
            valid = (sense == "=" and rhs == 0) or (sense == ">=" and 0 >= rhs) or (sense == "<=" and 0 <= rhs)
            if not valid:
                raise ArithmeticError(f"constant infeasible constraint {name}: 0 {sense} {rhs}")
            return
        constraints.append(Constraint(name, row, sense, rhs))

    for exponent in range(limit):
        row = matrix[exponent]
        count = counts[exponent]
        radius_sq = Fraction(exponent, 2 * p)
        radial = Fraction(2 * count, d * (d + 2)) * radius_sq**2
        tau = trace[exponent] + zdim * radial
        if tau < 0:
            raise ArithmeticError(f"negative raw trace mass at exponent {exponent}: {tau}")

        if exponent <= fixed_through:
            if any(row):
                raise ArithmeticError(f"nonzero free row at fixed exponent {exponent}")
            raw = base[exponent] + radial
            upper = tau / representation_dimension
            if not (0 <= raw <= upper):
                raise ArithmeticError(
                    f"fixed raw bound fails at exponent {exponent}: {raw} not in [0,{upper}]"
                )
            if count or tau:
                fixed_checks.append(
                    {
                        "exponent": exponent,
                        "raw_value": str(raw),
                        "raw_upper": str(upper),
                    }
                )
            continue

        if count == 0 or tau == 0:
            harmonic = -radial if count else Fraction()
            add_constraint(
                f"s{exponent}_eq", row, "=", harmonic - base[exponent]
            )
            continue

        shift = base[exponent] + radial
        add_constraint(f"s{exponent}_lo", row, ">=", -shift)
        add_constraint(
            f"s{exponent}_hi",
            row,
            "<=",
            tau / representation_dimension - shift,
        )

    return ExactModel(
        target_base=target_base,
        target=target,
        constraints=tuple(constraints),
        fixed_checks=tuple(fixed_checks),
    )


def build_shellwise_conserved_model(
    reductions: dict[
        str,
        tuple[
            list[Fraction],
            list[tuple[Fraction, ...]],
            Fraction,
            tuple[Fraction, ...],
        ],
    ],
    cases: list[dict[str, int | str]],
    distinguished_index: int,
    counts: list[int],
    trace: list[Fraction],
    p: int,
) -> tuple[ExactModel, tuple[dict[str, int | str], ...]]:
    """Couple all constituents by exact raw-mass conservation on each shell.

    Constituents of one case are interchangeable.  For a target interval it
    is therefore exact to retain the distinguished constituent and replace
    every other same-case constituent by their average.  Positivity is convex,
    so this symmetry quotient loses no feasible target endpoint.
    """
    if not 0 <= distinguished_index < len(cases):
        raise IndexError(distinguished_index)
    n = p * p + 1
    d = n // 2
    zdim = n * (n - 6) // 8
    fixed_through = 2 * (p + 3)

    widths = {len(reduction[3]) for reduction in reductions.values()}
    if len(widths) != 1:
        raise ValueError("all affine reductions must have one common width")
    block_width = widths.pop()
    representatives: list[dict[str, int | str]] = []
    for case_index, case in enumerate(cases):
        count = int(case["component_count"])
        if case_index == distinguished_index:
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
    total_width = block_width * len(representatives)
    distinguished_rep = next(
        index
        for index, representative in enumerate(representatives)
        if representative["role"] == "distinguished"
    )
    distinguished_channel = str(representatives[distinguished_rep]["channel"])
    _base, _matrix, target_base, target = reductions[distinguished_channel]
    objective = [Fraction()] * total_width
    objective[
        distinguished_rep * block_width : (distinguished_rep + 1) * block_width
    ] = target

    constraints: list[Constraint] = []
    fixed_checks: list[dict[str, str | int]] = []

    def extended_row(rep_index: int, row: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
        output = [Fraction()] * total_width
        output[rep_index * block_width : (rep_index + 1) * block_width] = row
        return tuple(output)

    def add_constraint(
        name: str, row: tuple[Fraction, ...], sense: str, rhs: Fraction
    ) -> None:
        if all(value == 0 for value in row):
            valid = (
                (sense == "=" and rhs == 0)
                or (sense == ">=" and 0 >= rhs)
                or (sense == "<=" and 0 <= rhs)
            )
            if not valid:
                raise ArithmeticError(
                    f"constant infeasible constraint {name}: 0 {sense} {rhs}"
                )
            return
        constraints.append(Constraint(name, row, sense, rhs))

    row_limits = [len(reductions[str(rep["channel"])][0]) for rep in representatives]
    limit = min(len(counts), len(trace), *row_limits)
    for exponent in range(limit):
        count = counts[exponent]
        radius_sq = Fraction(exponent, 2 * p)
        radial = Fraction(2 * count, d * (d + 2)) * radius_sq**2
        tau = trace[exponent] + zdim * radial
        if tau < 0:
            raise ArithmeticError(
                f"negative raw trace mass at exponent {exponent}: {tau}"
            )

        conservation = [Fraction()] * total_width
        conservation_constant = Fraction()
        raw_constants: list[Fraction] = []
        for rep_index, representative in enumerate(representatives):
            channel = str(representative["channel"])
            base, matrix, _target_base, _target = reductions[channel]
            row = matrix[exponent]
            raw_constant = base[exponent] + radial
            raw_constants.append(raw_constant)
            name = (
                f"s{exponent}_{representative['name']}_{representative['role']}"
            ).replace("-", "_")
            add_constraint(
                f"{name}_positive",
                extended_row(rep_index, row),
                ">=",
                -raw_constant,
            )
            weight = int(representative["multiplicity"]) * int(
                representative["representation_dimension"]
            )
            conservation_constant += weight * raw_constant
            offset = rep_index * block_width
            for column, coefficient in enumerate(row):
                conservation[offset + column] += weight * coefficient

        add_constraint(
            f"s{exponent}_mass",
            tuple(conservation),
            "=",
            tau - conservation_constant,
        )
        if exponent <= fixed_through:
            if any(
                any(reductions[str(rep["channel"])][1][exponent])
                for rep in representatives
            ):
                raise ArithmeticError(
                    f"nonzero free row at fixed exponent {exponent}"
                )
            if any(value < 0 for value in raw_constants):
                raise ArithmeticError(
                    f"negative fixed raw value at exponent {exponent}"
                )
            fixed_checks.append(
                {
                    "exponent": exponent,
                    "raw_mass": str(conservation_constant),
                    "trace_mass": str(tau),
                }
            )

    return (
        ExactModel(
            target_base=target_base,
            target=tuple(objective),
            constraints=tuple(constraints),
            fixed_checks=tuple(fixed_checks),
        ),
        tuple(representatives),
    )


def build_broad_channel_conserved_model(
    reductions: dict[
        str,
        tuple[
            list[Fraction],
            list[tuple[Fraction, ...]],
            Fraction,
            tuple[Fraction, ...],
        ],
    ],
    cases: list[dict[str, int | str]],
    distinguished_index: int,
    counts: list[int],
    broad_masses: dict[str, list[Fraction]],
    broad_targets: dict[str, Fraction],
    p: int,
) -> tuple[ExactModel, tuple[dict[str, int | str], ...]]:
    """Couple constituents by exact mass conservation in every broad channel.

    The square-circle operator splits ``Z`` into the three channels in
    ``CHANNELS``.  ``broad_masses[channel][s]`` is the exact trace of the raw
    positive shell operator on that whole eigenspace, while
    ``broad_targets[channel]`` is its dimension-normalized harmonic half-cusp
    trace.  Thus both the shell masses and transformed targets are conserved
    separately in each channel.

    As in :func:`build_shellwise_conserved_model`, constituents of one case
    are interchangeable.  Keeping the distinguished constituent and the
    average of its peers is an exact convex symmetry quotient for either
    endpoint of its target functional.
    """
    if not 0 <= distinguished_index < len(cases):
        raise IndexError(distinguished_index)
    if set(reductions) != set(CHANNELS):
        raise ValueError("affine reductions must contain exactly the broad channels")
    if set(broad_masses) != set(CHANNELS):
        raise ValueError("broad masses must contain exactly the broad channels")
    if set(broad_targets) != set(CHANNELS):
        raise ValueError("broad targets must contain exactly the broad channels")

    n = p * p + 1
    d = n // 2
    zdim = n * (n - 6) // 8
    fixed_through = 2 * (p + 3)
    widths = {len(reduction[3]) for reduction in reductions.values()}
    if len(widths) != 1:
        raise ValueError("all affine reductions must have one common width")
    block_width = widths.pop()
    if block_width < 1:
        raise ValueError("the common affine reduction has zero width")
    if any(
        len(row) != block_width
        for _base, matrix, _target_base, _target in reductions.values()
        for row in matrix
    ):
        raise ValueError("an affine coefficient row has the wrong width")

    representatives: list[dict[str, int | str]] = []
    for case_index, case in enumerate(cases):
        count = int(case["component_count"])
        if count < 1:
            raise ValueError(f"case {case['name']} has no constituents")
        if case_index == distinguished_index:
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

    total_width = block_width * len(representatives)
    distinguished_rep = next(
        index
        for index, representative in enumerate(representatives)
        if representative["role"] == "distinguished"
    )
    distinguished_channel = str(representatives[distinguished_rep]["channel"])
    _base, _matrix, target_base, target = reductions[distinguished_channel]
    objective = [Fraction()] * total_width
    objective[
        distinguished_rep * block_width : (distinguished_rep + 1) * block_width
    ] = target

    channel_dimensions = {
        channel: sum(
            int(case["component_count"])
            * int(case["representation_dimension"])
            for case in cases
            if str(case["channel"]) == channel
        )
        for channel in CHANNELS
    }
    if sum(channel_dimensions.values()) != zdim:
        raise ArithmeticError(
            f"broad channel dimensions sum to {sum(channel_dimensions.values())}, "
            f"not {zdim}"
        )
    represented_dimensions = {
        channel: sum(
            int(rep["multiplicity"]) * int(rep["representation_dimension"])
            for rep in representatives
            if str(rep["channel"]) == channel
        )
        for channel in CHANNELS
    }
    if represented_dimensions != channel_dimensions:
        raise ArithmeticError("the symmetry quotient changed a broad dimension")

    constraints: list[Constraint] = []
    fixed_checks: list[dict[str, str | int]] = []

    def extended_row(rep_index: int, row: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
        output = [Fraction()] * total_width
        output[rep_index * block_width : (rep_index + 1) * block_width] = row
        return tuple(output)

    def add_constraint(
        name: str, row: tuple[Fraction, ...], sense: str, rhs: Fraction
    ) -> None:
        if len(row) != total_width:
            raise ValueError(f"wrong constraint width for {name}")
        if all(value == 0 for value in row):
            valid = (
                (sense == "=" and rhs == 0)
                or (sense == ">=" and 0 >= rhs)
                or (sense == "<=" and 0 <= rhs)
            )
            if not valid:
                raise ArithmeticError(
                    f"constant infeasible constraint {name}: 0 {sense} {rhs}"
                )
            return
        constraints.append(Constraint(name, row, sense, rhs))

    row_limits = [len(reductions[str(rep["channel"])][0]) for rep in representatives]
    limit = min(
        len(counts),
        *(len(broad_masses[channel]) for channel in CHANNELS),
        *row_limits,
    )
    if limit < 1:
        raise ValueError("no common broad-channel coefficient prefix")
    if any(
        mass < 0
        for channel in CHANNELS
        for mass in broad_masses[channel][:limit]
    ):
        raise ArithmeticError("a supplied broad raw mass is negative")

    for exponent in range(limit):
        count = counts[exponent]
        radius_sq = Fraction(exponent, 2 * p)
        radial = Fraction(2 * count, d * (d + 2)) * radius_sq**2
        channel_rows = {
            channel: [Fraction()] * total_width for channel in CHANNELS
        }
        channel_constants = {channel: Fraction() for channel in CHANNELS}
        raw_constants: list[Fraction] = []

        for rep_index, representative in enumerate(representatives):
            channel = str(representative["channel"])
            base, matrix, _target_base, _target = reductions[channel]
            row = matrix[exponent]
            raw_constant = base[exponent] + radial
            raw_constants.append(raw_constant)
            name = (
                f"s{exponent}_{representative['name']}_{representative['role']}"
            ).replace("-", "_")
            add_constraint(
                f"{name}_positive",
                extended_row(rep_index, row),
                ">=",
                -raw_constant,
            )
            weight = int(representative["multiplicity"]) * int(
                representative["representation_dimension"]
            )
            channel_constants[channel] += weight * raw_constant
            offset = rep_index * block_width
            for column, coefficient in enumerate(row):
                channel_rows[channel][offset + column] += weight * coefficient

        for channel in CHANNELS:
            add_constraint(
                f"s{exponent}_{channel.replace('-', '_')}_mass",
                tuple(channel_rows[channel]),
                "=",
                broad_masses[channel][exponent] - channel_constants[channel],
            )

        if exponent <= fixed_through:
            if any(
                any(reductions[str(rep["channel"])][1][exponent])
                for rep in representatives
            ):
                raise ArithmeticError(
                    f"nonzero free row at fixed exponent {exponent}"
                )
            if any(value < 0 for value in raw_constants):
                raise ArithmeticError(
                    f"negative fixed raw value at exponent {exponent}"
                )
            for channel in CHANNELS:
                if channel_constants[channel] != broad_masses[channel][exponent]:
                    raise ArithmeticError(
                        f"fixed {channel} mass fails at exponent {exponent}"
                    )
            fixed_checks.append(
                {
                    "exponent": exponent,
                    **{
                        f"{channel}_mass": str(broad_masses[channel][exponent])
                        for channel in CHANNELS
                    },
                }
            )

    # This transformed-target identity is implied once the coefficient rows
    # span the full affine space, but retaining it makes the conservation used
    # in the endpoint proof explicit and independently checkable.
    for channel in CHANNELS:
        row = [Fraction()] * total_width
        constant = Fraction()
        for rep_index, representative in enumerate(representatives):
            if str(representative["channel"]) != channel:
                continue
            _base, _matrix, rep_target_base, rep_target = reductions[channel]
            weight = int(representative["multiplicity"]) * int(
                representative["representation_dimension"]
            )
            constant += weight * rep_target_base
            offset = rep_index * block_width
            for column, coefficient in enumerate(rep_target):
                row[offset + column] += weight * coefficient
        add_constraint(
            f"target_{channel.replace('-', '_')}_mass",
            tuple(row),
            "=",
            channel_dimensions[channel] * broad_targets[channel] - constant,
        )

    return (
        ExactModel(
            target_base=target_base,
            target=tuple(objective),
            constraints=tuple(constraints),
            fixed_checks=tuple(fixed_checks),
        ),
        tuple(representatives),
    )


def lp_number(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def linear_expression(coefficients: Iterable[Fraction]) -> str:
    terms: list[str] = []
    for index, coefficient in enumerate(coefficients, start=1):
        if coefficient == 0:
            continue
        magnitude = lp_number(abs(coefficient))
        term = f"{magnitude} y{index}"
        if not terms:
            terms.append(f"- {term}" if coefficient < 0 else term)
        else:
            terms.append(f" {'-' if coefficient < 0 else '+'} {term}")
    return "".join(terms) if terms else "0"


def write_lp(path: Path, model: ExactModel, sense: str) -> tuple[Fraction, ...]:
    if sense not in {"minimum", "maximum"}:
        raise ValueError(sense)
    objective = model.target if sense == "minimum" else tuple(-x for x in model.target)
    objective_expression = linear_expression(objective)
    if objective_expression == "0" and objective:
        # QSopt_ex's LP reader rejects a bare zero objective as a
        # coefficient without a variable.  An explicit zero coefficient is
        # the same objective and keeps feasibility-only models parseable.
        objective_expression = "0 y1"
    lines = ["Minimize", f" obj: {objective_expression}", "Subject To"]
    for constraint in model.constraints:
        lines.append(
            f" {constraint.name}: {linear_expression(constraint.coefficients)} "
            f"{constraint.sense} {lp_number(constraint.rhs)}"
        )
    lines.append("Bounds")
    lines.extend(f" y{index} free" for index in range(1, len(model.target) + 1))
    lines.append("End")
    path.write_text("\n".join(lines) + "\n")
    return objective


def parse_solution(path: Path) -> tuple[Fraction, dict[str, Fraction], dict[str, Fraction]]:
    text = path.read_text()
    if "status OPTIMAL" not in text:
        raise ArithmeticError(f"QSopt_ex did not report OPTIMAL in {path}")
    value_match = re.search(r"\bValue\s*=\s*([^\s]+)", text)
    if value_match is None:
        raise ValueError(f"objective value missing from {path}")
    value = Fraction(value_match.group(1))

    def section(start: str, end: str) -> dict[str, Fraction]:
        match = re.search(
            rf"{re.escape(start)}\n(.*?){re.escape(end)}", text, re.DOTALL
        )
        if match is None:
            return {}
        output: dict[str, Fraction] = {}
        for line in match.group(1).splitlines():
            item = re.match(r"\s*([^\s=]+)\s*=\s*([^\s]+)\s*$", line)
            if item:
                output[item.group(1)] = Fraction(item.group(2))
        return output

    return value, section("VARS:", "REDUCED COST:"), section("PI:", "SLACK:")


def verify_certificate(
    model: ExactModel,
    objective: tuple[Fraction, ...],
    solution_path: Path,
) -> dict[str, str | int]:
    solver_value, variables_by_name, pi_by_name = parse_solution(solution_path)
    variables = tuple(
        variables_by_name.get(f"y{index}", Fraction())
        for index in range(1, len(objective) + 1)
    )
    primal_value = sum((c * x for c, x in zip(objective, variables)), Fraction())
    if primal_value != solver_value:
        raise ArithmeticError(f"primal objective mismatch: {primal_value} != {solver_value}")

    for constraint in model.constraints:
        lhs = sum(
            (coefficient * variable for coefficient, variable in zip(constraint.coefficients, variables)),
            Fraction(),
        )
        valid = (
            (constraint.sense == "=" and lhs == constraint.rhs)
            or (constraint.sense == ">=" and lhs >= constraint.rhs)
            or (constraint.sense == "<=" and lhs <= constraint.rhs)
        )
        if not valid:
            raise ArithmeticError(
                f"primal constraint {constraint.name} fails: {lhs} {constraint.sense} {constraint.rhs}"
            )

    pi = {
        constraint.name: pi_by_name.get(constraint.name, Fraction())
        for constraint in model.constraints
    }
    for constraint in model.constraints:
        dual = pi[constraint.name]
        if constraint.sense == ">=" and dual < 0:
            raise ArithmeticError(f"negative dual on >= row {constraint.name}")
        if constraint.sense == "<=" and dual > 0:
            raise ArithmeticError(f"positive dual on <= row {constraint.name}")
    for index, expected in enumerate(objective):
        actual = sum(
            (
                pi[constraint.name] * constraint.coefficients[index]
                for constraint in model.constraints
            ),
            Fraction(),
        )
        if actual != expected:
            raise ArithmeticError(
                f"dual stationarity fails at y{index + 1}: {actual} != {expected}"
            )
    dual_value = sum(
        (pi[constraint.name] * constraint.rhs for constraint in model.constraints),
        Fraction(),
    )
    if dual_value != solver_value:
        raise ArithmeticError(f"dual objective mismatch: {dual_value} != {solver_value}")
    return {
        "solver_objective": str(solver_value),
        "nonzero_primal_variables": sum(value != 0 for value in variables),
        "nonzero_dual_variables": sum(value != 0 for value in pi.values()),
        "primal_constraints_verified": len(model.constraints),
        "dual_stationarity_equations_verified": len(objective),
    }


def run_qsopt(
    esolver: Path, library_directory: Path, lp_path: Path, solution_path: Path
) -> str:
    environment = os.environ.copy()
    previous = environment.get("LD_LIBRARY_PATH")
    environment["LD_LIBRARY_PATH"] = (
        str(library_directory)
        if not previous
        else f"{library_directory}:{previous}"
    )
    result = subprocess.run(
        [str(esolver), "-L", "-O", str(solution_path), str(lp_path)],
        text=True,
        capture_output=True,
        env=environment,
        check=False,
    )
    log = result.stdout + result.stderr
    if result.returncode != 0 or "Problem Solved Exactly" not in log:
        raise RuntimeError(f"QSopt_ex failed for {lp_path}:\n{log}")
    return log


def maximum_weighted_variance_exact(
    intervals: list[tuple[Fraction, Fraction, int]], total: Fraction
) -> dict[str, object]:
    weights = [Fraction(weight) for _lower, _upper, weight in intervals]
    lower = [row[0] for row in intervals]
    upper = [row[1] for row in intervals]
    total_weight = sum(weights, Fraction())
    if not (
        sum((w * x for w, x in zip(weights, lower)), Fraction())
        <= total
        <= sum((w * x for w, x in zip(weights, upper)), Fraction())
    ):
        raise ArithmeticError("trace target is outside the exact interval box")
    mean = total / total_weight
    # Ten of the fifteen p=11 constituents have the same interval and weight.
    # Enumerate upper-endpoint counts within identical groups, rather than all
    # 2^14 labelled endpoint assignments.  This is the same complete vertex
    # enumeration, reduced by permutation symmetry.
    grouped: dict[tuple[Fraction, Fraction, int], int] = {}
    for interval in intervals:
        grouped[interval] = grouped.get(interval, 0) + 1
    groups = [(*key, multiplicity) for key, multiplicity in grouped.items()]
    best: Fraction | None = None
    best_point: list[Fraction] | None = None
    vertices_checked = 0
    for pivot_group, (pivot_lower, pivot_upper, pivot_weight_int, pivot_count) in enumerate(groups):
        endpoint_counts = []
        for group_index, (_lo, _hi, _weight, multiplicity) in enumerate(groups):
            available = multiplicity - (1 if group_index == pivot_group else 0)
            endpoint_counts.append(range(available + 1))
        for upper_counts in itertools.product(*endpoint_counts):
            subtotal = Fraction()
            fixed_value = Fraction()
            point: list[Fraction] = []
            for group_index, (lo, hi, weight_int, multiplicity) in enumerate(groups):
                available = multiplicity - (1 if group_index == pivot_group else 0)
                high_count = upper_counts[group_index]
                low_count = available - high_count
                weight = Fraction(weight_int)
                subtotal += weight * (low_count * lo + high_count * hi)
                fixed_value += weight * (
                    low_count * (lo - mean) ** 2
                    + high_count * (hi - mean) ** 2
                )
                point.extend([lo] * low_count)
                point.extend([hi] * high_count)
            pivot_weight = Fraction(pivot_weight_int)
            pivot = (total - subtotal) / pivot_weight
            if not pivot_lower <= pivot <= pivot_upper:
                continue
            vertices_checked += 1
            point.append(pivot)
            value = fixed_value + pivot_weight * (pivot - mean) ** 2
            if best is None or value > best:
                best = value
                best_point = point
    if best is None or best_point is None:
        raise ArithmeticError("no exact interval vertex satisfies the trace target")
    return {
        "normalized_mean": mean,
        "normalized_variance_max": best,
        "maximizer": best_point,
        "interval_symmetry_groups": len(groups),
        "feasible_vertices_checked": vertices_checked,
    }


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
    parser.add_argument("--theta-output", type=Path, required=True)
    parser.add_argument("--trace-output", type=Path, required=True)
    parser.add_argument("--trace-target-export", type=Path, required=True)
    parser.add_argument("--exact-row-directory", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--solve", action="store_true")
    parser.add_argument(
        "--shellwise-conservation",
        action="store_true",
        help=(
            "also impose the exact weighted raw-mass identity across all "
            "constituents at every shell and use those tighter intervals"
        ),
    )
    parser.add_argument(
        "--esolver",
        type=Path,
        default=Path("/home/nick/.local/qsopt-ex/usr/bin/esolver"),
    )
    parser.add_argument(
        "--qsopt-library-directory",
        type=Path,
        default=Path("/home/nick/.local/qsopt-ex/usr/lib/x86_64-linux-gnu"),
    )
    args = parser.parse_args()
    if args.p != 11:
        raise ValueError("the cached exact affine reductions are specific to p=11")
    args.work_directory.mkdir(parents=True, exist_ok=True)

    counts = load_counts(args.theta_output)
    trace = parse_named_vector(args.trace_output, "TRACE_HARMONIC_COEFS")
    trace_target = load_trace_target(args.trace_target_export, trace)
    reductions = {}
    for channel in CHANNELS:
        path = (
            args.exact_row_directory
            / f"p11_affine_{channel}_qrows_exact_20260827.txt"
        )
        reductions[channel] = (*load_qrows(path), path)

    cases_output: list[dict[str, object]] = []
    intervals: list[tuple[Fraction, Fraction, int]] = []
    for case in component_cases(args.p):
        channel = str(case["channel"])
        representation_dimension = int(case["representation_dimension"])
        base, matrix, target_base, target, row_path = reductions[channel]
        model = build_model(
            base,
            matrix,
            target_base,
            target,
            counts,
            trace,
            args.p,
            representation_dimension,
        )
        case_name = str(case["name"])
        record: dict[str, object] = {
            **case,
            "exact_row_file": row_path,
            "exact_row_sha256": hashlib.sha256(row_path.read_bytes()).hexdigest(),
            "free_dimension": len(target),
            "constraint_count": len(model.constraints),
            "fixed_checks": model.fixed_checks,
            "target_base": target_base,
            "bounds": {},
        }
        endpoints: dict[str, Fraction] = {}
        for sense in ("minimum", "maximum"):
            stem = f"p11_{case_name}_{sense}".replace("-", "_")
            lp_path = args.work_directory / f"{stem}.lp"
            solution_path = args.work_directory / f"{stem}.sol"
            log_path = args.work_directory / f"{stem}.log"
            objective = write_lp(lp_path, model, sense)
            bound_record: dict[str, object] = {
                "lp_file": lp_path,
                "lp_sha256": hashlib.sha256(lp_path.read_bytes()).hexdigest(),
            }
            if args.solve:
                log_path.write_text(
                    run_qsopt(
                        args.esolver,
                        args.qsopt_library_directory,
                        lp_path,
                        solution_path,
                    )
                )
                verification = verify_certificate(model, objective, solution_path)
                solver_value = Fraction(str(verification["solver_objective"]))
                endpoint = (
                    model.target_base + solver_value
                    if sense == "minimum"
                    else model.target_base - solver_value
                )
                endpoints[sense] = endpoint
                bound_record.update(
                    {
                        "solution_file": solution_path,
                        "solution_sha256": hashlib.sha256(solution_path.read_bytes()).hexdigest(),
                        "log_file": log_path,
                        "target_endpoint": endpoint,
                        "certificate_verification": verification,
                    }
                )
            record["bounds"][sense] = bound_record
        if args.solve:
            if endpoints["minimum"] > endpoints["maximum"]:
                raise ArithmeticError(f"reversed exact interval for {case_name}")
            for _ in range(int(case["component_count"])):
                intervals.append(
                    (
                        endpoints["minimum"],
                        endpoints["maximum"],
                        representation_dimension,
                    )
                )
        cases_output.append(record)

    shellwise_output: list[dict[str, object]] = []
    shellwise_intervals: list[tuple[Fraction, Fraction, int]] = []
    if args.shellwise_conservation:
        pure_reductions = {
            channel: (base, matrix, target_base, target)
            for channel, (base, matrix, target_base, target, _path) in reductions.items()
        }
        cases = component_cases(args.p)
        for distinguished_index, case in enumerate(cases):
            model, representatives = build_shellwise_conserved_model(
                pure_reductions,
                cases,
                distinguished_index,
                counts,
                trace,
                args.p,
            )
            case_name = str(case["name"])
            record = {
                **case,
                "symmetry_representatives": representatives,
                "free_dimension": len(model.target),
                "constraint_count": len(model.constraints),
                "fixed_checks": model.fixed_checks,
                "target_base": model.target_base,
                "bounds": {},
            }
            endpoints: dict[str, Fraction] = {}
            for sense in ("minimum", "maximum"):
                stem = f"p11_shellwise_{case_name}_{sense}".replace("-", "_")
                lp_path = args.work_directory / f"{stem}.lp"
                solution_path = args.work_directory / f"{stem}.sol"
                log_path = args.work_directory / f"{stem}.log"
                objective = write_lp(lp_path, model, sense)
                bound_record: dict[str, object] = {
                    "lp_file": lp_path,
                    "lp_sha256": hashlib.sha256(lp_path.read_bytes()).hexdigest(),
                }
                if args.solve:
                    log_path.write_text(
                        run_qsopt(
                            args.esolver,
                            args.qsopt_library_directory,
                            lp_path,
                            solution_path,
                        )
                    )
                    verification = verify_certificate(model, objective, solution_path)
                    solver_value = Fraction(str(verification["solver_objective"]))
                    endpoint = (
                        model.target_base + solver_value
                        if sense == "minimum"
                        else model.target_base - solver_value
                    )
                    endpoints[sense] = endpoint
                    bound_record.update(
                        {
                            "solution_file": solution_path,
                            "solution_sha256": hashlib.sha256(
                                solution_path.read_bytes()
                            ).hexdigest(),
                            "log_file": log_path,
                            "target_endpoint": endpoint,
                            "certificate_verification": verification,
                        }
                    )
                record["bounds"][sense] = bound_record
            if args.solve:
                if endpoints["minimum"] > endpoints["maximum"]:
                    raise ArithmeticError(
                        f"reversed shellwise interval for {case_name}"
                    )
                for _ in range(int(case["component_count"])):
                    shellwise_intervals.append(
                        (
                            endpoints["minimum"],
                            endpoints["maximum"],
                            int(case["representation_dimension"]),
                        )
                    )
            shellwise_output.append(record)

    output: dict[str, object] = {
        "experiment": "r1_p11_trace_coupled_exact_lp",
        "status": "exact_qsopt_certified" if args.solve else "exact_lps_generated",
        "p": args.p,
        "coefficient_rows_used": min(len(counts), len(trace), 201),
        "trace_target": trace_target,
        "inputs": {
            "theta_output": args.theta_output,
            "theta_output_sha256": hashlib.sha256(args.theta_output.read_bytes()).hexdigest(),
            "trace_output": args.trace_output,
            "trace_output_sha256": hashlib.sha256(args.trace_output.read_bytes()).hexdigest(),
            "trace_target_export": args.trace_target_export,
            "trace_target_export_sha256": hashlib.sha256(
                args.trace_target_export.read_bytes()
            ).hexdigest(),
            "esolver": args.esolver,
            "qsopt_library_directory": args.qsopt_library_directory,
        },
        "cases": cases_output,
    }
    if args.shellwise_conservation:
        output["shellwise_conserved_cases"] = shellwise_output
        output["variance_interval_source"] = "shellwise_conserved_cases"
    else:
        output["variance_interval_source"] = "independent_channel_upper_bounds"
    if args.solve:
        p = args.p
        n = p * p + 1
        zdim = n * (n - 6) // 8
        lbar = Fraction(8 * (n - 2), n - 6)
        spherical = Fraction(8 * n, n + 4)
        variance_intervals = (
            shellwise_intervals if args.shellwise_conservation else intervals
        )
        variance = maximum_weighted_variance_exact(variance_intervals, trace_target)
        scale = -zdim * (lbar - spherical) / trace_target
        phi_variance = scale**2 * variance["normalized_variance_max"]
        phi_endpoint_values = [
            spherical - scale * endpoint
            for lower, upper, _weight in variance_intervals
            for endpoint in (lower, upper)
        ]
        phi_lower = min(phi_endpoint_values)
        phi_upper = max(phi_endpoint_values)
        exact_threshold = Fraction(n, 2) * (lbar - 6) ** 2
        strong_threshold = Fraction(2 * n)
        variance.update(
            {
                "poisson_scale_from_trace": scale,
                "Phi_frobenius_variance_max": phi_variance,
                "Phi_eigenvalue_lower_bound": phi_lower,
                "Phi_eigenvalue_upper_bound": phi_upper,
                "principal_spectral_floor_margin": phi_lower - 6,
                "closes_principal_spectral_floor": phi_lower >= 6,
                "R1_exact_threshold": exact_threshold,
                "R1_exact_margin": exact_threshold - phi_variance,
                "R1_strong_n_over_12_threshold": strong_threshold,
                "R1_strong_margin": strong_threshold - phi_variance,
                "closes_exact_R1": phi_variance <= exact_threshold,
                "closes_strong_R1": phi_variance <= strong_threshold,
            }
        )
        output["variance_bound"] = variance

    args.output.write_text(json.dumps(json_ready(output), indent=2) + "\n")
    print(json.dumps(json_ready(output), indent=2))


if __name__ == "__main__":
    main()
