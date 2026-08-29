#!/usr/bin/env python3
"""Exact affine edge-lift attack on the p=19 all-b2 slack-20 profile.

The profile has phase zero ``{0:5,16:5}`` and phase one ``{2:10}``.
All five phase-zero b=0 directions have pointwise slack zero.  Nine of the
phase-one b=2 directions have pointwise slack ``(t-1)^2``; the remaining
direction has scaled mean 38 and Proposition 15.697 proves its excess lift
is Boolean.  The model imposes every resulting coefficient identity together
with the exact edge, boundary, product-sign, profile, mean, and parallel-count
constraints.

The 3,420-form Boolean classification currently depends on the external
Filmus--Vinciguerra restriction theorem.  Consequently INFEASIBLE is a
profile exclusion conditional on that theorem; FEASIBLE is only an affine
necessary-condition witness.  UNKNOWN is not evidence.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from itertools import combinations

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15632 import field_direction_data, projective_directions  # noqa: E402


P = 19
Q = P * P
N = Q + 1
H_SIZE = 4 * P + 1
INFINITY_DEGREES = (0, 20, 38)
LOW_PHASE_ONE_ROWS = (
    (0, 4, 0),
    (20, 2, 2),
    (38, 2, 4),
)
ZERO_PHASE_ZERO_ROWS = (
    (0, 3, 0),
    (20, 1, 2),
    (38, 1, 4),
)


def mean_only_parallel_rows(mean: int) -> tuple[tuple[int, int], ...]:
    rows = []
    for infinity_degree in INFINITY_DEGREES:
        finite_edges = H_SIZE - infinity_degree
        for parallel in range(finite_edges + 1):
            signed_cross = infinity_degree + P * parallel - 3 * P - mean
            if abs(signed_cross) <= finite_edges - parallel:
                rows.append((infinity_degree, parallel))
    return tuple(rows)


ELEVATED_PHASE_ONE_ROWS = mean_only_parallel_rows(38)
HIGH_PHASE_ZERO_ROWS = mean_only_parallel_rows(40)


def canonical_pair_target(
    constant: int,
    linear: dict[int, int],
    quadratic: dict[tuple[int, int], int],
) -> tuple[tuple[int, ...], int]:
    """Pure-pair target and parallel offset for a slice quadratic."""
    pairs = tuple(combinations(range(P), 2))
    row_sums = {
        s: sum(
            quadratic.get((min(s, t), max(s, t)), 0)
            for t in range(P) if t != s
        )
        for s in range(P)
    }
    shifts = {s: row_sums[s] + 2 * linear.get(s, 0) for s in range(P)}
    targets = []
    for s, t in pairs:
        numerator = quadratic.get((s, t), 0) + shifts[s] + shifts[t]
        if numerator & 1:
            raise ArithmeticError("nonintegral canonical pair target")
        targets.append(numerator // 2)
    constants = {
        -sum(
            targets[pairs.index((min(s, t), max(s, t)))]
            for t in range(P) if t != s
        )
        - linear.get(s, 0)
        + 9 * shifts[s]
        for s in range(P)
    }
    if len(constants) != 1:
        raise ArithmeticError("canonical target has inconsistent slice gauge")
    gauge_constant = constants.pop()
    offset = (
        2 * constant + 3 - sum(targets) - 20 * gauge_constant
    )
    return tuple(targets), offset


def elevated_boolean_forms() -> tuple[tuple[str, tuple[int, ...], int], ...]:
    """All density-5/19 Boolean quadratic lifts after fixing odd pair 0,1."""
    forms = []

    def add_form(
        label: str,
        b_constant: int,
        b_linear: dict[int, int],
        b_quadratic: dict[tuple[int, int], int],
    ) -> None:
        constant = 1 + 2 * b_constant
        linear = {0: -1, 1: -1}
        quadratic = {(0, 1): 2}
        for index, value in b_linear.items():
            linear[index] = linear.get(index, 0) + 2 * value
        for pair, value in b_quadratic.items():
            pair = tuple(sorted(pair))
            quadratic[pair] = quadratic.get(pair, 0) + 2 * value
        target, offset = canonical_pair_target(constant, linear, quadratic)
        forms.append((label, target, offset))

    # B=x_i*x_j.
    for i, j in combinations(range(P), 2):
        add_form(f"pp:{i},{j}", 0, {}, {(i, j): 1})
    # B=x_i*(1-x_j), with i the positive literal.
    for i in range(P):
        for j in range(P):
            if i != j:
                add_form(f"pn:{i},{j}", 0, {i: 1}, {(i, j): -1})
    # B is one on the antipodal patterns 100 and 011 on (i,j,k).
    for i in range(P):
        for j, k in combinations((u for u in range(P) if u != i), 2):
            add_form(
                f"anti:{i},{j},{k}",
                0,
                {i: 1},
                {(j, k): 1, (i, j): -1, (i, k): -1},
            )
    if len(forms) != 3420:
        raise ArithmeticError("elevated Boolean form catalog changed")
    return tuple(forms)


ELEVATED_BOOLEAN_FORMS = elevated_boolean_forms()


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def conference_data() -> tuple[list[tuple[int, int]], list[int]]:
    _q, _mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(P)
    edges: list[tuple[int, int]] = []
    signs: list[int] = []
    for a in range(N):
        for b in range(a + 1, N):
            edges.append((a, b))
            signs.append(1 if a == 0 else int(chi((a - 1) - (b - 1))))
    return edges, signs


def and_var(model, left, right, name: str):
    value = model.new_bool_var(name)
    model.add_multiplication_equality(value, [left, right])
    return value


def solve(
    c_h: int,
    infinity_degree_value: int,
    seconds: float,
    workers: int,
    seed: int,
    form_family: str = "all",
    form_start: int = 0,
    form_end: int | None = None,
) -> dict:
    from ortools.sat.python import cp_model

    if infinity_degree_value not in INFINITY_DEGREES:
        raise ValueError("infinity degree is outside the exact aggregate table")
    forms = tuple(
        row
        for row in ELEVATED_BOOLEAN_FORMS
        if form_family == "all" or row[0].startswith(f"{form_family}:")
    )
    # A form's summed coefficient identity must fit in the available cross
    # edges.  This removes pp at I=28,38 before model construction.
    capacity_forms = []
    for row in forms:
        target_sum = sum(row[1])
        offset = row[2]
        if any(
            (parallel := offset + 9 * gauge - infinity_degree_value) >= 0
            and (cross := H_SIZE - infinity_degree_value - parallel) >= 0
            and abs(target_sum + 171 * gauge - 18 * infinity_degree_value)
            <= cross
            for gauge in range(-H_SIZE, H_SIZE + 1)
        ):
            capacity_forms.append(row)
    forms = tuple(capacity_forms)
    form_end = len(forms) if form_end is None else min(form_end, len(forms))
    if not 0 <= form_start < form_end <= len(forms):
        raise ValueError("empty or invalid Boolean form shard")
    forms = forms[form_start:form_end]
    started = time.time()
    edges, signs = conference_data()
    edge_index = {edge: index for index, edge in enumerate(edges)}
    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"e_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == H_SIZE)

    boundary = [model.new_bool_var(f"boundary_{u}") for u in range(Q)]
    model.add(sum(boundary) == 16)
    incident: list[list[object]] = [[] for _ in range(N)]
    for variable, (a, b) in zip(selected, edges):
        incident[a].append(variable)
        incident[b].append(variable)
    # Infinity is outside the all-finite boundary.
    model.add_bool_xor([*incident[0], model.new_constant(1)])
    for u in range(Q):
        model.add_bool_xor([*incident[u + 1], ~boundary[u]])

    negative = [variable for variable, sign in zip(selected, signs) if sign == -1]
    if c_h == -1:
        model.add_bool_xor(negative)
    else:
        model.add_bool_xor([*negative, model.new_constant(1)])

    infinity_degree = model.new_int_var(0, H_SIZE, "infinity_degree")
    model.add(
        infinity_degree
        == sum(selected[edge_index[(0, u + 1)]] for u in range(Q))
    )
    model.add(infinity_degree == infinity_degree_value)
    phase_one_elevated = []
    elevated_form_index = None
    phase_zero_high = []
    direction_rows = []
    for direction_index, direction in enumerate(projective_directions(P)):
        eps, labels = field_direction_data(P, direction)
        phase = int(-eps * c_h == -1)
        fibres = [[] for _ in range(P)]
        for u, label in enumerate(labels):
            fibres[label].append(u)
        odd = []
        for fibre, points in enumerate(fibres):
            value = model.new_bool_var(f"odd_{direction_index}_{fibre}")
            model.add_bool_xor([*(boundary[u] for u in points), ~value])
            odd.append(value)

        star = [
            sum(selected[edge_index[(0, u + 1)]] for u in points)
            for points in fibres
        ]
        model.add(sum(star) == infinity_degree)
        parallel_expression = sum(
            selected[index]
            for index, (a, b) in enumerate(edges)
            if a != 0 and labels[a - 1] == labels[b - 1]
        )
        parallel = model.new_int_var(0, H_SIZE, f"parallel_{direction_index}")
        model.add(parallel == parallel_expression)

        mean_terms = []
        cross_terms: dict[tuple[int, int], list[tuple[int, object]]] = {
            (s, t): [] for s in range(P) for t in range(s + 1, P)
        }
        for index, ((a, b), sign) in enumerate(zip(edges, signs)):
            if a == 0:
                coefficient = 1
            else:
                s, t = labels[a - 1], labels[b - 1]
                if s == t:
                    coefficient = P
                else:
                    coefficient = -eps * sign
                    key = (s, t) if s < t else (t, s)
                    cross_terms[key].append((eps * sign, selected[index]))
            mean_terms.append(coefficient * selected[index])
        scaled_mean = sum(mean_terms) - 3 * P

        if phase == 1:
            elevated = model.new_bool_var(f"elevated_b2_{direction_index}")
            model.add(sum(odd) == 2)
            model.add(scaled_mean == 18 + 20 * elevated)
            phase_one_elevated.append(elevated)
            if len(phase_one_elevated) == 1:
                # After the square torus fixes this direction, its affine
                # stabilizer induces AGL(1,19) on the fibres and is sharply
                # two-transitive.  Normalize the elevated odd pair.
                model.add(odd[0] == 1)
                model.add(odd[1] == 1)

            gauge = model.new_int_var(-H_SIZE, H_SIZE, f"gauge_{direction_index}")
            for infinity_row, parallel_row, gauge_row in LOW_PHASE_ONE_ROWS:
                if infinity_row == infinity_degree_value:
                    model.add(parallel == parallel_row).only_enforce_if(~elevated)
                    model.add(gauge == gauge_row).only_enforce_if(~elevated)
                    break
            if len(phase_one_elevated) == 1:
                elevated_form_index = model.new_int_var(
                    0, len(forms) - 1, "elevated_boolean_form"
                )
                offset = model.new_int_var(3, 9, "elevated_parallel_offset")
                model.add_element(
                    elevated_form_index,
                    [row[2] for row in forms],
                    offset,
                )
                model.add(
                    parallel == offset + 9 * gauge - infinity_degree
                ).only_enforce_if(elevated)
            else:
                model.add(gauge == 0).only_enforce_if(elevated)
            low_absolute_values = []
            for s in range(P):
                for t in range(s + 1, P):
                    odd_pair = and_var(
                        model, odd[s], odd[t], f"odd_pair_{direction_index}_{s}_{t}"
                    )
                    signed_cross = sum(
                        coefficient * variable
                        for coefficient, variable in cross_terms[s, t]
                    )
                    signed_cross_var = model.new_int_var(
                        -H_SIZE, H_SIZE, f"cross_{direction_index}_{s}_{t}"
                    )
                    absolute_cross = model.new_int_var(
                        0, H_SIZE, f"abs_cross_{direction_index}_{s}_{t}"
                    )
                    model.add(signed_cross_var == signed_cross)
                    model.add_abs_equality(absolute_cross, signed_cross_var)
                    low_absolute_values.append(absolute_cross)
                    model.add(
                        signed_cross
                        == odd_pair + gauge - star[s] - star[t]
                    ).only_enforce_if(~elevated)
                    if len(phase_one_elevated) == 1:
                        column = s * (2 * P - s - 1) // 2 + (t - s - 1)
                        target = model.new_int_var(
                            -5, 5, f"elevated_target_{s}_{t}"
                        )
                        model.add_element(
                            elevated_form_index,
                            [row[1][column] for row in forms],
                            target,
                        )
                        model.add(
                            signed_cross
                            == target + gauge - star[s] - star[t]
                        ).only_enforce_if(elevated)
            model.add(
                sum(low_absolute_values) <= H_SIZE - infinity_degree - parallel
            ).only_enforce_if(~elevated)
            role = "phase_one_b2_one_elevated"
        else:
            high = model.new_bool_var(f"phase_zero_b16_{direction_index}")
            model.add(sum(odd) == 16 * high)
            model.add(scaled_mean == 40 * high)
            phase_zero_high.append(high)

            gauge = model.new_int_var(-H_SIZE, H_SIZE, f"gauge_{direction_index}")
            model.add(gauge == 0).only_enforce_if(high)
            for infinity_row, parallel_row, gauge_row in ZERO_PHASE_ZERO_ROWS:
                if infinity_row == infinity_degree_value:
                    model.add(parallel == parallel_row).only_enforce_if(~high)
                    model.add(gauge == gauge_row).only_enforce_if(~high)
                    break
            model.add_allowed_assignments(
                [infinity_degree, parallel], HIGH_PHASE_ZERO_ROWS
            ).only_enforce_if(high)
            zero_absolute_values = []
            for s in range(P):
                for t in range(s + 1, P):
                    signed_cross = sum(
                        coefficient * variable
                        for coefficient, variable in cross_terms[s, t]
                    )
                    signed_cross_var = model.new_int_var(
                        -H_SIZE, H_SIZE, f"cross_{direction_index}_{s}_{t}"
                    )
                    absolute_cross = model.new_int_var(
                        0, H_SIZE, f"abs_cross_{direction_index}_{s}_{t}"
                    )
                    model.add(signed_cross_var == signed_cross)
                    model.add_abs_equality(absolute_cross, signed_cross_var)
                    zero_absolute_values.append(absolute_cross)
                    model.add(
                        signed_cross == gauge - star[s] - star[t]
                    ).only_enforce_if(~high)
            model.add(
                sum(zero_absolute_values) <= H_SIZE - infinity_degree - parallel
            ).only_enforce_if(~high)
            role = "phase_zero_b0_rigid_or_b16"

        direction_rows.append(
            {"direction": list(direction), "eps": int(eps), "phase": phase, "role": role}
        )

    model.add(sum(phase_one_elevated) == 1)
    # The square torus acts regularly on the ten directions of one type.
    model.add(phase_one_elevated[0] == 1)
    model.add(sum(phase_zero_high) == 5)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    solver.parameters.linearization_level = 0
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        "experiment": "p19_slack20_allb2_lift_cpsat",
        "status": "exact_affine_edge_lift_model",
        "p": P,
        "c_H": c_h,
        "fixed_infinity_degree": infinity_degree_value,
        "phase_zero_profile": {"0": 5, "16": 5},
        "phase_one_profile": {"2": 10},
        "edge_variables": len(selected),
        "active_rigid_coefficient_identities": 14 * (1 + P * (P - 1) // 2),
        "elevated_boolean_form_family": form_family,
        "elevated_boolean_catalog_conditional_on_external_restriction_theorem": True,
        "elevated_boolean_form_shard": [form_start, form_end],
        "elevated_boolean_form_catalog_size": len(forms),
        "direction_rows": direction_rows,
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "workers": workers,
        "seed": seed,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        assert elevated_form_index is not None
        out["elevated_boolean_form"] = forms[
            solver.value(elevated_form_index)
        ][0]
        out["boundary"] = [u for u, value in enumerate(boundary) if solver.value(value)]
        out["chosen_edges_H"] = [
            list(edge) for edge, variable in zip(edges, selected) if solver.value(variable)
        ]
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--infinity-degree", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15697001)
    parser.add_argument("--form-family", choices=("all", "pp", "pn", "anti"), default="all")
    parser.add_argument("--form-start", type=int, default=0)
    parser.add_argument("--form-end", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve(
        args.c_h,
        args.infinity_degree,
        args.seconds,
        args.workers,
        args.seed,
        args.form_family,
        args.form_start,
        args.form_end,
    )
    if args.output is not None:
        atomic_write(args.output, out)
    compact = {key: value for key, value in out.items() if key != "chosen_edges_H"}
    print(json.dumps(compact, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
