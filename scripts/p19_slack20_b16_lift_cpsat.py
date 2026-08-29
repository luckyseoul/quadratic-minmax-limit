#!/usr/bin/env python3
"""Exact affine edge-lift attack on the surviving p=19 b=16 profile.

The phase-one floors saturate their 200-unit type budget.  Nine b=2
directions therefore have pointwise slack ``(t-1)^2``.  The b=16 direction
has two integral normal-form orbits from its rank-169 three-layer kernel:
the t=9 values on the three empty fibres are either ``{0,2,2}`` or
``{0,0,4}``.

The model imposes the equivalent inter-fibre coefficient identities, the
exact phase profiles and means, the 77-edge count, edge-product sign, and
the odd-degree boundary.  INFEASIBLE is a rigorous profile exclusion;
FEASIBLE is only an affine necessary-condition witness.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15632 import field_direction_data, projective_directions  # noqa: E402
from e1_gmin_m4_prop15696 import AGGREGATE_ROWS  # noqa: E402


P = 19
Q = P * P
N = Q + 1
H_SIZE = 4 * P + 1
def phase_zero_parallel_rows() -> tuple[tuple[int, ...], ...]:
    rows = []
    for infinity_degree, _gauge, _parallel in AGGREGATE_ROWS:
        finite_edges = H_SIZE - infinity_degree
        for role, base_mean in ((0, 0), (2, 20), (16, 40)):
            for elevated in (0, 1):
                mean = base_mean + 20 * elevated
                for parallel in range(finite_edges + 1):
                    signed_cross = (
                        infinity_degree + P * parallel - 3 * P - mean
                    )
                    if abs(signed_cross) <= finite_edges - parallel:
                        rows.append(
                            (
                                infinity_degree,
                                int(role == 0),
                                int(role == 2),
                                int(role == 16),
                                elevated,
                                parallel,
                            )
                        )
    return tuple(rows)


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
    seconds: float,
    workers: int,
    seed: int,
    fixed_infinity_degree: int | None = None,
    b16_shape: str = "both",
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    edges, signs = conference_data()
    edge_index = {edge: index for index, edge in enumerate(edges)}
    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"e_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == H_SIZE)
    model.add(selected[edge_index[(0, 1)]] == 1)

    boundary = [model.new_bool_var(f"boundary_{u}") for u in range(Q)]
    model.add(sum(boundary) == 16)
    incident: list[list[object]] = [[] for _ in range(N)]
    for variable, (a, b) in zip(selected, edges):
        incident[a].append(variable)
        incident[b].append(variable)
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
    model.add_allowed_assignments(
        [infinity_degree], [[row[0]] for row in AGGREGATE_ROWS]
    )
    if fixed_infinity_degree is not None:
        if fixed_infinity_degree not in {row[0] for row in AGGREGATE_ROWS}:
            raise ValueError("fixed infinity degree is outside the aggregate table")
        model.add(infinity_degree == fixed_infinity_degree)

    phase_one_high = []
    phase_one_four_shape = []
    phase_zero_roles = {0: [], 2: [], 16: []}
    phase_zero_elevated = []
    rigid_identities = 0
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
            high = model.new_bool_var(f"b16_{direction_index}")
            model.add(sum(odd) == 2 + 14 * high)
            model.add(scaled_mean == 18 + 20 * high)
            phase_one_high.append(high)

            four_shape = model.new_bool_var(f"b16_shape_400_{direction_index}")
            model.add(four_shape <= high)
            phase_one_four_shape.append(four_shape)

            special = []
            other_empty = []
            for s in range(P):
                value = model.new_bool_var(f"special_empty_{direction_index}_{s}")
                model.add(value <= high)
                model.add(value + odd[s] <= 1)
                special.append(value)
                other = model.new_bool_var(f"other_empty_{direction_index}_{s}")
                model.add(other <= high)
                model.add(other + odd[s] <= 1)
                model.add(other + value <= 1)
                model.add(other >= high - odd[s] - value)
                other_empty.append(other)
            model.add(sum(special) == high)

            gauge = model.new_int_var(-H_SIZE, H_SIZE, f"gauge_{direction_index}")
            model.add(parallel == 4 + 9 * gauge - infinity_degree)
            model.add_allowed_assignments(
                [infinity_degree, gauge, parallel], AGGREGATE_ROWS
            )
            rigid_identities += 1
            for s in range(P):
                for t in range(s + 1, P):
                    odd_pair = and_var(
                        model, odd[s], odd[t], f"odd_pair_{direction_index}_{s}_{t}"
                    )
                    low_pair = model.new_bool_var(
                        f"low_pair_{direction_index}_{s}_{t}"
                    )
                    model.add(low_pair <= odd_pair)
                    model.add(low_pair + high <= 1)
                    model.add(low_pair >= odd_pair - high)

                    other_empty_pair = model.new_bool_var(
                        f"other_empty_pair_{direction_index}_{s}_{t}"
                    )
                    model.add(other_empty_pair <= high)
                    model.add(other_empty_pair + odd[s] <= 1)
                    model.add(other_empty_pair + odd[t] <= 1)
                    model.add(other_empty_pair + special[s] <= 1)
                    model.add(other_empty_pair + special[t] <= 1)
                    model.add(
                        other_empty_pair
                        >= high - odd[s] - odd[t] - special[s] - special[t]
                    )
                    target_022 = (
                        low_pair - special[s] - special[t] - other_empty_pair
                    )
                    special_odd_s = and_var(
                        model,
                        special[s],
                        odd[t],
                        f"special_odd_s_{direction_index}_{s}_{t}",
                    )
                    special_odd_t = and_var(
                        model,
                        special[t],
                        odd[s],
                        f"special_odd_t_{direction_index}_{s}_{t}",
                    )
                    has_other_empty = model.new_bool_var(
                        f"has_other_empty_{direction_index}_{s}_{t}"
                    )
                    model.add_max_equality(
                        has_other_empty, [other_empty[s], other_empty[t]]
                    )
                    target_400 = special_odd_s + special_odd_t - has_other_empty
                    signed_cross = sum(
                        coefficient * variable
                        for coefficient, variable in cross_terms[s, t]
                    )
                    model.add(
                        signed_cross
                        == target_022 + gauge - star[s] - star[t]
                    ).only_enforce_if(~four_shape)
                    model.add(
                        signed_cross
                        == target_400 + gauge - star[s] - star[t]
                    ).only_enforce_if(four_shape)
                    rigid_identities += 1
            role = "phase_one_b2_or_pointwise_b16"
        else:
            role_vars = {
                b: model.new_bool_var(f"phase_zero_b{b}_{direction_index}")
                for b in (0, 2, 16)
            }
            model.add_exactly_one(list(role_vars.values()))
            model.add(
                sum(odd)
                == 2 * role_vars[2] + 16 * role_vars[16]
            )
            elevated = model.new_bool_var(f"phase_zero_elevated_{direction_index}")
            model.add(
                scaled_mean
                == 20 * role_vars[2] + 40 * role_vars[16] + 20 * elevated
            )
            model.add_allowed_assignments(
                [
                    infinity_degree,
                    role_vars[0],
                    role_vars[2],
                    role_vars[16],
                    elevated,
                    parallel,
                ],
                phase_zero_parallel_rows(),
            )
            for b, value in role_vars.items():
                phase_zero_roles[b].append(value)
            phase_zero_elevated.append(elevated)
            role = "phase_zero_0x5_2x1_16x4"
        direction_rows.append(
            {"direction": list(direction), "eps": int(eps), "phase": phase, "role": role}
        )

    model.add(sum(phase_one_high) == 1)
    if b16_shape == "022":
        model.add(sum(phase_one_four_shape) == 0)
    elif b16_shape == "400":
        model.add(sum(phase_one_four_shape) == 1)
    elif b16_shape != "both":
        raise ValueError("b16_shape must be one of both, 022, 400")
    # The square torus fixes infinity and finite zero and acts regularly on
    # the ten directions of one quadratic type, so normalize the unique
    # b=16 phase-one direction to the first slot of that type.
    model.add(phase_one_high[0] == 1)
    model.add(sum(phase_zero_roles[0]) == 5)
    model.add(sum(phase_zero_roles[2]) == 1)
    model.add(sum(phase_zero_roles[16]) == 4)
    model.add(sum(phase_zero_elevated) == 1)

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
        "experiment": "p19_slack20_b16_lift_cpsat",
        "status": "exact_affine_edge_lift_model",
        "p": P,
        "c_H": c_h,
        "fixed_infinity_degree": fixed_infinity_degree,
        "phase_zero_profile": {"0": 5, "2": 1, "16": 4},
        "phase_one_profile": {"2": 9, "16": 1},
        "b16_shape": b16_shape,
        "edge_variables": len(selected),
        "rigid_coefficient_identities": rigid_identities,
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
        out["boundary"] = [u for u, value in enumerate(boundary) if solver.value(value)]
        out["chosen_edges_H"] = [
            list(edge) for edge, variable in zip(edges, selected) if solver.value(variable)
        ]
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15696001)
    parser.add_argument("--infinity-degree", type=int)
    parser.add_argument("--b16-shape", choices=("both", "022", "400"), default="both")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve(
        args.c_h,
        args.seconds,
        args.workers,
        args.seed,
        args.infinity_degree,
        args.b16_shape,
    )
    if args.output is not None:
        atomic_write(args.output, out)
    compact = {key: value for key, value in out.items() if key != "chosen_edges_H"}
    print(json.dumps(compact, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
