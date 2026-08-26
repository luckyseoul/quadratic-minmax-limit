#!/usr/bin/env python3
"""Exact finite classification of the p=5 negative two-point profiles.

For D={infinity,0} and c_H=-1, every directional slack has the form

    A_d(X) = 1 - x_0 + 2 B_d(X).

At p=5 each quadratic type has three directions and excess budget six.
The same-type directional-mean identity quantizes differences in units six,
so one type has exactly one of two profiles:

* ``unique``: lift masses (sum B) are (3,0,0), and the exceptional
  parallel count is one above the common baseline count;
* ``distributed``: all three lift masses are one, with equal parallel
  counts.

This model enumerates both choices for each type, every admissible parallel
count, and every placement of a unique exceptional direction.  It imposes
the exact 60 affine score equalities, edge counts, boundary XORs, and
negative edge-product parity.  INFEASIBLE is an exact finite certificate for
the corresponding profile; FEASIBLE is only an affine witness.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import field_direction_data, projective_directions  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def count_candidates() -> list[dict]:
    out = []
    for positive_profile, negative_profile in itertools.product(
        ("unique", "distributed"), repeat=2
    ):
        ep = int(positive_profile == "unique")
        en = int(negative_profile == "unique")
        for x in range(8):
            for y in range(8):
                finite_edges = 3 * (x + y) + ep + en
                infinity_edges = 21 - finite_edges
                if infinity_edges < 1 or infinity_edges % 2 == 0:
                    continue
                if (3 * y + en) % 2 != 1:
                    continue
                if infinity_edges - 1 > 2 * finite_edges:
                    continue
                if positive_profile == "unique" and x % 2:
                    continue
                if negative_profile == "unique" and y % 2:
                    continue
                out.append(
                    {
                        "positive_profile": positive_profile,
                        "negative_profile": negative_profile,
                        "positive_parallel_baseline": x,
                        "negative_parallel_baseline": y,
                        "finite_edges": finite_edges,
                        "infinity_edges": infinity_edges,
                    }
                )
    return out


def solve_case(
    candidate: dict,
    positive_exception: int | None,
    negative_exception: int | None,
    seconds: float,
    workers: int,
    seed: int,
    enforce_baseline_k: bool = False,
    enforce_pair_coefficients: bool = False,
) -> dict:
    from ortools.sat.python import cp_model

    p = 5
    q2 = p * p
    directions = projective_directions(p)
    data = [field_direction_data(p, direction) for direction in directions]
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    edges = tuple(itertools.combinations(range(q2), 2))
    edge_sign = [int(C[u + 1, v + 1]) for u, v in edges]
    edge_direction = []
    for u, v in edges:
        matches = [d for d, (_eps, labels) in enumerate(data) if labels[u] == labels[v]]
        if len(matches) != 1:
            raise AssertionError("each affine edge must have one direction")
        edge_direction.append(matches[0])

    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"edge_{u}_{v}") for u, v in edges]
    star = [model.new_bool_var(f"star_{u}") for u in range(q2)]
    model.add(sum(selected) == candidate["finite_edges"])
    model.add(sum(star) == candidate["infinity_edges"])

    required = {}
    lift_mass = {}
    for d, (eps, _labels) in enumerate(data):
        if eps == 1:
            profile = candidate["positive_profile"]
            baseline = candidate["positive_parallel_baseline"]
            exceptional = positive_exception
        else:
            profile = candidate["negative_profile"]
            baseline = candidate["negative_parallel_baseline"]
            exceptional = negative_exception
        if profile == "unique":
            required[d] = baseline + int(d == exceptional)
            lift_mass[d] = 3 if d == exceptional else 0
        else:
            required[d] = baseline
            lift_mass[d] = 1
        model.add(
            sum(selected[e] for e, de in enumerate(edge_direction) if de == d)
            == required[d]
        )

    # The finite graph boundary is star symmetric-difference {0}.
    for u in range(q2):
        incident = [
            selected[e]
            for e, (a, b) in enumerate(edges)
            if a == u or b == u
        ]
        if u == 0:
            model.add_bool_xor([star[u], *incident])
        else:
            model.add_bool_xor([~star[u], *incident])

    negative_edges = [selected[e] for e, sign in enumerate(edge_sign) if sign == -1]
    model.add_bool_xor(negative_edges)

    # Materialize the additive inter-fibre matrices in every zero-lift
    # direction.  These are redundant consequences of the ten exact score
    # rows below, but expose their strongest global linear coupling directly.
    if enforce_baseline_k:
        for d, (eps, labels) in enumerate(data):
            if lift_mass[d] != 0:
                continue
            special = labels[0]
            numerator = 2 * (candidate["infinity_edges"] + required[d] - 3)
            if numerator % (p - 1):
                raise AssertionError("baseline inter-fibre coefficient is not integral")
            twice_c = numerator // (p - 1)
            fibre_star = []
            for s in range(p):
                count = model.new_int_var(0, p, f"star_count_{d}_{s}")
                model.add(count == sum(star[u] for u in range(q2) if labels[u] == s))
                fibre_star.append(count)
            absolute_values = []
            for s, t in itertools.combinations(range(p), 2):
                signed_cross = sum(
                    edge_sign[e] * selected[e]
                    for e, (u, v) in enumerate(edges)
                    if {labels[u], labels[v]} == {s, t}
                )
                rhs = eps * (
                    twice_c
                    - (1 if s == special else 0)
                    - (1 if t == special else 0)
                    - fibre_star[s]
                    - fibre_star[t]
                )
                model.add(signed_cross == rhs)
                absolute = model.new_int_var(0, 2 * p, f"K_abs_{d}_{s}_{t}")
                model.add_abs_equality(absolute, rhs)
                absolute_values.append(absolute)
            model.add(
                sum(absolute_values) <= candidate["finite_edges"] - required[d]
            )

    for d, (eps, labels) in enumerate(data):
        special = labels[0]
        lifts = {}
        for chosen in itertools.combinations(range(p), 3):
            chosen_set = set(chosen)
            y = [1 if labels[u] in chosen_set else -1 for u in range(q2)]
            infinity_score = sum(eps * y[u] * star[u] for u in range(q2))
            finite_score = sum(
                edge_sign[e] * y[u] * y[v] * selected[e]
                for e, (u, v) in enumerate(edges)
            )
            B = model.new_int_var(0, lift_mass[d], f"lift_{d}_{'_'.join(map(str, chosen))}")
            baseline_score = 3 if special in chosen_set else 5
            model.add(eps * (infinity_score + finite_score) == baseline_score + 4 * B)
            lifts[chosen] = B
        model.add(sum(lifts.values()) == lift_mass[d])

        if enforce_pair_coefficients:
            # J(5,3) has ten points and the ten pair monomials form a basis.
            # Invert the lift values and equate pair coefficients of the
            # exact score polynomial.  This is redundant with the ten rows.
            vertex_mass = []
            for i in range(p):
                value = model.new_int_var(0, lift_mass[d], f"lift_U_{d}_{i}")
                model.add(value == sum(B for X, B in lifts.items() if i in X))
                vertex_mass.append(value)
            lift_coefficient_six = {}
            for i, j in itertools.combinations(range(p), 2):
                value = model.new_int_var(0, lift_mass[d], f"lift_T_{d}_{i}_{j}")
                model.add(
                    value
                    == sum(B for X, B in lifts.items() if i in X and j in X)
                )
                coefficient = model.new_int_var(
                    -6 * lift_mass[d],
                    8 * lift_mass[d],
                    f"lift_c6_{d}_{i}_{j}",
                )
                model.add(
                    coefficient
                    == 6 * value
                    - 3 * vertex_mass[i]
                    - 3 * vertex_mass[j]
                    + 2 * lift_mass[d]
                )
                lift_coefficient_six[i, j] = coefficient

            fibre_star = []
            for s in range(p):
                count = model.new_int_var(0, p, f"score_star_count_{d}_{s}")
                model.add(count == sum(star[u] for u in range(q2) if labels[u] == s))
                fibre_star.append(count)
            cross = {
                (s, t): sum(
                    edge_sign[e] * selected[e]
                    for e, (u, v) in enumerate(edges)
                    if {labels[u], labels[v]} == {s, t}
                )
                for s, t in itertools.combinations(range(p), 2)
            }
            row_sum = [
                sum(cross[tuple(sorted((s, t)))] for t in range(p) if t != s)
                for s in range(p)
            ]
            total_cross = sum(cross.values())
            constant_score = (
                -candidate["infinity_edges"] + required[d] + eps * total_cross
            )
            for s, t in itertools.combinations(range(p), 2):
                delta = int(s == special or t == special)
                model.add(
                    6 * (fibre_star[s] + fibre_star[t])
                    - 6 * eps * (row_sum[s] + row_sum[t])
                    + 24 * eps * cross[s, t]
                    + 2 * constant_score
                    == 10 - 6 * delta + 4 * lift_coefficient_six[s, t]
                )

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = seed
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        **candidate,
        "positive_exception": positive_exception,
        "negative_exception": negative_exception,
        "parallel_counts": required,
        "lift_masses": lift_mass,
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "wall_time_seconds": solver.wall_time,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "constraints": {
            "baseline_K": enforce_baseline_k,
            "pair_coefficients": enforce_pair_coefficients,
        },
    }
    if feasible:
        out["star_point_set"] = [u for u, value in enumerate(star) if solver.value(value)]
        out["finite_edge_set"] = [
            list(edge) for edge, value in zip(edges, selected) if solver.value(value)
        ]
    return out


def solve_all(seconds: float, workers: int) -> dict:
    started = time.time()
    p = 5
    data = [field_direction_data(p, direction) for direction in projective_directions(p)]
    by_type = {
        eps: [d for d, (kind, _labels) in enumerate(data) if kind == eps]
        for eps in (-1, 1)
    }
    rows = []
    for candidate_index, candidate in enumerate(count_candidates()):
        positive_choices = (
            by_type[1] if candidate["positive_profile"] == "unique" else [None]
        )
        negative_choices = (
            by_type[-1] if candidate["negative_profile"] == "unique" else [None]
        )
        for positive_exception, negative_exception in itertools.product(
            positive_choices, negative_choices
        ):
            row = solve_case(
                candidate,
                positive_exception,
                negative_exception,
                seconds,
                workers,
                1565000 + 100 * candidate_index + 10 * (positive_exception or 0) + (negative_exception or 0),
            )
            rows.append(row)
    return {
        "experiment": "p5_negative_profile_cpsat",
        "status": "exact_finite_affine_profile_classification",
        "p": p,
        "count_candidates": count_candidates(),
        "rows": rows,
        "case_count": len(rows),
        "all_decided": all(row["solver_status"] != "UNKNOWN" for row in rows),
        "all_infeasible": all(row["solver_status"] == "INFEASIBLE" for row in rows),
        "feasible_count": sum(row["feasible"] for row in rows),
        "unknown_count": sum(row["solver_status"] == "UNKNOWN" for row in rows),
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=30.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve_all(args.seconds, args.workers)
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
