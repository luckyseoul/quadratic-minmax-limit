#!/usr/bin/env python3
"""Exact finite scout for the Prop. 15.644 negative-product normal form.

This is a laboratory search, not a proof for general p.  It asks whether a
simple finite-edge graph can satisfy simultaneously

* I=2p-1 infinity edges and E=2p+2 finite edges;
* two exceptional directions of opposite quadratic type, with parallel
  counts (U,V)=(3,1) or (1,3);
* two parallel edges in every baseline direction;
* finite odd boundary S symmetric-difference {v};
* c_H=-1; and
* every exact baseline inter-fibre identity
      K_st = -eps_d (a_s+a_t),
  where a_s=|S cap fibre_s|+1_{s=s_d(v)}-2.

FEASIBLE supplies an exact finite witness.  INFEASIBLE is only a certificate
for the selected prime and exception pair/count split.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def solve(
    p: int,
    seconds: float,
    workers: int,
    *,
    positive_baseline: int | None = None,
    negative_baseline: int | None = None,
    exception_indices: tuple[int, int] | None = None,
    enforce_exception_slack: bool = False,
    enforce_boundary: bool = True,
    enforce_product: bool = True,
    enforce_k: bool = True,
    k_parity_only: bool = False,
) -> dict:
    from ortools.sat.python import cp_model

    if (positive_baseline is None) != (negative_baseline is None):
        raise ValueError("both baseline counts must be supplied together")
    effective_mode = positive_baseline is not None
    if effective_mode:
        assert positive_baseline is not None and negative_baseline is not None
        m = (p + 1) // 2
        finite_edge_count = m * (positive_baseline + negative_baseline) + 2
        infinity_edge_count = 4 * p + 1 - finite_edge_count
        count_splits = [(positive_baseline + 1, negative_baseline + 1)]
    else:
        finite_edge_count = 2 * p + 2
        infinity_edge_count = 2 * p - 1
        count_splits = [(3, 1), (1, 3)]

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
            raise AssertionError("each affine edge must have one projective direction")
        edge_direction.append(matches[0])

    opposite_pairs = [
        (a, b)
        for a, b in itertools.combinations(range(p + 1), 2)
        if data[a][0] != data[b][0]
    ]
    if exception_indices is not None:
        normalized = tuple(sorted(exception_indices))
        if normalized not in opposite_pairs:
            raise ValueError("exception indices must be an opposite-type pair")
        opposite_pairs = [normalized]
    rows = []
    for first, second in opposite_pairs:
        positive_exception = first if data[first][0] == 1 else second
        negative_exception = second if positive_exception == first else first
        for positive_count, negative_count in count_splits:
            model = cp_model.CpModel()
            selected = [model.new_bool_var(f"edge_{u}_{v}") for u, v in edges]
            star = [model.new_bool_var(f"star_{u}") for u in range(q2)]
            model.add(sum(selected) == finite_edge_count)
            model.add(sum(star) == infinity_edge_count)

            required = {
                d: (
                    positive_baseline
                    if data[d][0] == 1
                    else negative_baseline
                )
                for d in range(p + 1)
            } if effective_mode else {d: 2 for d in range(p + 1)}
            required[positive_exception] = positive_count
            required[negative_exception] = negative_count
            for d in range(p + 1):
                model.add(
                    sum(selected[e] for e, de in enumerate(edge_direction) if de == d)
                    == required[d]
                )

            # The finite boundary of the infinity star plus finite graph is {v},
            # normalized to field element zero.
            if enforce_boundary:
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

            # Infinity-edge signs are +1 in this normalization.  c_H=-1 says
            # that the number of selected negative finite edges is odd.
            if enforce_product:
                negative_edges = [selected[e] for e, sign in enumerate(edge_sign) if sign == -1]
                model.add_bool_xor(negative_edges)

            for d, (eps, labels) in enumerate(data):
                if d in (first, second):
                    continue
                special = labels[0]
                parallel_count = required[d]
                twice_c = 2 * (infinity_edge_count + parallel_count - 3) // (p - 1)
                if 2 * (infinity_edge_count + parallel_count - 3) % (p - 1):
                    raise AssertionError("baseline coefficient divisibility failed")
                fibre_star = []
                for s in range(p):
                    count = model.new_int_var(0, p, f"star_count_{d}_{s}")
                    model.add(count == sum(star[u] for u in range(q2) if labels[u] == s))
                    fibre_star.append(count)
                k_absolute_values = []
                if enforce_k:
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
                        absolute = model.new_int_var(0, 2 * p, f"K_abs_{d}_{s}_{t}")
                        model.add_abs_equality(absolute, rhs)
                        k_absolute_values.append(absolute)
                        if k_parity_only:
                            half = model.new_int_var(-2 * p, 2 * p, f"K_half_{d}_{s}_{t}")
                            model.add(signed_cross - rhs == 2 * half)
                        else:
                            model.add(signed_cross == rhs)
                    model.add(
                        sum(k_absolute_values)
                        <= finite_edge_count - parallel_count
                    )

            if enforce_exception_slack:
                middle_weight = (p + 1) // 2
                for d in (first, second):
                    eps, labels = data[d]
                    for chosen in itertools.combinations(range(p), middle_weight):
                        chosen_set = set(chosen)
                        y = [1 if labels[u] in chosen_set else -1 for u in range(q2)]
                        infinity_score = sum(eps * y[u] * star[u] for u in range(q2))
                        finite_score = sum(
                            edge_sign[e] * y[u] * y[v] * selected[e]
                            for e, (u, v) in enumerate(edges)
                        )
                        model.add(eps * (infinity_score + finite_score) >= 3)

            solver = cp_model.CpSolver()
            solver.parameters.max_time_in_seconds = seconds
            solver.parameters.num_search_workers = workers
            solver.parameters.random_seed = 15646 + p + first * (p + 1) + second
            status = solver.solve(model)
            feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
            row = {
                "exception_indices": [first, second],
                "exception_directions": [directions[first], directions[second]],
                "exception_types": [data[first][0], data[second][0]],
                "positive_negative_parallel_counts": [positive_count, negative_count],
                "solver_status": solver.status_name(status),
                "feasible": feasible,
                "wall_time_seconds": solver.wall_time,
            }
            if feasible:
                row["star_point_set"] = [u for u, var in enumerate(star) if solver.value(var)]
                row["finite_edges"] = [list(edges[e]) for e, var in enumerate(selected) if solver.value(var)]
                rows.append(row)
                return {
                    "experiment": "residual_negative_full_cpsat",
                    "status": "finite_full_normal_form_scout_only",
                    "mode": "effective_all_prime" if effective_mode else "prop15644",
                    "positive_negative_baselines": (
                        [positive_baseline, negative_baseline] if effective_mode else None
                    ),
                    "constraints": {
                        "boundary": enforce_boundary,
                        "product": enforce_product,
                        "K": enforce_k,
                        "K_parity_only": k_parity_only,
                        "exception_slack": enforce_exception_slack,
                    },
                    "p": p,
                    "rows": rows,
                    "found_feasible": True,
                    "all_cases_decided": False,
                }
            rows.append(row)

    return {
        "experiment": "residual_negative_full_cpsat",
        "status": "finite_full_normal_form_scout_only",
        "mode": "effective_all_prime" if effective_mode else "prop15644",
        "positive_negative_baselines": (
            [positive_baseline, negative_baseline] if effective_mode else None
        ),
        "p": p,
        "constraints": {
            "boundary": enforce_boundary,
            "product": enforce_product,
            "K": enforce_k,
            "K_parity_only": k_parity_only,
            "exception_slack": enforce_exception_slack,
        },
        "rows": rows,
        "found_feasible": False,
        "all_cases_decided": all(row["solver_status"] != "UNKNOWN" for row in rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=30.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--positive-baseline", type=int)
    parser.add_argument("--negative-baseline", type=int)
    parser.add_argument("--exception-indices", type=int, nargs=2)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--omit-boundary", action="store_true")
    parser.add_argument("--omit-product", action="store_true")
    parser.add_argument("--omit-k", action="store_true")
    parser.add_argument("--k-parity-only", action="store_true")
    parser.add_argument("--enforce-exception-slack", action="store_true")
    args = parser.parse_args()
    out = solve(
        args.p,
        args.seconds,
        args.workers,
        positive_baseline=args.positive_baseline,
        negative_baseline=args.negative_baseline,
        exception_indices=(
            tuple(args.exception_indices) if args.exception_indices is not None else None
        ),
        enforce_boundary=not args.omit_boundary,
        enforce_product=not args.omit_product,
        enforce_k=not args.omit_k,
        k_parity_only=args.k_parity_only,
        enforce_exception_slack=args.enforce_exception_slack,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
