#!/usr/bin/env python3
"""Can a residual affine separator have a two-vertex boundary?

Laboratory model only.  The selected variables are the ``4p`` edges of G;
the distinguished edge e=(0,1) is omitted, and H=G union {e}.  It imposes
all affine Max+/Max- score inequalities, ``|boundary(H)|=2``, and a chosen
Paley-sign product c_H.  INFEASIBLE is an exact finite certificate for this
branch; FEASIBLE is only an affine witness, not a full-shell residual.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from residual_affine_johnson_milp import (  # noqa: E402
    affine_shell,
    feature_rows,
    unique_rows,
)
from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)


def solve(
    p: int,
    c_h: int,
    time_limit: float,
    workers: int,
    infinity_boundary: bool = False,
    boundary_vertex: int | None = None,
) -> dict:
    from ortools.sat.python import cp_model

    if c_h not in (-1, 1):
        raise ValueError("c_h must be +/-1")
    started = time.time()
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    Yp = affine_shell(p, 1, C)
    Ym = affine_shell(p, -1, C)
    edges, Fp = feature_rows(Yp, C)
    edges_m, Fm = feature_rows(Ym, C)
    if edges_m != edges:
        raise RuntimeError("edge orders differ")
    edge_index = {edge: j for j, edge in enumerate(edges)}
    ei = edge_index[(0, 1)]
    k = 4 * p

    model = cp_model.CpModel()
    x = [model.new_bool_var(f"e_{a}_{b}") for a, b in edges]
    model.add(sum(x) == k)
    model.add(x[ei] == 0)

    boundary = [model.new_bool_var(f"boundary_{v}") for v in range(len(C))]
    for v in range(len(C)):
        incident = [x[j] for j, (a, b) in enumerate(edges) if a == v or b == v]
        model.add_modulo_equality(
            boundary[v], sum(incident) + (1 if v in (0, 1) else 0), 2
        )
    model.add(sum(boundary) == 2)
    if infinity_boundary:
        model.add(boundary[0] == 1)
    if boundary_vertex is not None:
        if not 1 <= boundary_vertex < len(C):
            raise ValueError("boundary_vertex must be a finite vertex index")
        model.add(boundary[0] == 1)
        model.add(boundary[boundary_vertex] == 1)

    negative = [
        x[j]
        for j, (a, b) in enumerate(edges)
        if j != ei and int(C[a, b]) == -1
    ]
    fixed_negative = 1 if int(C[0, 1]) == -1 else 0
    c_negative = model.new_bool_var("c_H_negative")
    model.add_modulo_equality(c_negative, sum(negative) + fixed_negative, 2)
    model.add(c_negative == (1 if c_h == -1 else 0))

    # Materialize Prop. 15.632's type-split parity costs.  These are
    # redundant consequences of the affine rows, but CP-SAT otherwise has
    # to rediscover their global coupling through thousands of local cuts.
    costs_by_type: dict[int, list] = {-1: [], 1: []}
    for j, direction in enumerate(projective_directions(p)):
        eps, labels = field_direction_data(p, direction)
        fibre_parities = []
        for s in range(p):
            parity = model.new_bool_var(f"fibre_parity_{j}_{s}")
            vertices = [
                boundary[1 + u] for u, label in enumerate(labels) if label == s
            ]
            model.add_modulo_equality(parity, sum(vertices), 2)
            fibre_parities.append(parity)
        b_value = model.new_int_var(0, p, f"odd_fibres_{j}")
        model.add(b_value == sum(fibre_parities))
        cost = model.new_int_var(0, 2 * p, f"parity_floor_{j}")
        table = []
        for infinity_bit in (0, 1):
            for c_bit in (0, 1):
                sign = -eps * (-1 if c_bit else 1)
                if infinity_bit:
                    sign *= eps
                for b in range(p + 1):
                    signed = sign * (-1 if b & 1 else 1)
                    phase = 0 if signed == 1 else 1
                    table.append(
                        [infinity_bit, c_bit, b, scaled_direction_floor(p, b, phase)]
                    )
        model.add_allowed_assignments(
            [boundary[0], c_negative, b_value, cost], table
        )
        costs_by_type[eps].append(cost)
    type_budget = (p + 1) ** 2 // 2
    for eps in (-1, 1):
        model.add(sum(costs_by_type[eps]) <= type_budget)

    # If D={infinity,v} and c_H=+1, every direction has parity x_{s_d(v)}.
    # Its baseline already consumes the complete type-split budget, hence
    # A_d=x_{s_d(v)} pointwise.  Materializing this proved equality is much
    # stronger propagation than asking CP-SAT to rediscover saturation from
    # the score inequalities and global means.
    n_saturated_equalities = 0
    if infinity_boundary and c_h == 1:
        for Y, F, eps in ((Yp, Fp, 1), (Ym, Fm, -1)):
            for y, row in zip(Y, F):
                fe = int(row[ei])
                selected_boundary = [
                    boundary[v]
                    for v in range(1, len(C))
                    if int(y[v]) == 1
                ]
                baseline = sum(selected_boundary)
                model.add(
                    eps * (sum(int(row[j]) * x[j] for j in range(len(x))) + fe)
                    == 3 + 2 * baseline
                )
                n_saturated_equalities += 1

    n_score_constraints = 0
    for F, shell in ((Fp, "plus"), (Fm, "minus")):
        fe = F[:, ei]
        for sign in (-1, 1):
            for row in unique_rows(F[fe == sign]):
                if shell == "plus":
                    bound = 4 if sign == -1 else 2
                    indices = np.flatnonzero(row < 0).tolist()
                    model.add(sum(x[j] for j in indices) <= (k - bound) // 2)
                else:
                    bound = -2 if sign == -1 else -4
                    indices = np.flatnonzero(row > 0).tolist()
                    model.add(sum(x[j] for j in indices) <= (k + bound) // 2)
                n_score_constraints += 1

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(time_limit)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = 15641 + p + (1 if c_h == 1 else 0)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 2
    status = solver.solve(model)
    status_name = solver.status_name(status)
    out = {
        "experiment": "residual_boundary_two_lift_cpsat",
        "status": "affine_model_not_full_shell",
        "p": p,
        "c_H": c_h,
        "solver_status": status_name,
        "feasible": status in (cp_model.OPTIMAL, cp_model.FEASIBLE),
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "n_variables": len(x),
        "n_score_constraints": n_score_constraints,
        "type_budget": type_budget,
        "infinity_boundary": infinity_boundary,
        "boundary_vertex": boundary_vertex,
        "saturated_pointwise_equalities": n_saturated_equalities,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time": solver.wall_time,
        "seconds": round(time.time() - started, 3),
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        xi = np.fromiter((solver.value(v) for v in x), dtype=np.int8)
        out["boundary"] = [v for v, bit in enumerate(boundary) if solver.value(bit)]
        out["chosen_edges_G"] = [list(edge) for edge, bit in zip(edges, xi) if bit]
        out["plus_scores"] = sorted({int(v) for v in (Fp @ xi).tolist()})
        out["minus_scores"] = sorted({int(v) for v in (Fm @ xi).tolist()})
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--infinity-boundary", action="store_true")
    parser.add_argument("--boundary-vertex", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve(
        args.p,
        args.c_h,
        args.time_limit,
        args.workers,
        infinity_boundary=args.infinity_boundary or args.boundary_vertex is not None,
        boundary_vertex=args.boundary_vertex,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
