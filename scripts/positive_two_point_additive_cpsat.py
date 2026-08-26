#!/usr/bin/env python3
"""Exact additive-coefficient model for the positive two-point boundary.

Assume ``|H|=4p+1``, ``boundary(H)={infinity,0}``, and ``c_H=+1``.
Propositions 15.642--15.643 force, in every projective direction ``d``,

    eps_d S_H(z) = 4 + z_{s_d(0)}

on the middle slice.  If ``I`` is the infinity-edge count, ``P_d`` the
parallel finite-edge count, ``n_s`` the infinity-star count in fibre ``s``,
and ``K_st`` the signed finite-edge count between fibres ``s,t``, coefficient
comparison is exactly

    I = 5 + q k0,
    P_d = q kd,
    K_st = eps_d (k0 + kd + delta_sj + delta_tj - n_s - n_t),

where ``q=(p-1)/2`` and ``sum_d kd=8-k0``.  The directional l1 inequality
from Proposition 15.643 restricts the allowed ``(k0,kd)`` pairs.

The model selects H directly, imposes its exact boundary and Paley-sign
product, and enforces every coefficient identity.  Translation symmetry
fixes the finite boundary point at zero.  Therefore INFEASIBLE is an exact
finite certificate for the selected prime/k0 branch.  FEASIBLE is only an
affine-boundary witness, not a full-shell residual witness.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)
from e1_gmin_m4_prop15643 import populated_direction_necessary  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def allowed_k0_values(p: int) -> list[int]:
    """Return k0 values compatible with counts, infinity parity, and l1."""
    q = (p - 1) // 2
    out = []
    for k0 in range(9):
        infinity_edges = 5 + q * k0
        target = 8 - k0
        if infinity_edges > 4 * p + 1 or infinity_edges % 2 == 0:
            continue
        allowed = [
            kd
            for kd in range(target + 1)
            if populated_direction_necessary(p, k0, kd)
        ]
        possible_sums = {0}
        for _ in projective_directions(p):
            possible_sums = {
                old + value
                for old in possible_sums
                for value in allowed
                if old + value <= target
            }
        if target in possible_sums:
            out.append(k0)
    return out


def _bounded_partitions(total: int, length: int, minimum: int, maximum: int):
    """Yield nondecreasing bounded partitions of total with fixed length."""
    if length == 0:
        if total == 0:
            yield ()
        return
    upper = min(maximum, total // length)
    for value in range(minimum, upper + 1):
        for tail in _bounded_partitions(
            total - value, length - 1, value, maximum
        ):
            yield (value, *tail)


def exact_l1_star_profiles(p: int, k0: int, kd: int) -> list[tuple[int, tuple[int, ...]]]:
    """Unordered fibre-count profiles surviving the exact directional l1 cut."""
    q = (p - 1) // 2
    infinity_edges = 5 + q * k0
    transverse_edges = q * (8 - k0 - kd)
    profiles = []
    for special_count in range(min(p, infinity_edges) + 1):
        for other_counts in _bounded_partitions(
            infinity_edges - special_count, p - 1, 0, p
        ):
            counts = (special_count, *other_counts)
            l1 = sum(
                abs(
                    k0 + kd + int(s == 0) + int(t == 0)
                    - counts[s] - counts[t]
                )
                for s, t in itertools.combinations(range(p), 2)
            )
            if l1 <= transverse_edges:
                profiles.append((special_count, other_counts))
    return profiles


def build_geometry(p: int) -> dict:
    q2 = p * p
    directions = projective_directions(p)
    data = [field_direction_data(p, direction) for direction in directions]
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    finite_edges = tuple(itertools.combinations(range(q2), 2))
    signs = tuple(int(C[u + 1, v + 1]) for u, v in finite_edges)
    edge_directions = []
    for u, v in finite_edges:
        matches = [d for d, (_eps, labels) in enumerate(data) if labels[u] == labels[v]]
        if len(matches) != 1:
            raise AssertionError("each affine edge must have exactly one direction")
        if signs[len(edge_directions)] != data[matches[0]][0]:
            raise AssertionError("parallel finite edge has the wrong Paley sign")
        edge_directions.append(matches[0])
    star_signs = tuple(int(C[0, u + 1]) for u in range(q2))
    return {
        "directions": directions,
        "direction_data": data,
        "finite_edges": finite_edges,
        "finite_signs": signs,
        "edge_directions": tuple(edge_directions),
        "star_signs": star_signs,
    }


def verify_witness(p: int, k0: int, geometry: dict, witness: dict) -> bool:
    q = (p - 1) // 2
    q2 = p * p
    star_set = set(witness["infinity_star_points"])
    finite_set = {tuple(edge) for edge in witness["finite_edges"]}
    finite_edges = geometry["finite_edges"]
    signs = geometry["finite_signs"]
    edge_directions = geometry["edge_directions"]
    data = geometry["direction_data"]
    if len(star_set) + len(finite_set) != 4 * p + 1:
        return False
    if len(star_set) != 5 + q * k0:
        return False
    degree = [int(u in star_set) for u in range(q2)]
    negative_count = sum(
        geometry["star_signs"][u] == -1 for u in star_set
    )
    selected_indices = []
    for e, (u, v) in enumerate(finite_edges):
        if (u, v) not in finite_set:
            continue
        selected_indices.append(e)
        degree[u] += 1
        degree[v] += 1
        negative_count += signs[e] == -1
    if {u for u, value in enumerate(degree) if value & 1} != {0}:
        return False
    if len(star_set) % 2 != 1 or negative_count % 2 != 0:
        return False
    kd_values = witness["direction_multiplicities"]
    if sum(kd_values) != 8 - k0:
        return False
    for d, (eps, labels) in enumerate(data):
        kd = kd_values[d]
        if not populated_direction_necessary(p, k0, kd):
            return False
        parallel = [e for e in selected_indices if edge_directions[e] == d]
        if len(parallel) != q * kd:
            return False
        counts = [sum(labels[u] == s for u in star_set) for s in range(p)]
        special = labels[0]
        for s, t in itertools.combinations(range(p), 2):
            signed_cross = sum(
                signs[e]
                for e in selected_indices
                if {labels[finite_edges[e][0]], labels[finite_edges[e][1]]} == {s, t}
            )
            expected = eps * (
                k0
                + kd
                + int(s == special)
                + int(t == special)
                - counts[s]
                - counts[t]
            )
            if signed_cross != expected:
                return False
    return True


def solve_case(
    p: int,
    k0_value: int,
    seconds: float,
    workers: int,
    seed: int,
    fixed_kd: tuple[int, ...] | None = None,
    star_zero: int | None = None,
    fixed_star_in: tuple[int, ...] = (),
    fixed_star_out: tuple[int, ...] = (),
) -> dict:
    from ortools.sat.python import cp_model

    if k0_value not in allowed_k0_values(p):
        raise ValueError(f"k0={k0_value} is not arithmetically allowed for p={p}")
    if fixed_kd is not None and len(fixed_kd) != p + 1:
        raise ValueError(f"fixed_kd must have {p + 1} entries")
    if star_zero not in (None, 0, 1):
        raise ValueError("star_zero must be 0, 1, or None")
    if set(fixed_star_in) & set(fixed_star_out):
        raise ValueError("a point cannot be forced both into and out of the star")
    if not all(0 <= u < p * p for u in (*fixed_star_in, *fixed_star_out)):
        raise ValueError("fixed star point lies outside F_{p^2}")
    started = time.time()
    geometry = build_geometry(p)
    q = (p - 1) // 2
    q2 = p * p
    data = geometry["direction_data"]
    finite_edges = geometry["finite_edges"]
    signs = geometry["finite_signs"]
    edge_directions = geometry["edge_directions"]

    model = cp_model.CpModel()
    star = [model.new_bool_var(f"star_{u}") for u in range(q2)]
    selected = [model.new_bool_var(f"edge_{u}_{v}") for u, v in finite_edges]
    infinity_edges = 5 + q * k0_value
    model.add(sum(star) == infinity_edges)
    if star_zero is not None:
        model.add(star[0] == star_zero)
    for u in fixed_star_in:
        model.add(star[u] == 1)
    for u in fixed_star_out:
        model.add(star[u] == 0)
    model.add(sum(selected) == 4 * p + 1 - infinity_edges)

    # Infinity has odd degree.  At finite vertices the finite graph boundary
    # is the infinity-star symmetric difference with {0}.
    model.add_bool_xor(star)
    incident_by_vertex = [[] for _ in range(q2)]
    for e, (u, v) in enumerate(finite_edges):
        incident_by_vertex[u].append(selected[e])
        incident_by_vertex[v].append(selected[e])
    for u in range(q2):
        literals = (
            [star[u], *incident_by_vertex[u]]
            if u == 0
            else [~star[u], *incident_by_vertex[u]]
        )
        model.add_bool_xor(literals)

    negative = [
        star[u] for u, sign in enumerate(geometry["star_signs"]) if sign == -1
    ] + [selected[e] for e, sign in enumerate(signs) if sign == -1]
    if not negative:
        raise AssertionError("Paley edge set unexpectedly has no negative edges")
    # add_bool_xor enforces odd parity; complementing one literal enforces
    # the required even number of negative selected edges, i.e. c_H=+1.
    model.add_bool_xor([~negative[0], *negative[1:]])

    kd_variables = []
    negative_type_kd = []
    aggregate_allowed_kd = [
        kd
        for kd in range(9 - k0_value)
        if populated_direction_necessary(p, k0_value, kd)
    ]
    allowed_kd = [
        kd
        for kd in aggregate_allowed_kd
        if exact_l1_star_profiles(p, k0_value, kd)
    ]
    if not allowed_kd:
        raise AssertionError("no exact-l1 direction profile survives")
    for d, (eps, labels) in enumerate(data):
        kd = model.new_int_var_from_domain(
            cp_model.Domain.from_values(allowed_kd), f"kd_{d}"
        )
        kd_variables.append(kd)
        if fixed_kd is not None:
            model.add(kd == fixed_kd[d])
        if eps == -1:
            negative_type_kd.append(kd)
        model.add(
            sum(selected[e] for e, direction in enumerate(edge_directions) if direction == d)
            == q * kd
        )
        fibre_star = []
        for s in range(p):
            count = model.new_int_var(0, p, f"star_count_{d}_{s}")
            model.add(count == sum(star[u] for u in range(q2) if labels[u] == s))
            fibre_star.append(count)
        special = labels[0]
        max_fibre_count = min(p, infinity_edges)
        count_indicators = {
            (s, value): model.new_bool_var(f"star_count_choice_{d}_{s}_{value}")
            for s in range(p)
            if s != special
            for value in range(max_fibre_count + 1)
        }
        for s in range(p):
            if s == special:
                continue
            model.add_exactly_one(
                count_indicators[s, value] for value in range(max_fibre_count + 1)
            )
            for value in range(max_fibre_count + 1):
                model.add(fibre_star[s] == value).only_enforce_if(
                    count_indicators[s, value]
                )
        histograms = []
        for value in range(max_fibre_count + 1):
            histogram = model.new_int_var(0, p - 1, f"star_hist_{d}_{value}")
            model.add(
                histogram
                == sum(
                    count_indicators[s, value]
                    for s in range(p)
                    if s != special
                )
            )
            histograms.append(histogram)
        profile_rows = set()
        for profile_kd in allowed_kd:
            for special_count, other_counts in exact_l1_star_profiles(
                p, k0_value, profile_kd
            ):
                histogram = tuple(
                    other_counts.count(value) for value in range(max_fibre_count + 1)
                )
                profile_rows.add((profile_kd, special_count, *histogram))
        model.add_allowed_assignments(
            [kd, fibre_star[special], *histograms], sorted(profile_rows)
        )
        cross_groups = {pair: [] for pair in itertools.combinations(range(p), 2)}
        for e, (u, v) in enumerate(finite_edges):
            s, t = labels[u], labels[v]
            if s != t:
                cross_groups[tuple(sorted((s, t)))].append(e)
        absolute_cross_values = []
        for (s, t), indices in cross_groups.items():
            signed_cross = sum(signs[e] * selected[e] for e in indices)
            cross_value = model.new_int_var(-p * p, p * p, f"K_{d}_{s}_{t}")
            model.add(cross_value == signed_cross)
            model.add(
                cross_value
                == eps * (
                    k0_value + kd + int(s == special) + int(t == special)
                    - fibre_star[s] - fibre_star[t]
                )
            )
            absolute = model.new_int_var(0, p * p, f"K_abs_{d}_{s}_{t}")
            model.add_abs_equality(absolute, cross_value)
            absolute_cross_values.append(absolute)
        # Every unit of |K_st| needs a selected transverse edge.  This is a
        # redundant consequence of the signed-cross equations, but exposes
        # the exact per-fibre l1 obstruction directly to presolve.
        model.add(
            sum(absolute_cross_values)
            <= 4 * p + 1 - infinity_edges - q * kd
        )
    model.add(sum(kd_variables) == 8 - k0_value)
    # A finite edge has the sign of its unique parallel direction.  Infinity
    # edges have sign +1 in the standard Paley normalization, so c_H=+1 also
    # gives this compact direction-count parity cut.
    if q % 2:
        model.add_modulo_equality(0, sum(negative_type_kd), 2)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = seed
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        "experiment": "positive_two_point_additive_cpsat",
        "status": "exact_finite_positive_boundary_relaxation",
        "p": p,
        "k0": k0_value,
        "infinity_edges": infinity_edges,
        "finite_edges_count": 4 * p + 1 - infinity_edges,
        "allowed_kd": allowed_kd,
        "aggregate_allowed_kd": aggregate_allowed_kd,
        "exact_l1_profile_counts": {
            str(kd): len(exact_l1_star_profiles(p, k0_value, kd))
            for kd in aggregate_allowed_kd
        },
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "workers": workers,
        "seed": seed,
        "fixed_direction_multiplicities": list(fixed_kd) if fixed_kd is not None else None,
        "fixed_star_zero": star_zero,
        "fixed_star_in": list(fixed_star_in),
        "fixed_star_out": list(fixed_star_out),
        "model": {
            "selected_object": "H",
            "edge_count": 4 * p + 1,
            "boundary": ["infinity", 0],
            "c_H": 1,
            "translation_symmetry_used": True,
            "all_additive_coefficient_rows": (p + 1) * p * (p - 1) // 2,
        },
    }
    if feasible:
        witness = {
            "infinity_star_points": [u for u, value in enumerate(star) if solver.value(value)],
            "finite_edges": [
                list(edge) for edge, value in zip(finite_edges, selected) if solver.value(value)
            ],
            "direction_multiplicities": [solver.value(kd) for kd in kd_variables],
        }
        out["witness"] = witness
        out["witness_verified"] = verify_witness(p, k0_value, geometry, witness)
        if not out["witness_verified"]:
            raise AssertionError("solver witness failed independent arithmetic verification")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--k0", type=int)
    parser.add_argument("--list-k0", action="store_true")
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=15651001)
    parser.add_argument("--fixed-kd", type=int, nargs="+")
    parser.add_argument("--star-zero", type=int, choices=(0, 1))
    parser.add_argument("--star-in", type=int, nargs="*", default=())
    parser.add_argument("--star-out", type=int, nargs="*", default=())
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    allowed = allowed_k0_values(args.p)
    if args.list_k0:
        print(json.dumps({"p": args.p, "allowed_k0": allowed}, indent=2))
        return
    if args.k0 is None:
        raise SystemExit("--k0 is required unless --list-k0 is used")
    out = solve_case(
        args.p,
        args.k0,
        args.seconds,
        args.workers,
        args.seed,
        tuple(args.fixed_kd) if args.fixed_kd is not None else None,
        args.star_zero,
        tuple(args.star_in),
        tuple(args.star_out),
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
