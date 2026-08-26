#!/usr/bin/env python3
"""Exact p=5 full-shell graph solver for one fixed even boundary.

The model selects the 21 edges of the residual graph ``H``, fixes the
distinguished edge, imposes every vertex-degree and Paley-product parity,
and materializes all 260 complete-shell slack equalities

    bad_count(y) = 9 - P(y) - 2 L(y),   0 <= L(y) <= 4.

It also adds the exact directional mean/floor identities from Proposition
15.632.  A SAT/FEASIBLE result is therefore a genuine full-shell residual
witness; an INFEASIBLE result is a rigorous fixed-boundary exclusion.
UNKNOWN has no mathematical force.  Every returned witness is audited from
the selected edge list independently of the CP-SAT variables.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


@functools.lru_cache(maxsize=2)
def shell_rows(eps: int) -> tuple[np.ndarray, np.ndarray]:
    data = geometry(5, "full")
    edges = data["edges"]
    C = data["C"]
    left = np.asarray([a for a, _b in edges], dtype=np.int16)
    right = np.asarray([b for _a, b in edges], dtype=np.int16)
    Y = data["shells"][eps]
    features = (Y[:, left] * Y[:, right] * C[left, right]).astype(np.int8)
    unique, indices, counts = np.unique(
        np.ascontiguousarray(features),
        axis=0,
        return_index=True,
        return_counts=True,
    )
    if unique.shape != (130, 325) or not np.all(counts == 2):
        raise AssertionError("unexpected p=5 full eigenshell structure")
    return Y[indices].astype(np.int8), (eps * unique).astype(np.int8)


def parity_vector(representatives: np.ndarray, eps: int, c_h: int, boundary: tuple[int, ...]) -> np.ndarray:
    products = np.prod(representatives[:, boundary].astype(np.int16), axis=1)
    return (-eps * c_h * products == -1).astype(np.int8)


def audit_witness(
    data: dict,
    c_h: int,
    boundary: tuple[int, ...],
    chosen_edges: list[list[int]],
) -> dict:
    edges = data["edges"]
    C = data["C"]
    chosen = {tuple(int(value) for value in edge) for edge in chosen_edges}
    degrees = [0] * int(data["n"])
    for a, b in chosen:
        degrees[a] += 1
        degrees[b] += 1
    observed_boundary = tuple(index for index, degree in enumerate(degrees) if degree & 1)
    product = math.prod(int(C[a, b]) for a, b in chosen)
    shell_audits = {}
    for eps in (-1, 1):
        representatives, normalized = shell_rows(eps)
        parity = parity_vector(representatives, eps, c_h, boundary)
        bad = (normalized < 0).astype(np.int8)
        selected = np.asarray([int(edge in chosen) for edge in edges], dtype=np.int16)
        bad_counts = bad.astype(np.int16) @ selected
        numerators = 9 - parity.astype(np.int16) - bad_counts
        lifts = numerators // 2
        shell_audits[str(eps)] = {
            "minimum_normalized_score": int((21 - 2 * bad_counts).min()),
            "parity_mass": int(parity.sum()),
            "lift_mass": int(lifts.sum()),
            "expected_lift_mass": int((78 - int(parity.sum())) // 2),
            "all_slack_parities": bool(np.all(numerators % 2 == 0)),
            "all_lifts_in_0_4": bool(np.all((0 <= lifts) & (lifts <= 4))),
        }
    valid = bool(
        len(chosen) == 21
        and (0, 1) in chosen
        and observed_boundary == boundary
        and product == c_h
        and all(
            row["minimum_normalized_score"] >= 3
            and row["all_slack_parities"]
            and row["all_lifts_in_0_4"]
            and row["lift_mass"] == row["expected_lift_mass"]
            for row in shell_audits.values()
        )
    )
    return {
        "valid": valid,
        "edge_count": len(chosen),
        "distinguished_edge_present": (0, 1) in chosen,
        "boundary": list(observed_boundary),
        "c_H": product,
        "shells": shell_audits,
    }


def boundary_edge_stabilizers(boundary: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """Full signed-semilinear edge stabilizer of the fixed problem.

    Besides multiplications and Frobenius, the setwise stabilizer of the
    distinguished edge contains signed inversion.  We admit both Paley
    isometries and anti-isometries when their switching/product factor on
    the requested odd-degree boundary preserves ``c_H``.  Anti-isometries
    merely exchange the two complete eigenshells, both of which are imposed.
    """
    q2, mul, _add, chi, frob, _norm, _ia, _ib = field_ctx(5)
    boundary_set = set(boundary)
    data = geometry(5, "full")
    edges = data["edges"]
    edge_index = {edge: index for index, edge in enumerate(edges)}
    permutations = set()

    def inverse(value: int) -> int:
        if value == 0:
            raise ZeroDivisionError("finite-field inverse of zero")
        out, base, exponent = 1, value, q2 - 2
        while exponent:
            if exponent & 1:
                out = mul(out, base)
            base = mul(base, base)
            exponent >>= 1
        return out

    for alpha in range(1, q2):
        sigma = int(chi(alpha))
        for use_frobenius in (False, True):
            semilinear = tuple(
                frob(value) if use_frobenius else value for value in range(q2)
            )
            for use_inversion in (False, True):
                if not use_inversion:
                    finite = tuple(mul(alpha, value) for value in semilinear)
                    vertex = (0, *(value + 1 for value in finite))
                    # C_{g(a),g(b)} = sigma*u_a*u_b*C_{a,b}, with the
                    # multiplier anti-isometry switched only at infinity.
                    switching = (sigma, *(1 for _value in range(q2)))
                else:
                    image = [1, 0]
                    image.extend(
                        1 + mul(alpha, inverse(semilinear[value]))
                        for value in range(1, q2)
                    )
                    vertex = tuple(image)
                    # Signed inversion contributes chi(x) at finite x!=0;
                    # the following multiplier contributes sigma at the old
                    # zero, whose inverse image is infinity.
                    switching = tuple(
                        sigma if old_vertex == 1
                        else int(chi(old_vertex - 1)) if old_vertex > 1
                        else 1
                        for old_vertex in range(q2 + 1)
                    )
                if {vertex[value] for value in boundary} != boundary_set:
                    continue
                product_factor = sigma  # sigma**21, since |H|=21
                for value in boundary:
                    product_factor *= switching[value]
                if product_factor != 1:
                    continue
                permutations.add(
                    tuple(
                        edge_index[tuple(sorted((vertex[a], vertex[b])))]
                        for a, b in edges
                    )
                )
    return tuple(sorted(permutations))


def add_lex_leader(model, selected, permutation: tuple[int, ...], label: int) -> int:
    """Impose selected <=lex selected[permutation] exactly."""
    if permutation == tuple(range(len(selected))):
        return 0
    prefix = model.new_bool_var(f"lex_{label}_prefix_0")
    model.add(prefix == 1)
    auxiliary = 1
    for index, image_index in enumerate(permutation):
        next_prefix = model.new_bool_var(f"lex_{label}_prefix_{index + 1}")
        # prefix=0 remains 0.  At prefix=1, (1,0) is forbidden; equality
        # keeps the prefix active and (0,1) resolves the comparison.
        model.add_allowed_assignments(
            [prefix, selected[index], selected[image_index], next_prefix],
            [
                [0, 0, 0, 0],
                [0, 0, 1, 0],
                [0, 1, 0, 0],
                [0, 1, 1, 0],
                [1, 0, 0, 1],
                [1, 0, 1, 0],
                [1, 1, 1, 1],
            ],
        )
        prefix = next_prefix
        auxiliary += 1
    return auxiliary


def solve(
    source_path: Path,
    orbit_index: int,
    seconds: float,
    workers: int,
    seed: int,
    symmetry_breaking: bool = False,
    shell_encoding: str = "lift",
    boundary_internal_edges: int | None = None,
    boundary_cross_edges: int | None = None,
    fixed_internal_edge_indices: tuple[int, ...] | None = None,
    boundary_cross_degrees: tuple[int, ...] | None = None,
    outside_cross_odd_vertices: int | None = None,
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    source = json.loads(source_path.read_text())
    if int(source["p"]) != 5:
        raise ValueError("source must be p=5")
    orbit = source["orbits"][orbit_index]
    boundary = tuple(sorted(int(value) for value in orbit["representative_vertices"]))
    if len(boundary) % 2 or len(set(boundary)) != len(boundary):
        raise ValueError("boundary must have distinct vertices and even size")
    c_h = int(source["c_H"])
    if shell_encoding not in {"lift", "xor"}:
        raise ValueError("shell_encoding must be lift or xor")
    data = geometry(5, "full")
    C = data["C"]
    edges = data["edges"]
    n = int(data["n"])
    edge_index = {edge: index for index, edge in enumerate(edges)}

    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"edge_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == 21)
    model.add(selected[edge_index[(0, 1)]] == 1)

    boundary_set = set(boundary)
    internal_indices = [
        index
        for index, (a, b) in enumerate(edges)
        if a in boundary_set and b in boundary_set
    ]
    cross_indices = [
        index
        for index, (a, b) in enumerate(edges)
        if (a in boundary_set) ^ (b in boundary_set)
    ]
    if boundary_internal_edges is not None:
        if not 0 <= boundary_internal_edges <= len(internal_indices):
            raise ValueError("boundary_internal_edges is outside its exact range")
        model.add(
            sum(selected[index] for index in internal_indices)
            == boundary_internal_edges
        )
    if boundary_cross_edges is not None:
        if not 0 <= boundary_cross_edges <= min(21, len(cross_indices)):
            raise ValueError("boundary_cross_edges is outside its exact range")
        model.add(
            sum(selected[index] for index in cross_indices) == boundary_cross_edges
        )
    if boundary_cross_degrees is not None:
        if len(boundary_cross_degrees) != len(boundary):
            raise ValueError("boundary_cross_degrees must follow the boundary order")
        if any(not 0 <= value <= n - len(boundary) for value in boundary_cross_degrees):
            raise ValueError("a boundary crossing degree is outside its exact range")
        if (
            boundary_cross_edges is not None
            and sum(boundary_cross_degrees) != boundary_cross_edges
        ):
            raise ValueError("boundary crossing degrees contradict their total")
        for vertex, degree in zip(boundary, boundary_cross_degrees):
            model.add(
                sum(
                    selected[index]
                    for index in cross_indices
                    if vertex in edges[index]
                )
                == degree
            )
    if outside_cross_odd_vertices is not None:
        outside = [vertex for vertex in range(n) if vertex not in boundary_set]
        if not 0 <= outside_cross_odd_vertices <= len(outside):
            raise ValueError("outside_cross_odd_vertices is outside its exact range")
        cross_parities = []
        for vertex in outside:
            parity = model.new_bool_var(f"outside_cross_parity_{vertex}")
            literals = [
                selected[index]
                for index in cross_indices
                if vertex in edges[index]
            ]
            model.add_bool_xor([*literals, ~parity])
            cross_parities.append(parity)
        model.add(sum(cross_parities) == outside_cross_odd_vertices)
    fixed_internal_set: set[int] | None = None
    if fixed_internal_edge_indices is not None:
        fixed_internal_set = set(int(index) for index in fixed_internal_edge_indices)
        if len(fixed_internal_set) != len(fixed_internal_edge_indices):
            raise ValueError("fixed_internal_edge_indices contains duplicates")
        if not fixed_internal_set <= set(internal_indices):
            raise ValueError("fixed_internal_edge_indices contains a non-internal edge")
        if edge_index[(0, 1)] not in fixed_internal_set:
            raise ValueError("the fixed internal pattern must contain edge (0,1)")
        if (
            boundary_internal_edges is not None
            and len(fixed_internal_set) != boundary_internal_edges
        ):
            raise ValueError("fixed internal pattern contradicts its count")
        for index in internal_indices:
            model.add(selected[index] == int(index in fixed_internal_set))

    incident = [[] for _ in range(n)]
    for index, (a, b) in enumerate(edges):
        incident[a].append(selected[index])
        incident[b].append(selected[index])
    for vertex in range(n):
        literals = incident[vertex]
        if vertex in boundary_set:
            model.add_bool_xor(literals)
        else:
            model.add_bool_xor([~literals[0], *literals[1:]])

    stabilizers = boundary_edge_stabilizers(boundary) if symmetry_breaking else ()
    if fixed_internal_set is not None:
        stabilizers = tuple(
            permutation
            for permutation in stabilizers
            if {permutation[index] for index in fixed_internal_set}
            == fixed_internal_set
        )
    symmetry_auxiliary = sum(
        add_lex_leader(model, selected, permutation, index)
        for index, permutation in enumerate(stabilizers)
    )

    negative = [
        selected[index]
        for index, (a, b) in enumerate(edges)
        if int(C[a, b]) == -1
    ]
    if c_h == -1:
        model.add_bool_xor(negative)
    else:
        model.add_bool_xor([~negative[0], *negative[1:]])

    shell_metadata = {}
    shell_lifts = {}
    circle_flip_pair_counts = {}
    cross_index_set = set(cross_indices)
    flip_sign = np.ones(len(edges), dtype=np.int8)
    flip_sign[cross_indices] = -1
    for eps in (-1, 1):
        representatives, normalized = shell_rows(eps)
        parity = parity_vector(representatives, eps, c_h, boundary)
        parity_mass = int(parity.sum())
        if parity_mass > 78:
            return {
                "experiment": "p5_full_shell_fixed_boundary_cpsat",
                "status": "exact_full_shell_edge_and_slack_model",
                "source": str(source_path),
                "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                "orbit_index": orbit_index,
                "orbit_size": int(orbit["size"]),
                "boundary": list(boundary),
                "c_H": c_h,
                "solver_status": "PARITY_MASS_INFEASIBLE",
                "finite_infeasibility_certificate": True,
                "feasible": False,
                "elapsed_seconds": time.time() - started,
            }
        lift_mass = (78 - parity_mass) // 2
        bad_masks = normalized < 0
        if shell_encoding == "lift":
            lifts = [model.new_int_var(0, 4, f"lift_{eps}_{row}") for row in range(130)]
            for row, bad_mask in enumerate(bad_masks):
                model.add(
                    sum(selected[index] for index in np.flatnonzero(bad_mask).tolist())
                    == 9 - int(parity[row]) - 2 * lifts[row]
                )
            model.add(sum(lifts) == lift_mass)
            shell_lifts[eps] = lifts
        else:
            shell_lifts[eps] = []
            for row, bad_mask in enumerate(bad_masks):
                bad_literals = [
                    selected[index] for index in np.flatnonzero(bad_mask).tolist()
                ]
                model.add(sum(bad_literals) <= 9)
                if int(parity[row]) == 0:
                    model.add_bool_xor(bad_literals)
                else:
                    model.add_bool_xor([~bad_literals[0], *bad_literals[1:]])
        row_lookup = {row.tobytes(): index for index, row in enumerate(normalized)}
        flip_pairs = []
        for first, row in enumerate(normalized):
            second = row_lookup.get((row * flip_sign).tobytes())
            if second is None or first >= second:
                continue
            if int(parity[first]) != int(parity[second]):
                raise AssertionError("circle flip changed fixed-boundary shell parity")
            flip_pairs.append((first, second))
            # Adding the paired score inequality explicitly exposes the
            # cancellation of all crossing edges to presolve.
            model.add(
                sum(
                    int(row[index]) * selected[index]
                    for index in range(len(edges))
                    if index not in cross_index_set
                )
                >= 3
            )
            if shell_encoding == "lift":
                noncross_bad = [
                    index
                    for index in np.flatnonzero(bad_masks[first]).tolist()
                    if index not in cross_index_set
                ]
                model.add(
                    2 * sum(selected[index] for index in noncross_bad)
                    + sum(selected[index] for index in cross_indices)
                    + 2 * shell_lifts[eps][first]
                    + 2 * shell_lifts[eps][second]
                    == 18 - 2 * int(parity[first])
                )
                model.add(
                    sum(
                        (1 if bool(bad_masks[first, index]) else -1)
                        * selected[index]
                        for index in cross_indices
                    )
                    + 2 * shell_lifts[eps][first]
                    - 2 * shell_lifts[eps][second]
                    == 0
                )
        circle_flip_pair_counts[str(eps)] = len(flip_pairs)
        shell_metadata[str(eps)] = {
            "rows": 130,
            "parity_mass": parity_mass,
            "lift_mass": lift_mass,
        }

    # Exact directional means and floor bounds, valid for any fixed even
    # boundary and redundant only at the level of integer consequences.
    direction_data = [
        field_direction_data(5, direction) for direction in projective_directions(5)
    ]
    half_means = {-1: [], 1: []}
    direction_rows = []
    finite_boundary = tuple(vertex - 1 for vertex in boundary if vertex != 0)
    infinity_value = int(0 in boundary_set)
    for direction_index, (eps, labels) in enumerate(direction_data):
        counts = [0] * 5
        for value in finite_boundary:
            counts[labels[value]] += 1
        odd_fibres = sum(count & 1 for count in counts)
        sign = -eps * c_h
        if infinity_value:
            sign *= eps
        if odd_fibres & 1:
            sign *= -1
        phase = int(sign == -1)
        floor = scaled_direction_floor(5, odd_fibres, phase)
        coefficients = []
        for a, b in edges:
            if a == 0:
                coefficient = 1
            else:
                la, lb = labels[a - 1], labels[b - 1]
                coefficient = 5 if la == lb else -eps * int(C[a, b])
            coefficients.append(coefficient)
        half_mean = model.new_int_var(0, 105, f"half_mean_{direction_index}")
        model.add(
            2 * half_mean
            == sum(
                coefficient * selected[index]
                for index, coefficient in enumerate(coefficients)
            )
            - 15
        )
        model.add(2 * half_mean >= floor)
        half_means[eps].append(half_mean)
        direction_rows.append(
            {"eps": eps, "odd_fibres": odd_fibres, "phase": phase, "floor": floor}
        )
    for eps in (-1, 1):
        model.add(2 * sum(half_means[eps]) == 18)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    status_name = solver.status_name(status)
    feasible = status_name in {"OPTIMAL", "FEASIBLE"}
    result = {
        "experiment": "p5_full_shell_fixed_boundary_cpsat",
        "status": "exact_full_shell_edge_and_slack_model",
        "source": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "orbit_index": orbit_index,
        "orbit_size": int(orbit["size"]),
        "boundary": list(boundary),
        "boundary_size": len(boundary),
        "infinity_value": infinity_value,
        "c_H": c_h,
        "solver_status": status_name,
        "finite_infeasibility_certificate": status_name == "INFEASIBLE",
        "feasible": feasible,
        "edge_variables": len(selected),
        "shell_lift_variables": sum(len(values) for values in shell_lifts.values()),
        "shell_encoding": shell_encoding,
        "shell_rows": shell_metadata,
        "circle_flip_pair_counts": circle_flip_pair_counts,
        "direction_rows": direction_rows,
        "workers": workers,
        "seed": seed,
        "symmetry_breaking": symmetry_breaking,
        "boundary_stabilizer_size": len(stabilizers) if symmetry_breaking else None,
        "symmetry_auxiliary_variables": symmetry_auxiliary,
        "boundary_internal_edges": boundary_internal_edges,
        "boundary_cross_edges": boundary_cross_edges,
        "fixed_internal_edges": (
            [list(edges[index]) for index in sorted(fixed_internal_set)]
            if fixed_internal_set is not None
            else None
        ),
        "boundary_cross_degrees": (
            {
                str(vertex): int(degree)
                for vertex, degree in zip(boundary, boundary_cross_degrees)
            }
            if boundary_cross_degrees is not None
            else None
        ),
        "outside_cross_odd_vertices": outside_cross_odd_vertices,
        "conflicts": int(solver.num_conflicts),
        "branches": int(solver.num_branches),
        "wall_time_seconds": float(solver.wall_time),
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        chosen_edges = [
            list(edge)
            for edge, variable in zip(edges, selected)
            if solver.value(variable)
        ]
        result["chosen_edges_H"] = chosen_edges
        result["witness_audit"] = audit_witness(data, c_h, boundary, chosen_edges)
        if not result["witness_audit"]["valid"]:
            raise AssertionError("CP-SAT witness failed independent audit")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=15656001)
    parser.add_argument("--symmetry-breaking", action="store_true")
    parser.add_argument("--shell-encoding", choices=("lift", "xor"), default="lift")
    parser.add_argument("--boundary-internal-edges", type=int)
    parser.add_argument("--boundary-cross-edges", type=int)
    parser.add_argument("--fixed-internal-edge-indices", type=int, nargs="+")
    parser.add_argument("--boundary-cross-degrees", type=int, nargs="+")
    parser.add_argument("--outside-cross-odd-vertices", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = solve(
        args.source,
        args.orbit_index,
        args.seconds,
        args.workers,
        args.seed,
        args.symmetry_breaking,
        args.shell_encoding,
        args.boundary_internal_edges,
        args.boundary_cross_edges,
        (
            tuple(args.fixed_internal_edge_indices)
            if args.fixed_internal_edge_indices is not None
            else None
        ),
        (
            tuple(args.boundary_cross_degrees)
            if args.boundary_cross_degrees is not None
            else None
        ),
        args.outside_cross_odd_vertices,
    )
    atomic_write(args.output, result)
    print(json.dumps({key: value for key, value in result.items() if key != "chosen_edges_H"}, indent=2))


if __name__ == "__main__":
    main()
