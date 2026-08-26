#!/usr/bin/env python3
"""SCIP MILP for exact p=5 full-shell boundary/profile cases.

The model uses binary edge variables, integer quotients for all parity
equations, and the 260 normalized full-shell score inequalities.  Optional
internal, crossing-degree, and exact crossing-pattern restrictions support
the hard Miquelian-circle divide-and-conquer attack.  A feasible solution is
audited directly; SCIP ``infeasible`` is treated as a finite exact exclusion.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
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
from p5_full_shell_fixed_boundary_cpsat import (  # noqa: E402
    atomic_write,
    audit_witness,
    boundary_edge_stabilizers,
    parity_vector,
    shell_rows,
)
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def solve_case(
    source_path: Path,
    orbit_index: int,
    seconds: float,
    workers: int,
    boundary_internal_edges: int | None = None,
    boundary_cross_edges: int | None = None,
    fixed_internal_edge_indices: tuple[int, ...] | None = None,
    boundary_cross_degrees: tuple[int, ...] | None = None,
    fixed_cross_edge_indices: tuple[int, ...] | None = None,
    required_cross_edge_indices: tuple[int, ...] | None = None,
    enumerate_crossing_patterns: bool = False,
    enumerate_crossing_vertices: tuple[int, ...] | None = None,
    crossing_start: int = 0,
    crossing_stop: int | None = None,
    seconds_per_crossing: float = 1.0,
) -> dict:
    from pyscipopt import Model, SCIP_PARAMSETTING, quicksum

    started = time.time()
    source = json.loads(source_path.read_text())
    if int(source["p"]) != 5:
        raise ValueError("source must be p=5")
    orbit = source["orbits"][orbit_index]
    boundary = tuple(sorted(int(value) for value in orbit["representative_vertices"]))
    if len(boundary) % 2 or len(set(boundary)) != len(boundary):
        raise ValueError("boundary must have distinct vertices and even size")
    boundary_set = set(boundary)
    c_h = int(source["c_H"])
    data = geometry(5, "full")
    edges = data["edges"]
    C = data["C"]
    n = int(data["n"])
    edge_index = {edge: index for index, edge in enumerate(edges)}
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
    outside = [vertex for vertex in range(n) if vertex not in boundary_set]

    fixed_internal_set = (
        set(int(value) for value in fixed_internal_edge_indices)
        if fixed_internal_edge_indices is not None
        else None
    )
    if fixed_internal_set is not None:
        if len(fixed_internal_set) != len(fixed_internal_edge_indices):
            raise ValueError("fixed internal edge indices contain duplicates")
        if not fixed_internal_set <= set(internal_indices):
            raise ValueError("fixed internal edge indices contain a non-internal edge")
        distinguished = edge_index[(0, 1)]
        if distinguished in internal_indices and distinguished not in fixed_internal_set:
            raise ValueError("fixed internal pattern must contain edge (0,1)")
        if (
            boundary_internal_edges is not None
            and len(fixed_internal_set) != boundary_internal_edges
        ):
            raise ValueError("fixed internal pattern contradicts its count")
    fixed_cross_set = (
        set(int(value) for value in fixed_cross_edge_indices)
        if fixed_cross_edge_indices is not None
        else None
    )
    if fixed_cross_set is not None:
        if len(fixed_cross_set) != len(fixed_cross_edge_indices):
            raise ValueError("fixed crossing edge indices contain duplicates")
        if not fixed_cross_set <= set(cross_indices):
            raise ValueError("fixed crossing edge indices contain a non-crossing edge")
        distinguished = edge_index[(0, 1)]
        if distinguished in cross_indices and distinguished not in fixed_cross_set:
            raise ValueError("fixed crossing pattern must contain edge (0,1)")
        if (
            boundary_cross_edges is not None
            and len(fixed_cross_set) != boundary_cross_edges
        ):
            raise ValueError("fixed crossing pattern contradicts its count")
    if enumerate_crossing_patterns and fixed_cross_set is not None:
        raise ValueError("enumeration and a fixed crossing pattern are exclusive")
    required_cross_set = (
        set(int(value) for value in required_cross_edge_indices)
        if required_cross_edge_indices is not None
        else set()
    )
    if len(required_cross_set) != len(required_cross_edge_indices or ()):
        raise ValueError("required crossing edge indices contain duplicates")
    if not required_cross_set <= set(cross_indices):
        raise ValueError("required crossing edge indices contain a non-crossing edge")
    if fixed_cross_set is not None and not required_cross_set <= fixed_cross_set:
        raise ValueError("required crossing edges contradict the fixed crossing pattern")
    if boundary_cross_degrees is not None:
        if len(boundary_cross_degrees) != len(boundary):
            raise ValueError("boundary crossing degrees must follow boundary order")
        if any(not 0 <= value <= len(outside) for value in boundary_cross_degrees):
            raise ValueError("a boundary crossing degree is outside its exact range")
        if (
            boundary_cross_edges is not None
            and sum(boundary_cross_degrees) != boundary_cross_edges
        ):
            raise ValueError("boundary crossing degrees contradict their total")

    model = Model("p5_full_shell_fixed_boundary")
    model.hideOutput(True)
    model.setRealParam("limits/time", float(seconds))
    model.setIntParam("parallel/maxnthreads", max(1, int(workers)))
    model.setPresolve(SCIP_PARAMSETTING.AGGRESSIVE)
    model.setSeparating(SCIP_PARAMSETTING.AGGRESSIVE)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    selected = [model.addVar(vtype="B", name=f"edge_{a}_{b}") for a, b in edges]
    model.addCons(quicksum(selected) == 21)
    model.addCons(selected[edge_index[(0, 1)]] == 1)

    if boundary_internal_edges is not None:
        model.addCons(
            quicksum(selected[index] for index in internal_indices)
            == boundary_internal_edges
        )
    if boundary_cross_edges is not None:
        model.addCons(
            quicksum(selected[index] for index in cross_indices)
            == boundary_cross_edges
        )
    if fixed_internal_set is not None:
        for index in internal_indices:
            model.addCons(selected[index] == int(index in fixed_internal_set))
    if boundary_cross_degrees is not None:
        for vertex, degree in zip(boundary, boundary_cross_degrees):
            model.addCons(
                quicksum(
                    selected[index]
                    for index in cross_indices
                    if vertex in edges[index]
                )
                == degree
            )
    if fixed_cross_set is not None:
        for index in cross_indices:
            model.addCons(selected[index] == int(index in fixed_cross_set))
    for index in required_cross_set:
        model.addCons(selected[index] == 1)

    for vertex in range(n):
        incident = [
            selected[index]
            for index, edge in enumerate(edges)
            if vertex in edge
        ]
        quotient = model.addVar(
            vtype="I", lb=0, ub=(n - 1) // 2, name=f"degree_half_{vertex}"
        )
        model.addCons(
            quicksum(incident) - 2 * quotient == int(vertex in boundary_set)
        )
    negative = [
        selected[index]
        for index, (a, b) in enumerate(edges)
        if int(C[a, b]) == -1
    ]
    sign_half = model.addVar(
        vtype="I", lb=0, ub=len(negative) // 2, name="negative_edge_half"
    )
    model.addCons(quicksum(negative) - 2 * sign_half == int(c_h == -1))

    shell_metadata = {}
    circle_flip_pair_counts = {}
    cross_index_set = set(cross_indices)
    flip_sign = np.ones(len(edges), dtype=np.int8)
    flip_sign[cross_indices] = -1
    score_constraints = 0
    for eps in (-1, 1):
        representatives, normalized = shell_rows(eps)
        parity = parity_vector(representatives, eps, c_h, boundary)
        parity_mass = int(parity.sum())
        bad_masks = normalized < 0
        lifts = []
        for row_index, bad_mask in enumerate(bad_masks):
            lift = model.addVar(vtype="I", lb=0, ub=4, name=f"lift_{eps}_{row_index}")
            bad_indices = np.flatnonzero(bad_mask).tolist()
            model.addCons(
                quicksum(selected[index] for index in bad_indices) + 2 * lift
                == 9 - int(parity[row_index])
            )
            lifts.append(lift)
            score_constraints += 1
        lift_mass = (78 - parity_mass) // 2
        model.addCons(quicksum(lifts) == lift_mass)

        row_lookup = {row.tobytes(): index for index, row in enumerate(normalized)}
        flip_pairs = []
        for first, row in enumerate(normalized):
            second = row_lookup.get((row * flip_sign).tobytes())
            if second is None or first >= second:
                continue
            if int(parity[first]) != int(parity[second]):
                raise AssertionError("circle flip changed shell parity")
            flip_pairs.append((first, second))
            model.addCons(
                quicksum(
                    int(row[index]) * selected[index]
                    for index in range(len(edges))
                    if index not in cross_index_set
                )
                >= 3
            )
        circle_flip_pair_counts[str(eps)] = len(flip_pairs)
        shell_metadata[str(eps)] = {
            "rows": 130,
            "parity_mass": parity_mass,
            "lift_mass": lift_mass,
        }

    direction_rows = []
    half_means = {-1: [], 1: []}
    finite_boundary = tuple(vertex - 1 for vertex in boundary if vertex != 0)
    infinity_value = int(0 in boundary_set)
    for direction_index, direction in enumerate(projective_directions(5)):
        eps, labels = field_direction_data(5, direction)
        counts = [0] * 5
        for value in finite_boundary:
            counts[labels[value]] += 1
        odd_fibres = sum(value & 1 for value in counts)
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
        half_mean = model.addVar(
            vtype="I", lb=floor // 2, ub=105, name=f"half_mean_{direction_index}"
        )
        model.addCons(
            2 * half_mean
            == quicksum(
                coefficient * selected[index]
                for index, coefficient in enumerate(coefficients)
            )
            - 15
        )
        half_means[eps].append(half_mean)
        direction_rows.append(
            {"eps": eps, "odd_fibres": odd_fibres, "phase": phase, "floor": floor}
        )
    for eps in (-1, 1):
        model.addCons(quicksum(half_means[eps]) == 9)

    common_output = {
        "experiment": "p5_full_shell_fixed_boundary_scip",
        "source": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "orbit_index": orbit_index,
        "orbit_size": int(orbit["size"]),
        "boundary": list(boundary),
        "c_H": c_h,
        "edge_variables": len(selected),
        "score_constraints": score_constraints,
        "shells": shell_metadata,
        "circle_flip_pair_counts": circle_flip_pair_counts,
        "direction_rows": direction_rows,
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
        "fixed_cross_edges": (
            [list(edges[index]) for index in sorted(fixed_cross_set)]
            if fixed_cross_set is not None
            else None
        ),
        "required_cross_edges": [
            list(edges[index]) for index in sorted(required_cross_set)
        ],
        "workers": workers,
    }

    if enumerate_crossing_patterns:
        if boundary_cross_degrees is None:
            raise ValueError(
                "crossing-pattern enumeration requires boundary crossing degrees"
            )
        crossing_groups = {}
        for vertex, degree in zip(boundary, boundary_cross_degrees):
            group = tuple(
                index for index in cross_indices if vertex in edges[index]
            )
            crossing_groups[vertex] = frozenset(group)
        enumerated_vertices = (
            tuple(boundary)
            if enumerate_crossing_vertices is None
            else tuple(int(value) for value in enumerate_crossing_vertices)
        )
        if (
            not enumerated_vertices
            or len(set(enumerated_vertices)) != len(enumerated_vertices)
            or not set(enumerated_vertices) <= boundary_set
        ):
            raise ValueError(
                "enumerated crossing vertices must be a nonempty boundary subset"
            )
        degree_by_vertex = dict(zip(boundary, boundary_cross_degrees))
        choices = [
            tuple(
                itertools.combinations(
                    tuple(crossing_groups[vertex]),
                    degree_by_vertex[vertex],
                )
            )
            for vertex in enumerated_vertices
        ]
        raw_total = math.prod(len(values) for values in choices)

        # Keep only exact signed-semilinear symmetries preserving the fixed
        # internal pattern and the ordered crossing-degree profile.
        group_to_vertex = {
            group: vertex for vertex, group in crossing_groups.items()
        }
        enumerated_vertex_set = set(enumerated_vertices)
        stabilizers = []
        for permutation in boundary_edge_stabilizers(boundary):
            if fixed_internal_set is not None and {
                permutation[index] for index in fixed_internal_set
            } != fixed_internal_set:
                continue
            valid = True
            for vertex, group in crossing_groups.items():
                image = frozenset(permutation[index] for index in group)
                target = group_to_vertex.get(image)
                if target is None or degree_by_vertex[target] != degree_by_vertex[vertex]:
                    valid = False
                    break
                if (vertex in enumerated_vertex_set) != (
                    target in enumerated_vertex_set
                ):
                    valid = False
                    break
            if valid and required_cross_set and {
                permutation[index] for index in required_cross_set
            } != required_cross_set:
                valid = False
            if valid:
                stabilizers.append(permutation)
        if not stabilizers:
            raise AssertionError("crossing-profile stabilizer lost the identity")

        representative_weights: dict[tuple[int, ...], int] = {}
        for per_vertex in itertools.product(*choices):
            pattern = tuple(sorted(index for group in per_vertex for index in group))
            canonical = min(
                tuple(sorted(permutation[index] for index in pattern))
                for permutation in stabilizers
            )
            representative_weights[canonical] = (
                representative_weights.get(canonical, 0) + 1
            )
        if sum(representative_weights.values()) != raw_total:
            raise AssertionError("crossing symmetry quotient lost raw patterns")
        representatives = sorted(representative_weights)
        start = max(0, int(crossing_start))
        stop = len(representatives) if crossing_stop is None else min(
            len(representatives), int(crossing_stop)
        )
        if not 0 <= start <= stop:
            raise ValueError("invalid crossing representative shard range")

        infeasible_representatives = 0
        infeasible_raw_patterns = 0
        unknown_cases = []
        witness = None
        node_total = 0
        solve_started = time.time()
        transformed = False
        enumerated_cross_indices = set().union(
            *(crossing_groups[vertex] for vertex in enumerated_vertices)
        )
        for representative_index in range(start, stop):
            if transformed:
                model.freeTransform()
            pattern = representatives[representative_index]
            pattern_set = set(pattern)
            for index in enumerated_cross_indices:
                value = float(index in pattern_set)
                model.chgVarLb(selected[index], value)
                model.chgVarUb(selected[index], value)
            model.setRealParam("limits/time", float(seconds_per_crossing))
            model.optimize()
            transformed = True
            case_status = str(model.getStatus())
            node_total += int(model.getNNodes())
            weight = representative_weights[pattern]
            if case_status == "infeasible":
                infeasible_representatives += 1
                infeasible_raw_patterns += weight
                continue
            n_solutions = int(model.getNSols())
            selected_crossing_edges = [list(edges[index]) for index in pattern]
            if n_solutions == 0:
                unknown_cases.append(
                    {
                        "representative_index": representative_index,
                        "orbit_weight": weight,
                        "solver_status": case_status,
                        "selected_crossing_edges": selected_crossing_edges,
                    }
                )
                continue
            solution = model.getBestSol()
            chosen_edges = [
                list(edge)
                for edge, variable in zip(edges, selected)
                if model.getSolVal(solution, variable) > 0.5
            ]
            witness = {
                "representative_index": representative_index,
                "orbit_weight": weight,
                "selected_crossing_edges": selected_crossing_edges,
                "chosen_edges_H": chosen_edges,
                "witness_audit": audit_witness(data, c_h, boundary, chosen_edges),
            }
            if not witness["witness_audit"]["valid"]:
                raise AssertionError("SCIP crossing witness failed independent audit")
            break
        attempted = (
            infeasible_representatives + len(unknown_cases) + int(witness is not None)
        )
        complete_range = witness is None and attempted == stop - start
        shard_infeasible = complete_range and stop > start and not unknown_cases
        full_profile_infeasible = (
            shard_infeasible
            and start == 0
            and stop == len(representatives)
            and not required_cross_set
        )
        output = {
            **common_output,
            "status": "exact_scip_crossing_orbit_shard",
            "solver_status": (
                "feasible"
                if witness is not None
                else "infeasible"
                if shard_infeasible
                else "unknown"
            ),
            "finite_infeasibility_certificate": shard_infeasible,
            "full_profile_infeasibility_certificate": full_profile_infeasible,
            "feasible": witness is not None,
            "crossing_raw_pattern_count": raw_total,
            "enumerated_crossing_vertices": list(enumerated_vertices),
            "crossing_stabilizer_size": len(stabilizers),
            "crossing_representative_count": len(representatives),
            "crossing_start": start,
            "crossing_stop": stop,
            "crossing_attempted": attempted,
            "infeasible_representatives": infeasible_representatives,
            "infeasible_raw_patterns": infeasible_raw_patterns,
            "unknown_case_count": len(unknown_cases),
            "unknown_cases": unknown_cases,
            "seconds_per_crossing": seconds_per_crossing,
            "node_total": node_total,
            "solve_seconds": time.time() - solve_started,
            "elapsed_seconds": time.time() - started,
        }
        if witness is not None:
            output.update(witness)
        return output

    model.optimize()
    status = str(model.getStatus())
    n_solutions = int(model.getNSols())
    feasible = n_solutions > 0
    output = {
        **common_output,
        "status": "exact_profile_edge_parity_and_shell_lift_milp",
        "solver_status": status,
        "finite_infeasibility_certificate": status == "infeasible",
        "feasible": feasible,
        "n_solutions": n_solutions,
        "nodes": int(model.getNNodes()),
        "gap": float(model.getGap()) if feasible else None,
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        solution = model.getBestSol()
        chosen_edges = [
            list(edge)
            for edge, variable in zip(edges, selected)
            if model.getSolVal(solution, variable) > 0.5
        ]
        output["chosen_edges_H"] = chosen_edges
        output["witness_audit"] = audit_witness(data, c_h, boundary, chosen_edges)
        if not output["witness_audit"]["valid"]:
            raise AssertionError("SCIP witness failed independent audit")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--boundary-internal-edges", type=int)
    parser.add_argument("--boundary-cross-edges", type=int)
    parser.add_argument("--fixed-internal-edge-indices", type=int, nargs="+")
    parser.add_argument("--boundary-cross-degrees", type=int, nargs="+")
    parser.add_argument("--fixed-cross-edge-indices", type=int, nargs="+")
    parser.add_argument("--required-cross-edge-indices", type=int, nargs="+")
    parser.add_argument("--enumerate-crossing-patterns", action="store_true")
    parser.add_argument("--enumerate-crossing-vertices", type=int, nargs="+")
    parser.add_argument("--crossing-start", type=int, default=0)
    parser.add_argument("--crossing-stop", type=int)
    parser.add_argument("--seconds-per-crossing", type=float, default=1.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = solve_case(
        args.source,
        args.orbit_index,
        args.seconds,
        args.workers,
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
        (
            tuple(args.fixed_cross_edge_indices)
            if args.fixed_cross_edge_indices is not None
            else None
        ),
        (
            tuple(args.required_cross_edge_indices)
            if args.required_cross_edge_indices is not None
            else None
        ),
        args.enumerate_crossing_patterns,
        (
            tuple(args.enumerate_crossing_vertices)
            if args.enumerate_crossing_vertices is not None
            else None
        ),
        args.crossing_start,
        args.crossing_stop,
        args.seconds_per_crossing,
    )
    atomic_write(args.output, result)
    print(json.dumps({key: value for key, value in result.items() if key != "chosen_edges_H"}, indent=2))


if __name__ == "__main__":
    main()
