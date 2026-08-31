#!/usr/bin/env python3
"""Native-XOR CryptoMiniSat model for one p=5 full-shell boundary case.

The Boolean edge model is equivalent to the full-shell feasibility model in
``p5_full_shell_fixed_boundary_cpsat.py``.  Cardinality constraints are
encoded to CNF, while every vertex, Paley-product, and shell parity equation
is passed to CryptoMiniSat as a native XOR clause.  This lets CryptoMiniSat
apply Gaussian elimination instead of rediscovering the parity structure
through Tseitin chains.

SAT witnesses are audited independently from their selected edge lists.
UNSAT is an exact finite exclusion; a timeout/UNKNOWN has no mathematical
force.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p5_full_shell_fixed_boundary_cpsat import (  # noqa: E402
    atomic_write,
    audit_witness,
    parity_vector,
    shell_rows,
)
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


CARDINALITY_ENCODINGS = {
    "seqcounter": "seqcounter",
    "sortnetwrk": "sortnetwrk",
    "cardnetwrk": "cardnetwrk",
    "totalizer": "totalizer",
    "mtotalizer": "mtotalizer",
    "kmtotalizer": "kmtotalizer",
}


def _imports():
    """Load PySAT from the project venv and the distro CMS binding."""
    try:
        from pysat.card import CardEnc, EncType
    except ModuleNotFoundError as exc:  # pragma: no cover - deployment guard
        raise RuntimeError(
            "python-sat is required; run with /home/nick/.venvs/mo-exact/bin/python"
        ) from exc
    try:
        from pycryptosat import Solver
    except ModuleNotFoundError:
        # Ubuntu installs the CPython binding here, while isolated venvs omit
        # distro packages from sys.path even when the ABI is identical.
        distro = "/usr/lib/python3/dist-packages"
        if distro not in sys.path:
            sys.path.append(distro)
        from pycryptosat import Solver
    return CardEnc, EncType, Solver


def solve_case(
    source_path: Path,
    orbit_index: int,
    seconds: float,
    threads: int,
    card_encoding: str,
    boundary_internal_edges: int | None = None,
    boundary_cross_edges: int | None = None,
    fixed_internal_edge_indices: tuple[int, ...] | None = None,
    boundary_cross_degrees: tuple[int, ...] | None = None,
    outside_cross_odd_vertices: int | None = None,
    enumerate_crossing_patterns: bool = False,
    crossing_start: int = 0,
    crossing_stop: int | None = None,
    seconds_per_crossing: float = 1.0,
) -> dict:
    CardEnc, EncType, Solver = _imports()
    encoding = getattr(EncType, CARDINALITY_ENCODINGS[card_encoding])
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
    C = data["C"]
    edges = data["edges"]
    n = int(data["n"])
    edge_index = {edge: index for index, edge in enumerate(edges)}
    edge_literals = list(range(1, len(edges) + 1))
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

    if boundary_internal_edges is not None and not (
        0 <= boundary_internal_edges <= len(internal_indices)
    ):
        raise ValueError("boundary_internal_edges is outside its exact range")
    if boundary_cross_edges is not None and not (
        0 <= boundary_cross_edges <= min(21, len(cross_indices))
    ):
        raise ValueError("boundary_cross_edges is outside its exact range")
    if boundary_cross_degrees is not None:
        if len(boundary_cross_degrees) != len(boundary):
            raise ValueError("boundary_cross_degrees must follow boundary order")
        if any(not 0 <= value <= len(outside) for value in boundary_cross_degrees):
            raise ValueError("a boundary crossing degree is outside its exact range")
        if (
            boundary_cross_edges is not None
            and sum(boundary_cross_degrees) != boundary_cross_edges
        ):
            raise ValueError("boundary crossing degrees contradict their total")
    fixed_internal_set: set[int] | None = None
    if fixed_internal_edge_indices is not None:
        fixed_internal_set = set(int(value) for value in fixed_internal_edge_indices)
        if len(fixed_internal_set) != len(fixed_internal_edge_indices):
            raise ValueError("fixed internal edge indices contain duplicates")
        if not fixed_internal_set <= set(internal_indices):
            raise ValueError("fixed internal edge indices contain a non-internal edge")
        if edge_index[(0, 1)] not in fixed_internal_set:
            raise ValueError("the fixed internal pattern must contain edge (0,1)")
        if (
            boundary_internal_edges is not None
            and len(fixed_internal_set) != boundary_internal_edges
        ):
            raise ValueError("fixed internal pattern contradicts its edge count")
    if outside_cross_odd_vertices is not None and not (
        0 <= outside_cross_odd_vertices <= len(outside)
    ):
        raise ValueError("outside_cross_odd_vertices is outside its exact range")

    # Allocate all non-encoding variables before CardEnc starts assigning IDs.
    top_id = len(edge_literals)
    outside_parity_literals: dict[int, int] = {}
    if outside_cross_odd_vertices is not None:
        for vertex in outside:
            top_id += 1
            outside_parity_literals[vertex] = top_id

    solver = Solver(verbose=0, threads=max(1, int(threads)))
    clause_count = 0
    xor_count = 0
    xor_literals = 0
    cardinality_count = 0

    def add_clauses(clauses: list[list[int]]) -> None:
        nonlocal clause_count
        if clauses:
            solver.add_clauses(clauses)
            clause_count += len(clauses)

    def add_unit(literal: int) -> None:
        add_clauses([[int(literal)]])

    def add_atmost(literals: list[int], bound: int) -> None:
        nonlocal top_id, cardinality_count
        bound = int(bound)
        if bound < 0:
            add_clauses([[]])
            cardinality_count += 1
            return
        if bound >= len(literals):
            return
        encoded = CardEnc.atmost(
            lits=literals,
            bound=bound,
            top_id=top_id,
            encoding=encoding,
        )
        top_id = max(top_id, int(encoded.nv))
        add_clauses(encoded.clauses)
        cardinality_count += 1

    def add_exact(literals: list[int], bound: int) -> None:
        nonlocal top_id, cardinality_count
        bound = int(bound)
        if not 0 <= bound <= len(literals):
            add_clauses([[]])
            cardinality_count += 1
            return
        encoded = CardEnc.equals(
            lits=literals,
            bound=bound,
            top_id=top_id,
            encoding=encoding,
        )
        top_id = max(top_id, int(encoded.nv))
        add_clauses(encoded.clauses)
        cardinality_count += 1

    def add_xor(literals: list[int], rhs: int) -> None:
        nonlocal xor_count, xor_literals
        if not literals:
            if rhs:
                add_clauses([[]])
            return
        solver.add_xor_clause([int(value) for value in literals], bool(rhs))
        xor_count += 1
        xor_literals += len(literals)

    add_exact(edge_literals, 21)
    add_unit(edge_literals[edge_index[(0, 1)]])
    if boundary_internal_edges is not None:
        add_exact(
            [edge_literals[index] for index in internal_indices],
            boundary_internal_edges,
        )
    if boundary_cross_edges is not None:
        add_exact(
            [edge_literals[index] for index in cross_indices],
            boundary_cross_edges,
        )
    if fixed_internal_set is not None:
        for index in internal_indices:
            add_unit(
                edge_literals[index] if index in fixed_internal_set else -edge_literals[index]
            )
    if boundary_cross_degrees is not None:
        for vertex, degree in zip(boundary, boundary_cross_degrees):
            add_exact(
                [
                    edge_literals[index]
                    for index in cross_indices
                    if vertex in edges[index]
                ],
                degree,
            )
    if outside_cross_odd_vertices is not None:
        for vertex in outside:
            crossing = [
                edge_literals[index]
                for index in cross_indices
                if vertex in edges[index]
            ]
            # p_v = XOR(crossing), written XOR(crossing,p_v)=0.
            add_xor([*crossing, outside_parity_literals[vertex]], 0)
        add_exact(
            list(outside_parity_literals.values()),
            outside_cross_odd_vertices,
        )

    incident = [[] for _ in range(n)]
    for index, (a, b) in enumerate(edges):
        incident[a].append(edge_literals[index])
        incident[b].append(edge_literals[index])
    for vertex in range(n):
        add_xor(incident[vertex], int(vertex in boundary_set))
    negative = [
        edge_literals[index]
        for index, (a, b) in enumerate(edges)
        if int(C[a, b]) == -1
    ]
    add_xor(negative, int(c_h == -1))

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
        if parity_mass > 78:
            return {
                "experiment": "p5_full_shell_fixed_boundary_cryptominisat",
                "solver_status": "PARITY_MASS_INFEASIBLE",
                "finite_infeasibility_certificate": True,
                "feasible": False,
                "elapsed_seconds": time.time() - started,
            }
        bad_masks = normalized < 0
        for row_index, bad_mask in enumerate(bad_masks):
            bad = [
                edge_literals[index]
                for index in np.flatnonzero(bad_mask).tolist()
            ]
            add_atmost(bad, 9)
            add_xor(bad, 1 - int(parity[row_index]))
            score_constraints += 1

        row_lookup = {row.tobytes(): index for index, row in enumerate(normalized)}
        flip_pairs = []
        for first, row in enumerate(normalized):
            second = row_lookup.get((row * flip_sign).tobytes())
            if second is None or first >= second:
                continue
            if int(parity[first]) != int(parity[second]):
                raise AssertionError("circle flip changed fixed-boundary shell parity")
            flip_pairs.append((first, second))
            if boundary_cross_edges is not None:
                # Adding the two score inequalities cancels all crossing
                # coefficients: 2*bad_noncross + cross_count <= 18.
                noncross_bad = [
                    edge_literals[index]
                    for index in np.flatnonzero(bad_masks[first]).tolist()
                    if index not in cross_index_set
                ]
                add_atmost(noncross_bad, (18 - boundary_cross_edges) // 2)
        circle_flip_pair_counts[str(eps)] = len(flip_pairs)
        shell_metadata[str(eps)] = {
            "rows": 130,
            "parity_mass": parity_mass,
            "lift_mass": (78 - parity_mass) // 2,
        }

    build_seconds = time.time() - started
    common_output = {
        "experiment": "p5_full_shell_fixed_boundary_cryptominisat",
        "source": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "orbit_index": orbit_index,
        "orbit_size": int(orbit["size"]),
        "boundary": list(boundary),
        "c_H": c_h,
        "solver": "cryptominisat-5.11-native-xor",
        "edge_variables": len(edge_literals),
        "total_variables": top_id,
        "clauses": clause_count,
        "native_xor_constraints": xor_count,
        "native_xor_literals": xor_literals,
        "cardinality_constraints": cardinality_count,
        "cardinality_encoding": card_encoding,
        "score_constraints": score_constraints,
        "shells": shell_metadata,
        "circle_flip_pair_counts": circle_flip_pair_counts,
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
        "threads": threads,
        "build_seconds": build_seconds,
    }

    if enumerate_crossing_patterns:
        if boundary_cross_degrees is None:
            raise ValueError(
                "crossing-pattern enumeration requires boundary crossing degrees"
            )
        choices = []
        for vertex, degree in zip(boundary, boundary_cross_degrees):
            incident_crossing = [
                edge_literals[index]
                for index in cross_indices
                if vertex in edges[index]
            ]
            choices.append(tuple(itertools.combinations(incident_crossing, degree)))
        total_patterns = math.prod(len(values) for values in choices)
        start = max(0, int(crossing_start))
        stop = total_patterns if crossing_stop is None else min(
            total_patterns, int(crossing_stop)
        )
        if not 0 <= start <= stop:
            raise ValueError("invalid crossing-pattern shard range")
        infeasible_count = 0
        unknown_cases = []
        witness = None
        iterator = itertools.islice(itertools.product(*choices), start, stop)
        solve_started = time.time()
        for case_index, per_vertex in enumerate(iterator, start=start):
            assumptions = sorted(
                literal for selected_group in per_vertex for literal in selected_group
            )
            satisfiable, assignment = solver.solve(
                assumptions=assumptions,
                time_limit=float(seconds_per_crossing),
            )
            if satisfiable is False:
                infeasible_count += 1
                continue
            selected_crossing_edges = [
                list(edges[literal - 1]) for literal in assumptions
            ]
            if satisfiable is None:
                unknown_cases.append(
                    {
                        "case_index": case_index,
                        "selected_crossing_edges": selected_crossing_edges,
                    }
                )
                continue
            chosen_edges = [
                list(edge)
                for edge, literal in zip(edges, edge_literals)
                if assignment[literal]
            ]
            witness = {
                "case_index": case_index,
                "selected_crossing_edges": selected_crossing_edges,
                "chosen_edges_H": chosen_edges,
                "witness_audit": audit_witness(data, c_h, boundary, chosen_edges),
            }
            if not witness["witness_audit"]["valid"]:
                raise AssertionError("CryptoMiniSat crossing witness failed audit")
            break
        attempted = infeasible_count + len(unknown_cases) + int(witness is not None)
        complete_range = witness is None and attempted == stop - start
        shard_infeasible = complete_range and not unknown_cases
        output = {
            **common_output,
            "status": "exact_incremental_crossing_pattern_shard",
            "solver_status": (
                "SATISFIABLE"
                if witness is not None
                else "UNSATISFIABLE"
                if shard_infeasible
                else "UNKNOWN"
            ),
            "finite_infeasibility_certificate": shard_infeasible,
            "feasible": witness is not None,
            "crossing_pattern_total": total_patterns,
            "crossing_start": start,
            "crossing_stop": stop,
            "crossing_attempted": attempted,
            "crossing_infeasible": infeasible_count,
            "unknown_case_count": len(unknown_cases),
            "unknown_cases": unknown_cases,
            "seconds_per_crossing": seconds_per_crossing,
            "solve_seconds": time.time() - solve_started,
            "elapsed_seconds": time.time() - started,
        }
        if witness is not None:
            output.update(witness)
        return output

    satisfiable, assignment = solver.solve(time_limit=float(seconds))
    status = (
        "SATISFIABLE"
        if satisfiable is True
        else "UNSATISFIABLE"
        if satisfiable is False
        else "UNKNOWN"
    )
    output = {
        **common_output,
        "status": "exact_native_xor_full_shell_cnf_model",
        "solver_status": status,
        "finite_infeasibility_certificate": satisfiable is False,
        "feasible": satisfiable is True,
        "solve_seconds": time.time() - started - build_seconds,
        "elapsed_seconds": time.time() - started,
    }
    if satisfiable is True:
        chosen_edges = [
            list(edge)
            for edge, literal in zip(edges, edge_literals)
            if assignment[literal]
        ]
        output["chosen_edges_H"] = chosen_edges
        output["witness_audit"] = audit_witness(data, c_h, boundary, chosen_edges)
        if not output["witness_audit"]["valid"]:
            raise AssertionError("CryptoMiniSat witness failed independent audit")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--card-encoding",
        choices=tuple(CARDINALITY_ENCODINGS),
        default="seqcounter",
    )
    parser.add_argument("--boundary-internal-edges", type=int)
    parser.add_argument("--boundary-cross-edges", type=int)
    parser.add_argument("--fixed-internal-edge-indices", type=int, nargs="+")
    parser.add_argument("--boundary-cross-degrees", type=int, nargs="+")
    parser.add_argument("--outside-cross-odd-vertices", type=int)
    parser.add_argument("--enumerate-crossing-patterns", action="store_true")
    parser.add_argument("--crossing-start", type=int, default=0)
    parser.add_argument("--crossing-stop", type=int)
    parser.add_argument("--seconds-per-crossing", type=float, default=1.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = solve_case(
        args.source,
        args.orbit_index,
        args.seconds,
        args.threads,
        args.card_encoding,
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
        args.enumerate_crossing_patterns,
        args.crossing_start,
        args.crossing_stop,
        args.seconds_per_crossing,
    )
    atomic_write(args.output, result)
    print(json.dumps({key: value for key, value in result.items() if key != "chosen_edges_H"}, indent=2))


if __name__ == "__main__":
    main()
