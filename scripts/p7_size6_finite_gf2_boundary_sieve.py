#!/usr/bin/env python3
"""Exact GF(2) boundary-parity sieve for p=7 six-finite orbits.

The common binary edge system contains edge count, distinguished edge,
all 50 degree parities, Paley-product parity, and all 280 affine bad-edge
counts.  For a fixed boundary the score slack satisfies

    A_X = |X cap B_d| + phase_d  (mod 2),
    bad_X = 13 - A_X             (mod 2).

Thus the complete right-hand side modulo two is fixed before choosing any
directional scaled means or high-mean slack catalogs.  A nonzero left-null
syndrome rigorously excludes the boundary.  A zero syndrome is only a
necessary condition and is retained.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import field_direction_data, projective_directions  # noqa: E402
from p7_no_infinity_unsaturated_cpsat import atomic_write  # noqa: E402
from p7_unsaturated_gf2_catalog_filter import equation_matrix, row_dependencies  # noqa: E402


POINTS = tuple(itertools.combinations(range(7), 4))


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rhs_mask(c_h: int, boundary: tuple[int, ...]) -> int:
    mask = 0

    def set_rhs(index: int, value: int) -> None:
        nonlocal mask
        if value & 1:
            mask |= 1 << index

    set_rhs(0, 29)
    set_rhs(1, 1)
    boundary_set = set(boundary)
    for vertex in range(50):
        set_rhs(2 + vertex, int(vertex in boundary_set))
    set_rhs(52, int(c_h == -1))
    for direction_index, direction in enumerate(projective_directions(7)):
        eps, labels = field_direction_data(7, direction)
        counts = [0] * 7
        for vertex in boundary:
            counts[labels[vertex - 1]] += 1
        B = {index for index, count in enumerate(counts) if count & 1}
        phase = int(-eps * c_h == -1)
        for point_index, point in enumerate(POINTS):
            slack_parity = (sum(value in B for value in point) + phase) & 1
            bad_parity = 1 ^ slack_parity
            set_rhs(53 + 35 * direction_index + point_index, bad_parity)
    return mask


def run(source: Path, c_h: int) -> dict:
    started = time.time()
    payload = json.loads(source.read_text())
    if (
        int(payload.get("p", 0)) != 7
        or int(payload.get("c_H", 0)) != c_h
        or int(payload.get("boundary_size", 0)) != 6
        or int(payload.get("infinity_value", -1)) != 0
        or payload.get("status") != "complete_boundary_only_enumeration"
    ):
        raise ValueError("source is not the complete p=7 six-finite orbit quotient")
    rows, metadata, variables = equation_matrix()
    rank, dependencies = row_dependencies(rows)
    if rank + len(dependencies) != len(rows):
        raise AssertionError("GF(2) dependency dimensions are inconsistent")
    for dependency in dependencies:
        combination = 0
        bits = dependency
        while bits:
            least = bits & -bits
            combination ^= rows[least.bit_length() - 1]
            bits ^= least
        if combination:
            raise AssertionError("recorded GF(2) witness is not left-null")

    survivors = []
    rejected_orbits = 0
    rejected_boundary_size_sum = 0
    deep_total = 0
    deep_survivors = 0
    for orbit_index, orbit in enumerate(payload["orbits"]):
        boundary = tuple(int(value) for value in orbit["representative_vertices"])
        rhs = rhs_mask(c_h, boundary)
        syndrome = sum(
            (((dependency & rhs).bit_count() & 1) << index)
            for index, dependency in enumerate(dependencies)
        )
        costs = {int(key): int(value) for key, value in orbit["type_costs"].items()}
        deep = any(value not in (24, 32) for value in costs.values())
        deep_total += int(deep)
        if syndrome:
            rejected_orbits += 1
            rejected_boundary_size_sum += int(orbit["size"])
        else:
            deep_survivors += int(deep)
            survivors.append(
                {
                    "orbit_index": orbit_index,
                    "orbit_size": int(orbit["size"]),
                    "representative_vertices": list(boundary),
                    "type_floor_sums": {str(key): value for key, value in costs.items()},
                }
            )
    return {
        "experiment": "p7_size6_finite_gf2_boundary_sieve",
        "status": "complete_exact_gf2_boundary_rhs_sieve",
        "p": 7,
        "c_H": c_h,
        "source": str(source),
        "source_sha256": source_hash(source),
        "source_orbits": int(payload["orbit_count"]),
        "source_boundary_size_sum": int(payload["orbit_size_sum"]),
        "gf2_variables": variables,
        "gf2_equations": len(rows),
        "gf2_rank": rank,
        "gf2_dependency_dimension": len(dependencies),
        "left_null_audit": True,
        "metadata_kinds": {
            kind: sum(row["kind"] == kind for row in metadata)
            for kind in sorted({row["kind"] for row in metadata})
        },
        "rejected_orbits": rejected_orbits,
        "rejected_boundary_size_sum": rejected_boundary_size_sum,
        "surviving_orbits": len(survivors),
        "deep_deficit_orbits": deep_total,
        "deep_deficit_surviving_orbits": deep_survivors,
        "survivors": survivors,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), default=-1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = run(args.source, args.c_h)
    if args.output is not None:
        atomic_write(args.output, out)
    print(json.dumps({key: value for key, value in out.items() if key != "survivors"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
