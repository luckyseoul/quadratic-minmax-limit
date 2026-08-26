#!/usr/bin/env python3
"""Independent audit of the p=7 unsaturated mod-seven certificate.

This audit intentionally differs from ``p7_unsaturated_mod7_batch.py``:

* it builds score rows by vectorized sign products;
* it discovers left dependencies by incremental row-span reduction;
* it obtains catalog right sides from interpolated target coefficients,
  rather than relabelling slack-value vectors;
* it independently reconstructs all boundary/elevation coverage keys.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from p7_no_infinity_unsaturated_cpsat import (  # noqa: E402
    atomic_write,
    direction_target_options,
)
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


POINTS = tuple(itertools.combinations(range(7), 4))
PAIRS = tuple(itertools.combinations(range(7), 2))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_equation_matrix() -> np.ndarray:
    data = geometry(7, "affine")
    C = data["C"]
    edges = data["edges"]
    left = np.asarray([a for a, _b in edges], dtype=np.int16)
    right = np.asarray([b for _a, b in edges], dtype=np.int16)
    rows = [np.ones(len(edges), dtype=np.int16)]
    fixed = np.zeros(len(edges), dtype=np.int16)
    fixed[edges.index((0, 1))] = 1
    rows.append(fixed)
    for direction in projective_directions(7):
        eps, labels = field_direction_data(7, direction)
        for X in POINTS:
            selected = np.zeros(7, dtype=np.int8)
            selected[list(X)] = 1
            y = np.empty(50, dtype=np.int8)
            y[0] = eps
            y[1:] = 2 * selected[np.asarray(labels, dtype=np.int16)] - 1
            feature = eps * y[left] * y[right] * C[left, right]
            rows.append((feature < 0).astype(np.int16))
    matrix = np.stack(rows)
    if matrix.shape != (282, 1225):
        raise AssertionError(f"unexpected matrix shape {matrix.shape}")
    return matrix


def incremental_dependencies(
    matrix: np.ndarray, modulus: int = 7
) -> tuple[int, np.ndarray]:
    """Find row dependencies by online span reduction with witnesses."""
    basis: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    dependencies = []
    m = matrix.shape[0]
    for row_index, original in enumerate(matrix.astype(np.int64) % modulus):
        row = original.copy()
        witness = np.zeros(m, dtype=np.int64)
        witness[row_index] = 1
        while np.any(row):
            pivot = int(np.flatnonzero(row)[0])
            if pivot not in basis:
                inverse = pow(int(row[pivot]), -1, modulus)
                basis[pivot] = (
                    row * inverse % modulus,
                    witness * inverse % modulus,
                )
                break
            basis_row, basis_witness = basis[pivot]
            factor = int(row[pivot])
            row = (row - factor * basis_row) % modulus
            witness = (witness - factor * basis_witness) % modulus
        if not np.any(row):
            dependencies.append(witness)
    return len(basis), np.stack(dependencies)


def elevation_cases(orbit: dict) -> tuple[tuple[int, ...], ...]:
    costs = {int(key): int(value) for key, value in orbit["type_costs"].items()}
    choices = []
    for eps in (-1, 1):
        if costs[eps] == 24:
            directions = tuple(
                index
                for index, row in enumerate(orbit["direction_rows"])
                if int(row["eps"]) == eps
            )
            if len(directions) != 4:
                raise AssertionError("quadratic type does not have four directions")
            choices.append(directions)
        elif costs[eps] != 32:
            raise AssertionError("unexpected type cost")
    return tuple(tuple(sorted(case)) for case in itertools.product(*choices))


def direction_rows(c_h: int, boundary: tuple[int, ...]) -> tuple[list[dict], dict[int, int]]:
    rows = []
    type_floors = {-1: 0, 1: 0}
    for direction in projective_directions(7):
        eps, labels = field_direction_data(7, direction)
        counts = [0] * 7
        for vertex in boundary:
            counts[labels[vertex - 1]] += 1
        B = {index for index, count in enumerate(counts) if count & 1}
        phase = int(-eps * c_h == -1)
        floor = scaled_direction_floor(7, len(B), phase)
        type_floors[eps] += floor
        rows.append(
            {
                "eps": eps,
                "B": B,
                "phase": phase,
                "floor": floor,
            }
        )
    return rows, type_floors


def target_bad_counts(
    odd_fibres: int,
    phase: int,
    B: set[int],
    type_floor_sum: int,
    elevated: bool,
) -> np.ndarray:
    options = direction_target_options(
        odd_fibres, phase, B, type_floor_sum, elevated
    )
    coefficients = np.asarray(
        [[int(value) for value in option[1:]] for option in options],
        dtype=np.int64,
    )
    evaluation = np.empty((35, 22), dtype=np.int64)
    evaluation[:, 0] = 1
    for point_index, X in enumerate(POINTS):
        X_set = set(X)
        evaluation[point_index, 1:] = [
            1 if ((s in X_set) == (t in X_set)) else -1 for s, t in PAIRS
        ]
    scores = evaluation @ coefficients.T
    if np.any(scores < 3) or np.any(scores > 29) or np.any(scores % 2 == 0):
        raise AssertionError("interpolated catalog has an invalid score")
    return ((29 - scores) // 2).astype(np.int16)


def count_matches(base: np.ndarray, contributions: list[np.ndarray]) -> int:
    target = (-base.astype(np.int16)) % 7
    if len(contributions) == 1:
        return int(np.count_nonzero(np.all(contributions[0] == target[:, None], axis=0)))
    if len(contributions) != 2:
        raise AssertionError("expected one or two non-singleton catalogs")
    first, second = contributions
    counts = Counter(bytes(first[:, index]) for index in range(first.shape[1]))
    return sum(
        counts.get(
            bytes(((target - second[:, index].astype(np.int16)) % 7).astype(np.uint8)),
            0,
        )
        for index in range(second.shape[1])
    )


def audit(source: Path, certificate: Path) -> dict:
    started = time.time()
    source_data = json.loads(source.read_text())
    cert = json.loads(certificate.read_text())
    if cert.get("source_sha256") != sha256(source):
        raise AssertionError("certificate source hash mismatch")
    c_h = int(cert["c_H"])
    matrix = independent_equation_matrix()
    rank, dependencies = incremental_dependencies(matrix, 7)
    if rank != 147 or dependencies.shape != (135, 282):
        raise AssertionError("independent mod-seven dimensions disagree")
    if np.any(dependencies @ (matrix % 7) % 7):
        raise AssertionError("independent dependency matrix is invalid")

    cert_rows = {
        (
            int(row["orbit_index"]),
            tuple(int(value) for value in row["elevated_directions"]),
        ): row
        for row in cert["rows"]
    }
    if len(cert_rows) != len(cert["rows"]):
        raise AssertionError("duplicate certificate case")

    base_edge = (
        dependencies[:, :2] @ np.asarray([29, 1], dtype=np.int64) % 7
    ).astype(np.uint8)
    target_cache: dict[tuple, np.ndarray] = {}
    contribution_cache: dict[tuple, np.ndarray] = {}
    expected_keys = set()
    recomputed_nonzero = []
    pattern_counts: Counter[tuple[int, ...]] = Counter()
    unsaturated_orbits = 0
    boundary_size_sum = 0

    for orbit_index, orbit in enumerate(source_data["orbits"]):
        if all(int(value) == 32 for value in orbit["type_costs"].values()):
            continue
        unsaturated_orbits += 1
        boundary_size_sum += int(orbit["size"])
        boundary = tuple(int(value) for value in orbit["representative_vertices"])
        rows, type_floors = direction_rows(c_h, boundary)
        for elevated in elevation_cases(orbit):
            key = (orbit_index, elevated)
            expected_keys.add(key)
            base = base_edge.astype(np.int16).copy()
            variable = []
            for direction_index, row in enumerate(rows):
                target_key = (
                    len(row["B"]),
                    int(row["phase"]),
                    tuple(sorted(row["B"])),
                    int(type_floors[int(row["eps"])]),
                    direction_index in set(elevated),
                )
                if target_key not in target_cache:
                    target_cache[target_key] = target_bad_counts(
                        target_key[0],
                        target_key[1],
                        set(target_key[2]),
                        target_key[3],
                        target_key[4],
                    )
                contribution_key = (direction_index, target_key)
                if contribution_key not in contribution_cache:
                    block = dependencies[
                        :,
                        2 + 35 * direction_index : 2 + 35 * (direction_index + 1),
                    ]
                    contribution_cache[contribution_key] = (
                        block @ target_cache[target_key] % 7
                    ).astype(np.uint8)
                contribution = contribution_cache[contribution_key]
                if contribution.shape[1] == 1:
                    base = (base + contribution[:, 0]) % 7
                else:
                    variable.append(contribution)
            pattern = tuple(sorted((item.shape[1] for item in variable), reverse=True))
            pattern_counts[pattern] += 1
            matches = count_matches(base.astype(np.uint8), variable)
            if matches:
                recomputed_nonzero.append(
                    {"orbit_index": orbit_index, "elevated": list(elevated), "count": matches}
                )
            recorded = cert_rows.get(key)
            if recorded is None:
                raise AssertionError(f"missing certificate case {key}")
            if int(recorded["mod7_consistent_catalog_tuples"]) != matches:
                raise AssertionError(f"syndrome count mismatch for {key}")

    missing = sorted(expected_keys - set(cert_rows))
    extra = sorted(set(cert_rows) - expected_keys)
    passed = bool(
        not missing
        and not extra
        and not recomputed_nonzero
        and unsaturated_orbits == 518
        and boundary_size_sum == 23520
        and len(expected_keys) == 2408
    )
    return {
        "experiment": "p7_unsaturated_mod7_audit",
        "status": "independent_target_coefficient_and_row_span_audit",
        "passed": passed,
        "source": str(source),
        "source_sha256": sha256(source),
        "certificate": str(certificate),
        "certificate_sha256": sha256(certificate),
        "modulus": 7,
        "independent_rank": rank,
        "independent_dependency_dimension": int(len(dependencies)),
        "left_null_audit": True,
        "unsaturated_orbits": unsaturated_orbits,
        "unsaturated_boundary_size_sum": boundary_size_sum,
        "expected_cases": len(expected_keys),
        "certificate_cases": len(cert_rows),
        "missing_cases": len(missing),
        "extra_cases": len(extra),
        "recomputed_nonzero_cases": len(recomputed_nonzero),
        "catalog_pattern_counts": {
            "x".join(map(str, pattern)): count
            for pattern, count in sorted(pattern_counts.items())
        },
        "target_cache_entries": len(target_cache),
        "contribution_cache_entries": len(contribution_cache),
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = audit(args.source, args.certificate)
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
