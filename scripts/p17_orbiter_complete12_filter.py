#!/usr/bin/env python3
"""Independently identify complete 12-arcs in Orbiter's PG(2,17) census."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import time
from itertools import combinations
from pathlib import Path


P = 17
N = P * P + P + 1


def atomic_write(path: Path, payload: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def orbiter_unrank(rank: int, q: int = P, length: int = 3) -> tuple[int, ...]:
    """Port of PG_element_unrank_modified_lint from Orbiter build 3361."""
    original = rank
    if rank < length:
        return tuple(int(i == rank) for i in range(length))
    rank -= length
    if rank == 0:
        return (1,) * length
    rank -= 1
    level = 1
    q_power = q
    geometric_sum = 1
    while level < length:
        if rank >= q_power - 1:
            rank -= q_power - 1
            geometric_sum += q_power
            q_power *= q
            level += 1
            continue
        vector = [0] * length
        vector[level] = 1
        rank += 1
        if level == length - 1 and rank >= geometric_sum:
            rank += 1
        index = 0
        while rank:
            remainder = rank % q
            vector[index] = remainder
            index += 1
            rank = (rank - remainder) // q
        return tuple(vector)
    raise ValueError(f"rank {original} is outside PG({length - 1},{q})")


def determinant(a: tuple[int, ...], b: tuple[int, ...], c: tuple[int, ...]) -> int:
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    ) % P


def read_representatives(path: Path) -> list[dict[str, object]]:
    rows = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(line for line in stream if line.strip() != "END")
        for row in reader:
            representative = tuple(int(value) for value in row["OrbitRep"].split(","))
            rows.append(
                {
                    "orbit_index": int(row["OrbitIdx"]),
                    "node": int(row["Node"]),
                    "representative": representative,
                    "stabilizer_order": int(row["StabOrder"]),
                    "orbit_length": int(row["OrbitLength"]),
                }
            )
    return rows


def addable_points(
    representative: tuple[int, ...],
    points: tuple[tuple[int, ...], ...],
) -> tuple[int, ...]:
    chosen = set(representative)
    blocked = set(chosen)
    for left, right in combinations(representative, 2):
        a, b = points[left], points[right]
        blocked.update(
            rank
            for rank, point in enumerate(points)
            if determinant(a, b, point) == 0
        )
    return tuple(rank for rank in range(N) if rank not in blocked)


def classify(path: Path) -> dict[str, object]:
    points = tuple(orbiter_unrank(rank) for rank in range(N))
    if len(set(points)) != N:
        raise ArithmeticError("Orbiter point unranking is not injective")
    rows = read_representatives(path)
    if len(rows) != 629:
        raise ArithmeticError("Orbiter p17 12-arc orbit count changed")

    complete = []
    extendible = []
    for row in rows:
        representative = row["representative"]
        if len(representative) != 12 or len(set(representative)) != 12:
            raise ArithmeticError("malformed 12-arc representative")
        if any(
            determinant(points[a], points[b], points[c]) == 0
            for a, b, c in combinations(representative, 3)
        ):
            raise ArithmeticError("Orbiter representative is not a 12-arc")
        addable = addable_points(representative, points)
        record = {
            **row,
            "representative": list(representative),
            "addable_points": list(addable),
            "complete": not addable,
        }
        (complete if not addable else extendible).append(record)

    if len(complete) != 553 or len(extendible) != 76:
        raise ArithmeticError(
            f"complete/extendible split changed: {len(complete)}/{len(extendible)}"
        )
    return {
        "experiment": "p17_orbiter_complete12_independent_filter",
        "orbiter_build": 3361,
        "p": P,
        "projective_point_count": N,
        "source_csv": str(path.resolve()),
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "twelve_arc_orbit_count": len(rows),
        "complete_twelve_arc_orbit_count": len(complete),
        "extendible_twelve_arc_orbit_count": len(extendible),
        "complete_representatives": complete,
        "extendible_representatives": extendible,
        "independent_collinearity_audit": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = classify(args.csv)
    if args.output is not None:
        atomic_write(args.output, result)
    print(
        f"12-arcs={result['twelve_arc_orbit_count']} "
        f"complete={result['complete_twelve_arc_orbit_count']} "
        f"extendible={result['extendible_twelve_arc_orbit_count']}"
    )


if __name__ == "__main__":
    main()
