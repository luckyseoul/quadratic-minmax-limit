#!/usr/bin/env python3
"""Exact PGL class generator for complete 13-/14-arcs in PG(2,17).

Every arc contains an ordered quadrangle.  PGL(3,17) acts transitively and
freely on ordered quadrangles, so fix the four canonical points e1,e2,e3 and
(1,1,1).  After each complete arc is found, enumerate all normalized images
obtained by mapping every ordered quadrangle of that arc to the canonical
one, deduplicate them, and block every image.  An eventual INFEASIBLE result
therefore proves that the accumulated representatives exhaust the PGL
classes within the requested size.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from collections import Counter
from itertools import permutations
from pathlib import Path


P = 17
Point = tuple[int, int, int]
CANONICAL_QUADRANGLE: tuple[Point, ...] = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 1),
)


def projective_points_or_lines() -> tuple[Point, ...]:
    return tuple(
        [(1, y, z) for y in range(P) for z in range(P)]
        + [(0, 1, z) for z in range(P)]
        + [(0, 0, 1)]
    )


def incident(point: Point, line: Point) -> bool:
    return sum(a * b for a, b in zip(point, line)) % P == 0


def normalize(vector: tuple[int, int, int]) -> Point:
    for value in vector:
        if value % P:
            inverse = pow(value % P, -1, P)
            return tuple(entry * inverse % P for entry in vector)  # type: ignore[return-value]
    raise ValueError("zero vector has no projective normalization")


def matrix_inverse(matrix: tuple[Point, Point, Point]) -> tuple[Point, Point, Point]:
    a, b, c = matrix
    determinant = (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    ) % P
    if determinant == 0:
        raise ValueError("singular matrix")
    scale = pow(determinant, -1, P)
    cofactors = (
        (
            b[1] * c[2] - b[2] * c[1],
            a[2] * c[1] - a[1] * c[2],
            a[1] * b[2] - a[2] * b[1],
        ),
        (
            b[2] * c[0] - b[0] * c[2],
            a[0] * c[2] - a[2] * c[0],
            a[2] * b[0] - a[0] * b[2],
        ),
        (
            b[0] * c[1] - b[1] * c[0],
            a[1] * c[0] - a[0] * c[1],
            a[0] * b[1] - a[1] * b[0],
        ),
    )
    return tuple(
        tuple(value * scale % P for value in row) for row in cofactors
    )  # type: ignore[return-value]


def row_times_matrix(row: Point, matrix: tuple[Point, Point, Point]) -> Point:
    return tuple(
        sum(row[k] * matrix[k][column] for k in range(3)) % P
        for column in range(3)
    )  # type: ignore[return-value]


def normalized_orbit_images(arc: tuple[Point, ...]) -> set[tuple[Point, ...]]:
    """All images of one PGL orbit containing the canonical quadrangle."""
    images: set[tuple[Point, ...]] = set()
    for first, second, third, fourth in permutations(arc, 4):
        inverse = matrix_inverse((first, second, third))
        coordinates = row_times_matrix(fourth, inverse)
        if any(value == 0 for value in coordinates):
            raise ArithmeticError("four arc points failed the quadrangle test")
        column_scales = tuple(pow(value, -1, P) for value in coordinates)
        transformation = tuple(
            tuple(inverse[row][column] * column_scales[column] % P for column in range(3))
            for row in range(3)
        )
        image = tuple(sorted(normalize(row_times_matrix(point, transformation)) for point in arc))
        if not set(CANONICAL_QUADRANGLE) <= set(image):
            raise ArithmeticError("normalized orbit image lost the quadrangle")
        images.add(image)
    return images


def secant_index_histogram(arc: tuple[Point, ...]) -> dict[int, int]:
    points = projective_points_or_lines()
    selected = set(arc)
    secants = [
        line
        for line in points
        if sum(point in selected for point in points if incident(point, line)) == 2
    ]
    return dict(
        sorted(
            Counter(
                sum(incident(point, line) for line in secants)
                for point in points
                if point not in selected
            ).items()
        )
    )


def classify(
    size: int,
    seconds: float,
    workers: int,
    use_published_class_count: bool = False,
) -> dict[str, object]:
    from ortools.sat.python import cp_model

    if size not in (13, 14):
        raise ValueError("size must be 13 or 14")
    points = projective_points_or_lines()
    index = {point: i for i, point in enumerate(points)}
    lines = points
    line_points = [
        [i for i, point in enumerate(points) if incident(point, line)]
        for line in lines
    ]
    point_lines = [
        [j for j, row in enumerate(line_points) if i in row]
        for i in range(len(points))
    ]
    model = cp_model.CpModel()
    chosen = [model.new_bool_var(f"x_{i}") for i in range(len(points))]
    model.add(sum(chosen) == size)
    for point in CANONICAL_QUADRANGLE:
        model.add(chosen[index[point]] == 1)
    secants = []
    for j, row in enumerate(line_points):
        occupancy = sum(chosen[i] for i in row)
        model.add(occupancy <= 2)
        secant = model.new_bool_var(f"secant_{j}")
        model.add(occupancy == 2).only_enforce_if(secant)
        model.add(occupancy <= 1).only_enforce_if(~secant)
        secants.append(secant)
    for i in range(len(points)):
        model.add(chosen[i] + sum(secants[j] for j in point_lines[i]) >= 1)

    expected_classes = {13: 8, 14: 1}[size]
    representatives = []
    solve_rows = []
    started = time.time()
    terminal_status = "PUBLISHED_CLASS_COUNT_REACHED"
    solve_count = expected_classes if use_published_class_count else expected_classes + 1
    for class_index in range(solve_count):
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = seconds
        solver.parameters.num_search_workers = workers
        solver.parameters.random_seed = 15702000 + 100 * size + class_index
        solver.parameters.symmetry_level = 2
        status = solver.solve(model)
        solve_rows.append(
            {
                "class_index": class_index,
                "solver_status": solver.status_name(status),
                "wall_time_seconds": solver.wall_time,
                "branches": solver.num_branches,
                "conflicts": solver.num_conflicts,
            }
        )
        if status == cp_model.INFEASIBLE:
            break
        if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            raise RuntimeError("class generation did not finish exactly")
        arc = tuple(points[i] for i, variable in enumerate(chosen) if solver.value(variable))
        images = normalized_orbit_images(arc)
        ordered_quadrangles = size * (size - 1) * (size - 2) * (size - 3)
        if ordered_quadrangles % len(images):
            raise ArithmeticError("normalized orbit size does not divide quadrangle count")
        stabilizer_order = ordered_quadrangles // len(images)
        representatives.append(
            {
                "representative": [list(point) for point in arc],
                "normalized_orbit_image_count": len(images),
                "inferred_pgl_stabilizer_order": stabilizer_order,
                "outside_secant_index_histogram": secant_index_histogram(arc),
            }
        )
        print(
            f"size={size} class={class_index + 1}/{expected_classes} "
            f"orbit_images={len(images)} stabilizer={stabilizer_order}",
            flush=True,
        )
        for image in images:
            model.add(sum(chosen[index[point]] for point in image) <= size - 1)

    if not use_published_class_count:
        terminal_status = str(solve_rows[-1]["solver_status"])
    if len(representatives) != expected_classes or (
        not use_published_class_count and terminal_status != "INFEASIBLE"
    ):
        raise ArithmeticError("complete-arc class exhaustion changed")
    observed_orders = sorted(int(row["inferred_pgl_stabilizer_order"]) for row in representatives)
    expected_orders = {13: [1, 2, 2, 2, 2, 3, 4, 6], 14: [8]}[size]
    if observed_orders != expected_orders:
        raise ArithmeticError("complete-arc stabilizer fingerprint changed")
    return {
        "experiment": "p17_complete_arc_class_generator",
        "p": P,
        "arc_size": size,
        "normalization": [list(point) for point in CANONICAL_QUADRANGLE],
        "class_count": len(representatives),
        "representatives": representatives,
        "solve_rows": solve_rows,
        "final_status": terminal_status,
        "classification_basis": (
            "inequivalent representatives match Sticker's published class "
            "count and stabilizer-order fingerprint"
            if use_published_class_count
            else "normalized PGL orbit blocking followed by exact INFEASIBLE"
        ),
        "elapsed_seconds": time.time() - started,
        "proved_by_local_terminal_infeasibility": not use_published_class_count,
        "proved_conditional_on_published_class_count": use_published_class_count,
    }


def atomic_write(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--size", type=int, choices=(13, 14), required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument(
        "--use-published-class-count",
        action="store_true",
        help="stop after the published number of inequivalent classes",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = classify(
        args.size,
        args.seconds,
        args.workers,
        args.use_published_class_count,
    )
    if args.output is not None:
        atomic_write(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
