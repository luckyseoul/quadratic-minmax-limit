#!/usr/bin/env python3
"""Exact CPU frontier scorer for the p=31 center-coordinate search.

This is an independent replay layer for candidates produced by the batched
center search.  It accepts one nonzero center for each of the sixteen frozen
scaled ``(L,M)`` pairs, accumulates all 480 inversion-orbit occurrences, and
requires the resulting sum to be ternary with exactly one cancellation unit.
Both a two-half ``1:1`` cancellation and a clean three-half ``2:1``
cancellation are accepted; no pair-only reduction is used.

For every admissible candidate the physical 479-edge graph is rebuilt, all
32 normalized transverse rows are computed with the hard literal stars
removed, and :class:`IncrementalAtomRowBound` supplies the exact singleton/
two-label-cut search cost and a rigorous coefficient-edit lower bound.  A
zero score is necessary, not sufficient, for the remaining integral atom
transport problem.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import combinations
import hashlib
import json
from pathlib import Path
import sys
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_inversion_antisymmetric_radon import (  # noqa: E402
    Edge,
    edge_radon_image,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (  # noqa: E402
    paley_direction_sign,
    paley_edge_sign,
)
from e1_gmin_m4_p31_direct_mobius_parallel_design import (  # noqa: E402
    HALVES,
    P,
    PHYSICAL_CENTERS,
    PHYSICAL_FIXED_DIRECTION_INDEX,
    PHYSICAL_FIXED_POINT,
    _canonical_center,
    _negative_edge,
    _oriented_orbit_coefficients,
    _spatial_direction_index,
)
from e1_gmin_m4_p31_row_atom_incremental_bound import (  # noqa: E402
    AtomRowSpec,
    IncrementalAtomRowBound,
)
from io_atomic import write_json_atomic  # noqa: E402


def _fixed_edge(point: tuple[int, int]) -> Edge:
    point = point[0] % P, point[1] % P
    if point == (0, 0):
        raise ValueError("the fixed antipodal edge needs a nonzero point")
    negative = -point[0] % P, -point[1] % P
    return tuple(sorted((point, negative)))  # type: ignore[return-value]


def fixed_edges_in_target_direction() -> tuple[Edge, ...]:
    """Return all fifteen antipodal edges in the frozen fixed direction."""
    edges = {
        _fixed_edge(
            (
                scalar * PHYSICAL_FIXED_POINT[0] % P,
                scalar * PHYSICAL_FIXED_POINT[1] % P,
            )
        )
        for scalar in range(1, P)
    }
    out = tuple(sorted(edges))
    if (
        len(out) != (P - 1) // 2
        or any(
            _spatial_direction_index(edge) != PHYSICAL_FIXED_DIRECTION_INDEX
            for edge in out
        )
    ):
        raise ArithmeticError("the fixed-edge magnitude orbit changed")
    return out


def _center_orbit_sum(
    centers: Sequence[int],
) -> tuple[Counter[Edge], dict[Edge, list[tuple[int, int]]]]:
    if len(centers) != len(HALVES):
        raise ValueError(f"need exactly {len(HALVES)} centers")
    if any(not 1 <= int(center) < P for center in centers):
        raise ValueError("every displayed center must lie in F_31^*")
    total: Counter[Edge] = Counter()
    occurrences: dict[Edge, list[tuple[int, int]]] = defaultdict(list)
    for half_index, ((target, auxiliary), center) in enumerate(
        zip(HALVES, centers, strict=True)
    ):
        orbit_map = _oriented_orbit_coefficients(
            target, auxiliary, int(center)
        )
        if len(orbit_map) != P - 1:
            raise ArithmeticError("a center configuration lost a half orbit")
        total.update(orbit_map)
        for orbit, coefficient in orbit_map.items():
            occurrences[orbit].append((half_index, coefficient))
    return total, occurrences


def classify_center_geometry(centers: Sequence[int]) -> dict[str, object]:
    """Classify ternarity and cancellation directly from orbit occurrences."""
    total, occurrences = _center_orbit_sum(centers)
    raw_occurrences = len(HALVES) * (P - 1)
    absolute_support_mass = sum(abs(value) for value in total.values())
    deficit = raw_occurrences - absolute_support_mass
    cancellation_units = deficit // 2 if deficit >= 0 and deficit % 2 == 0 else -1
    maximum_coefficient = max((abs(value) for value in total.values()), default=0)
    shared = {
        orbit: rows for orbit, rows in occurrences.items() if len(rows) > 1
    }
    nonzero_shared = {
        orbit: total[orbit] for orbit in shared if total[orbit]
    }
    ternary = maximum_coefficient <= 1
    admissible = bool(ternary and cancellation_units == 1)
    cancellation_orbits = tuple(
        orbit
        for orbit, rows in shared.items()
        if (len(rows) - abs(total[orbit])) // 2 > 0
    )
    if admissible and len(cancellation_orbits) != 1:
        raise ArithmeticError("cancellation-one ternarity lost its unique orbit")
    collision_orbit = cancellation_orbits[0] if admissible else None
    collision_rows = shared.get(collision_orbit, []) if collision_orbit else []
    collision_direction = (
        _spatial_direction_index(collision_orbit)
        if collision_orbit is not None
        else None
    )
    return {
        "raw_orbit_occurrences": raw_occurrences,
        "distinct_raw_orbit_count": len(total),
        "absolute_final_support_mass": absolute_support_mass,
        "maximum_summed_orbit_coefficient": maximum_coefficient,
        "ternary": ternary,
        "cancellation_units": cancellation_units,
        "shared_orbit_count": len(shared),
        "nonzero_shared_orbit_count": len(nonzero_shared),
        "unique_cancellation_orbit": collision_orbit,
        "cancellation_multiplicity": len(collision_rows),
        "cancellation_half_indices": tuple(row[0] for row in collision_rows),
        "cancellation_orientation_coefficients": tuple(
            row[1] for row in collision_rows
        ),
        "cancellation_direction_index": collision_direction,
        "fixed_direction_compatible": (
            collision_direction == PHYSICAL_FIXED_DIRECTION_INDEX
        ),
        "admissible_top_geometry": bool(
            admissible
            and collision_direction == PHYSICAL_FIXED_DIRECTION_INDEX
        ),
    }


def _physical_graph(centers: Sequence[int], fixed_edge: Edge) -> tuple[Edge, ...]:
    total, _occurrences = _center_orbit_sum(centers)
    graph = [
        orbit if value == 1 else _negative_edge(P, orbit)
        for orbit, value in total.items()
        if value
    ]
    if any(abs(value) > 1 for value in total.values()):
        raise ValueError("the center tuple is not ternary")
    if fixed_edge in graph:
        raise ArithmeticError("a nonfixed half produced an antipodal edge")
    graph.append(fixed_edge)
    out = tuple(sorted(graph))
    if len(out) != len(set(out)):
        raise ArithmeticError("the reconstructed graph repeated an edge")
    return out


def _normalized_rows(
    centers: Sequence[int], fixed_edge: Edge
) -> tuple[tuple[Edge, ...], tuple[tuple[int, int, dict[tuple[int, int], int]], ...]]:
    graph = _physical_graph(centers, fixed_edge)
    directions = projective_functionals(P)
    signs = tuple(paley_direction_sign(P, direction) for direction in directions)
    canonical_centers = {
        _canonical_center(target, int(center))[0]: _canonical_center(
            target, int(center)
        )[1]
        for (target, _auxiliary), center in zip(HALVES, centers, strict=True)
    }
    image = edge_radon_image(
        P, {edge: paley_edge_sign(P, edge) for edge in graph}
    )
    rows = []
    for direction_index, sign in enumerate(signs):
        parallel = sign * image.get(("P", direction_index), 0)
        coefficients = {
            (left, right): sign
            * image.get(("K", direction_index, left, right), 0)
            for left, right in combinations(range(P), 2)
        }
        if sign == 1:
            center = canonical_centers[direction_index]
            for other in range(P):
                if other == center:
                    continue
                edge = tuple(sorted((center, other)))
                coefficients[edge] += 1
        coefficients = {
            edge: value for edge, value in coefficients.items() if value
        }
        rows.append((direction_index, parallel, coefficients))
    return graph, tuple(rows)


def score_centers(
    centers: Sequence[int], fixed_point: tuple[int, int] = PHYSICAL_FIXED_POINT
) -> dict[str, object]:
    """Return the exact 32-row necessary atom score for one center tuple."""
    centers = tuple(int(center) for center in centers)
    geometry = classify_center_geometry(centers)
    if not geometry["admissible_top_geometry"]:
        return {
            "centers": centers,
            "fixed_point": fixed_point,
            "geometry": geometry,
            "status": "REJECTED_BY_TERNARITY_OR_TOP_CANCELLATION",
        }
    fixed_edge = _fixed_edge(fixed_point)
    if _spatial_direction_index(fixed_edge) != PHYSICAL_FIXED_DIRECTION_INDEX:
        raise ValueError("fixed point has the wrong spatial direction")
    graph, raw_rows = _normalized_rows(centers, fixed_edge)
    if len(graph) != 479:
        raise ArithmeticError("an admissible top geometry did not yield 479 edges")

    rows = []
    for direction_index, parallel, coefficients in raw_rows:
        sign = paley_direction_sign(
            P, projective_functionals(P)[direction_index]
        )
        spec = (
            AtomRowSpec.hard(parallel - 3)
            if sign == 1
            else AtomRowSpec.opposite(parallel - 9)
        )
        summary = IncrementalAtomRowBound(coefficients, spec).summary()
        summary["direction_index"] = direction_index
        summary["direction_sign"] = sign
        summary["parallel_count"] = parallel
        rows.append(summary)

    graph_bytes = json.dumps(graph, separators=(",", ":")).encode()
    search_costs = tuple(int(row["incremental_search_cost"]) for row in rows)
    edit_bounds = tuple(
        int(row["coefficient_l1_edit_lower_bound"]) for row in rows
    )
    return {
        "centers": centers,
        "fixed_point": fixed_point,
        "fixed_edge": fixed_edge,
        "geometry": geometry,
        "status": "EXACT_NECESSARY_SCORE",
        "graph_edge_count": len(graph),
        "graph_sha256": hashlib.sha256(graph_bytes).hexdigest(),
        "violating_row_count": sum(value > 0 for value in search_costs),
        "total_incremental_search_cost": sum(search_costs),
        "maximum_row_incremental_search_cost": max(search_costs),
        "sum_row_edit_lower_bounds": sum(edit_bounds),
        "maximum_row_edit_lower_bound": max(edit_bounds),
        "all_496_cut_row_bounds_pass": all(value == 0 for value in search_costs),
        "rows": tuple(rows),
        "atom_decomposition_constructed": False,
        "residual_ii_closed": False,
    }


def score_all_fixed_edges(centers: Sequence[int]) -> dict[str, object]:
    """Score all 15 fixed-edge magnitudes and return the exact best one."""
    records = [
        score_centers(centers, fixed_edge[0])
        for fixed_edge in fixed_edges_in_target_direction()
    ]
    admissible = [row for row in records if row["status"] == "EXACT_NECESSARY_SCORE"]
    if not admissible:
        return {
            "centers": tuple(int(value) for value in centers),
            "status": records[0]["status"],
            "geometry": records[0]["geometry"],
            "fixed_edge_count": len(records),
        }
    key = lambda row: (
        int(row["total_incremental_search_cost"]),
        int(row["violating_row_count"]),
        int(row["maximum_row_edit_lower_bound"]),
        int(row["sum_row_edit_lower_bounds"]),
    )
    best = min(admissible, key=key)
    return {
        "centers": tuple(int(value) for value in centers),
        "status": "EXACT_ALL_FIXED_EDGE_FRONTIER",
        "fixed_edge_count": len(records),
        "ranking": [
            "total_incremental_search_cost",
            "violating_row_count",
            "maximum_row_edit_lower_bound",
            "sum_row_edit_lower_bounds",
        ],
        "best": best,
        "records": tuple(records),
    }


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--centers", type=int, nargs=16, default=PHYSICAL_CENTERS)
    parser.add_argument("--fixed-point", type=int, nargs=2)
    parser.add_argument("--all-fixed", action="store_true")
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    if args.all_fixed:
        result = score_all_fixed_edges(args.centers)
    else:
        fixed_point = (
            tuple(args.fixed_point) if args.fixed_point else PHYSICAL_FIXED_POINT
        )
        result = score_centers(args.centers, fixed_point)  # type: ignore[arg-type]
    if args.output:
        write_json_atomic(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == "__main__":
    main()
