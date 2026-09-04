#!/usr/bin/env python3
"""Exact coupled congruence survival at the p=31 hard-fixed endpoint.

This module audits the first integral layer left open by the fractional
five-type semimetric decomposition.  It constructs one actual 479-edge
ledger graph with the exact top parallel profile, one fixed antipodal edge
in hard direction 1, no doubled nonfixed inversion orbit, and a common
choice of hard centres.  Every one of its 32 normalized transverse
residuals belongs to the *signed* integer atom lattice with the required
atom-count coordinates.

The algebra is all-n: compact atoms

    K(a,b;c) = e_ab - e_ac - e_bc

generate exactly the integer edge vectors having even signed degree at
every label.  Modulo two, their supports are triangles and span the cycle
space.  Once parity is removed, the identity

    K(u,w;v) + K(v,w;u) = -2 e_uv

corrects every remaining even cell.  Also

    K(a,c;b) + K(b,c;a) = -2 e_ab = e_ab  (mod 3),

so the mod-three audit has no additional coordinate obstruction.

This is deliberately a lattice/congruence result, not a nonnegative atom
decomposition.  The explicit replay uses negative compact multiplicities;
therefore Chvatal--Gomory cuts that genuinely use atom nonnegativity remain
available, and no residual-(ii) witness is asserted.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, islice
import json
import random
from typing import Iterable, Mapping

from e1_gmin_m4_inversion_antisymmetric_radon import projective_functionals
from e1_gmin_m4_mobius_half_symmetric import paley_edge_sign
from e1_gmin_m4_p31_top_j1_f3_case_split import (
    DIRECTION_SIGNS,
    HARD,
    TARGET_PROFILE,
)


P = 31
FIXED_DIRECTION = 1
CONSTRUCTION_SEED = 7_042_026
DIRECTIONS = tuple(projective_functionals(P))

Point = tuple[int, int]
PhysicalEdge = tuple[Point, Point]
LabelEdge = tuple[int, int]
CompactAtom = tuple[int, int, int]


def _add(first: Point, second: Point) -> Point:
    return (
        (first[0] + second[0]) % P,
        (first[1] + second[1]) % P,
    )


def _scale(scalar: int, point: Point) -> Point:
    return scalar * point[0] % P, scalar * point[1] % P


def _negative(point: Point) -> Point:
    return -point[0] % P, -point[1] % P


def _physical_edge(first: Point, second: Point) -> PhysicalEdge:
    if first == second:
        raise ValueError("a physical edge cannot be a loop")
    return tuple(sorted((first, second)))  # type: ignore[return-value]


def _negative_edge(edge: PhysicalEdge) -> PhysicalEdge:
    return _physical_edge(_negative(edge[0]), _negative(edge[1]))


def _orbit_key(edge: PhysicalEdge) -> PhysicalEdge:
    return min(edge, _negative_edge(edge))


def _is_fixed(edge: PhysicalEdge) -> bool:
    return edge[1] == _negative(edge[0])


def _functional_value(functional: Point, point: Point) -> int:
    return (
        functional[0] * point[0] + functional[1] * point[1]
    ) % P


def _spatial_direction_index(edge: PhysicalEdge) -> int:
    difference = (
        (edge[0][0] - edge[1][0]) % P,
        (edge[0][1] - edge[1][1]) % P,
    )
    hits = tuple(
        index
        for index, functional in enumerate(DIRECTIONS)
        if _functional_value(functional, difference) == 0
    )
    if len(hits) != 1:
        raise ArithmeticError("a physical edge lost its spatial direction")
    return hits[0]


def _kernel_vector(direction_index: int) -> Point:
    first, second = DIRECTIONS[direction_index]
    return second % P, -first % P


def _label_edge(first: int, second: int) -> LabelEdge:
    if first == second:
        raise ValueError("a label edge cannot be a loop")
    return tuple(sorted((first, second)))  # type: ignore[return-value]


def _compact_key(first: int, second: int, distinguished: int) -> CompactAtom:
    if len({first, second, distinguished}) != 3:
        raise ValueError("a compact atom needs three distinct labels")
    left, right = sorted((first, second))
    return left, right, distinguished


def compact_atom_vector(atom: CompactAtom) -> dict[LabelEdge, int]:
    """Return ``K(a,b;c)=e_ab-e_ac-e_bc`` in label-edge coordinates."""
    first, second, distinguished = atom
    if len({first, second, distinguished}) != 3:
        raise ValueError("a compact atom needs three distinct labels")
    return {
        _label_edge(first, second): 1,
        _label_edge(first, distinguished): -1,
        _label_edge(second, distinguished): -1,
    }


def _add_scaled(
    target: Counter[LabelEdge],
    source: Mapping[LabelEdge, int],
    multiplier: int,
) -> None:
    for edge, coefficient in source.items():
        target[edge] += multiplier * coefficient
        if not target[edge]:
            del target[edge]


def signed_compact_lattice_decomposition(
    coefficients: Mapping[LabelEdge, int],
    label_count: int = P,
) -> Counter[CompactAtom]:
    """Exactly decompose any integer edge vector with even signed degrees.

    Coefficients of the returned compact atoms are signed integers.  This is
    a constructive sufficiency proof for the compact-atom lattice, not a
    nonnegative-semigroup routine.
    """
    if label_count < 3:
        raise ValueError("at least three labels are required")
    cleaned = {
        _label_edge(*edge): int(value)
        for edge, value in coefficients.items()
        if value
    }
    if any(
        not 0 <= label < label_count
        for edge in cleaned
        for label in edge
    ):
        raise ValueError("a label is out of range")

    degrees = [0] * label_count
    for (first, second), coefficient in cleaned.items():
        degrees[first] += coefficient
        degrees[second] += coefficient
    if any(degree % 2 for degree in degrees):
        raise ValueError("the compact lattice requires even signed degrees")

    atoms: Counter[CompactAtom] = Counter()
    rebuilt: Counter[LabelEdge] = Counter()

    # Triangles through root zero generate the complete-graph cycle space
    # modulo two.  Even degree makes the remaining root edges vanish too.
    for (first, second), coefficient in sorted(cleaned.items()):
        if first != 0 and second != 0 and coefficient % 2:
            atom = _compact_key(first, second, 0)
            atoms[atom] += 1
            _add_scaled(rebuilt, compact_atom_vector(atom), 1)

    residual: Counter[LabelEdge] = Counter(cleaned)
    _add_scaled(residual, rebuilt, -1)
    if any(coefficient % 2 for coefficient in residual.values()):
        raise ArithmeticError("cycle-space reduction left an odd cell")

    # The two-atom identity independently supplies twice any label edge.
    for (first, second), coefficient in sorted(residual.items()):
        if not coefficient:
            continue
        third = next(
            label
            for label in range(label_count)
            if label not in (first, second)
        )
        multiplier = -(coefficient // 2)
        atoms[_compact_key(first, third, second)] += multiplier
        atoms[_compact_key(second, third, first)] += multiplier

    atoms = Counter(
        {atom: coefficient for atom, coefficient in atoms.items() if coefficient}
    )
    replay: Counter[LabelEdge] = Counter()
    for atom, coefficient in atoms.items():
        _add_scaled(replay, compact_atom_vector(atom), coefficient)
    if dict(replay) != cleaned:
        raise ArithmeticError("signed compact decomposition failed replay")
    if sum(atoms.values()) != -sum(cleaned.values()):
        raise ArithmeticError("the compact atom-count coordinate changed")
    return atoms


class _LedgerGraphBuilder:
    """Small deterministic constructor; it performs no optimization search."""

    def __init__(self) -> None:
        self.edges: set[PhysicalEdge] = set()
        self.orbits: set[PhysicalEdge] = set()
        self.direction_counts = [0] * (P + 1)

    def allowed(
        self, edge: PhysicalEdge, *, allow_fixed: bool = False
    ) -> bool:
        return bool(
            edge not in self.edges
            and _orbit_key(edge) not in self.orbits
            and (allow_fixed or not _is_fixed(edge))
        )

    def insert(
        self, edge: PhysicalEdge, *, allow_fixed: bool = False
    ) -> None:
        if not self.allowed(edge, allow_fixed=allow_fixed):
            raise ValueError("attempted to insert a forbidden edge")
        self.edges.add(edge)
        self.orbits.add(_orbit_key(edge))
        self.direction_counts[_spatial_direction_index(edge)] += 1

    def remove(self, edge: PhysicalEdge) -> None:
        self.edges.remove(edge)
        self.orbits.remove(_orbit_key(edge))
        self.direction_counts[_spatial_direction_index(edge)] -= 1


def _target_vertex_boundary(centers: Mapping[int, int]) -> set[Point]:
    """Return XOR of the sixteen hard affine center lines."""
    boundary: set[Point] = set()
    for x_coordinate in range(P):
        for y_coordinate in range(P):
            point = x_coordinate, y_coordinate
            incidence = sum(
                _functional_value(DIRECTIONS[index], point)
                == centers[index]
                for index in HARD
            )
            if incidence % 2:
                boundary.add(point)
    return boundary


def _boundary(edges: Iterable[PhysicalEdge]) -> set[Point]:
    boundary: set[Point] = set()
    for first, second in edges:
        boundary.symmetric_difference_update((first, second))
    return boundary


def _add_balanced_matching(
    builder: _LedgerGraphBuilder,
    vertices: set[Point],
    seed: int,
) -> None:
    random_source = random.Random(seed)
    remaining = set(vertices)
    while remaining:
        first = min(remaining)
        candidates: list[tuple[int, float, Point, PhysicalEdge]] = []
        for second in sorted(remaining):
            if second == first:
                continue
            edge = _physical_edge(first, second)
            if builder.allowed(edge):
                direction = _spatial_direction_index(edge)
                candidates.append(
                    (
                        builder.direction_counts[direction],
                        random_source.random(),
                        second,
                        edge,
                    )
                )
        if not candidates:
            raise ArithmeticError("the deterministic boundary matching failed")
        _count, _tie_break, second, edge = min(candidates)
        builder.insert(edge)
        remaining.remove(first)
        remaining.remove(second)


def _add_direction_triangle(
    builder: _LedgerGraphBuilder, direction: int
) -> None:
    kernel = _kernel_vector(direction)
    for x_coordinate in range(P):
        for y_coordinate in range(P):
            base = x_coordinate, y_coordinate
            for step in range(1, P):
                points = (
                    base,
                    _add(base, _scale(step, kernel)),
                    _add(base, _scale(2 * step, kernel)),
                )
                if len(set(points)) != 3:
                    continue
                edges = tuple(
                    _physical_edge(points[index], points[(index + 1) % 3])
                    for index in range(3)
                )
                if (
                    len({_orbit_key(edge) for edge in edges}) == 3
                    and all(builder.allowed(edge) for edge in edges)
                ):
                    for edge in edges:
                        builder.insert(edge)
                    return
    raise ArithmeticError("no clean direction triangle was found")


def _add_direction_four_cycle(
    builder: _LedgerGraphBuilder, direction: int
) -> None:
    kernel = _kernel_vector(direction)
    step_patterns = ((1, 2, 3), (1, 3, 7), (2, 5, 11))
    for x_coordinate in range(P):
        for y_coordinate in range(P):
            base = x_coordinate, y_coordinate
            for steps in step_patterns:
                points = (base,) + tuple(
                    _add(base, _scale(step, kernel)) for step in steps
                )
                if len(set(points)) != 4:
                    continue
                edges = tuple(
                    _physical_edge(points[index], points[(index + 1) % 4])
                    for index in range(4)
                )
                if (
                    len({_orbit_key(edge) for edge in edges}) == 4
                    and all(builder.allowed(edge) for edge in edges)
                ):
                    for edge in edges:
                        builder.insert(edge)
                    return
    raise ArithmeticError("no clean direction four-cycle was found")


def _subdivide_direction_edge(
    builder: _LedgerGraphBuilder,
    direction: int,
    protected_edge: PhysicalEdge,
) -> None:
    """Replace one edge by a three-edge path, preserving boundary and +2 count."""
    candidates = sorted(
        edge
        for edge in builder.edges
        if _spatial_direction_index(edge) == direction
        and edge != protected_edge
    )
    kernel = _kernel_vector(direction)
    for old_edge in candidates:
        builder.remove(old_edge)
        first, last = old_edge
        for first_step in range(1, P):
            middle_first = _add(first, _scale(first_step, kernel))
            if middle_first in (first, last):
                continue
            for second_step in range(1, P):
                middle_second = _add(first, _scale(second_step, kernel))
                if middle_second in (first, middle_first, last):
                    continue
                path = (
                    _physical_edge(first, middle_first),
                    _physical_edge(middle_first, middle_second),
                    _physical_edge(middle_second, last),
                )
                if (
                    len({_orbit_key(edge) for edge in path}) == 3
                    and all(builder.allowed(edge) for edge in path)
                ):
                    for edge in path:
                        builder.insert(edge)
                    return
        builder.insert(old_edge, allow_fixed=_is_fixed(old_edge))
    raise ArithmeticError("no clean three-edge subdivision was found")


def construct_coupled_ledger_graph() -> tuple[
    tuple[PhysicalEdge, ...], dict[int, int], tuple[int, ...]
]:
    """Construct the exact common graph used by the congruence replay."""
    centers = {index: 1 for index in HARD}
    target_boundary = _target_vertex_boundary(centers)
    builder = _LedgerGraphBuilder()

    kernel = _kernel_vector(FIXED_DIRECTION)
    fixed_edge = _physical_edge(kernel, _negative(kernel))
    builder.insert(fixed_edge, allow_fixed=True)

    matching_boundary = set(target_boundary)
    matching_boundary.symmetric_difference_update(fixed_edge)
    _add_balanced_matching(builder, matching_boundary, CONSTRUCTION_SEED)

    # A direction triangle has zero boundary and toggles only that direction's
    # edge-count parity.
    for direction, quota in enumerate(TARGET_PROFILE):
        if builder.direction_counts[direction] % 2 != quota % 2:
            _add_direction_triangle(builder, direction)
    counts_before_fill = tuple(builder.direction_counts)
    if any(
        count > quota
        for count, quota in zip(counts_before_fill, TARGET_PROFILE, strict=True)
    ):
        raise ArithmeticError("the deterministic parity graph exceeded a quota")

    # Within a direction fibre, a four-cycle seeds zero boundary if necessary;
    # replacing an edge by a length-three path then adds two edges without
    # changing boundary, direction, fixed count, or inversion-orbit cleanliness.
    for direction, quota in enumerate(TARGET_PROFILE):
        if builder.direction_counts[direction] == 0 and quota:
            _add_direction_four_cycle(builder, direction)
        while builder.direction_counts[direction] < quota:
            if quota - builder.direction_counts[direction] < 2:
                raise ArithmeticError("a direction fill lost its parity")
            _subdivide_direction_edge(builder, direction, fixed_edge)

    edges = tuple(sorted(builder.edges))
    if tuple(builder.direction_counts) != TARGET_PROFILE:
        raise ArithmeticError("the exact parallel profile changed")
    if _boundary(edges) != target_boundary:
        raise ArithmeticError("the graph boundary changed during filling")
    if sum(_is_fixed(edge) for edge in edges) != 1:
        raise ArithmeticError("the graph lost its unique fixed edge")
    if len({_orbit_key(edge) for edge in edges}) != len(edges):
        raise ArithmeticError("the graph acquired a doubled inversion orbit")
    return edges, centers, counts_before_fill


def _projected_boundary(
    edges: Iterable[PhysicalEdge], direction: int
) -> set[int]:
    functional = DIRECTIONS[direction]
    boundary: set[int] = set()
    for first, second in edges:
        first_label = _functional_value(functional, first)
        second_label = _functional_value(functional, second)
        if first_label != second_label:
            boundary.symmetric_difference_update((first_label, second_label))
    return boundary


def normalized_row_coefficients(
    edges: Iterable[PhysicalEdge],
    direction: int,
    centers: Mapping[int, int],
) -> dict[LabelEdge, int]:
    """Return the normalized graph row, plus the literal star when hard."""
    functional = DIRECTIONS[direction]
    epsilon = DIRECTION_SIGNS[direction]
    coefficients: Counter[LabelEdge] = Counter()
    for edge in edges:
        first = _functional_value(functional, edge[0])
        second = _functional_value(functional, edge[1])
        if first != second:
            coefficients[_label_edge(first, second)] += (
                epsilon * paley_edge_sign(P, edge)
            )
    if direction in HARD:
        center = centers[direction]
        for other in range(P):
            if other != center:
                coefficients[_label_edge(center, other)] += 1
    return {
        edge: coefficient
        for edge, coefficient in coefficients.items()
        if coefficient
    }


def _canonical_json_hash(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _row_lattice_replay(
    coefficients: Mapping[LabelEdge, int], direction: int
) -> tuple[dict[str, object], tuple[tuple[CompactAtom, int], ...]]:
    hard = direction in HARD
    quota = TARGET_PROFILE[direction]
    required_compact = quota - (3 if hard else 9)
    positive_triangles = (
        () if hard else tuple(islice(combinations(range(P), 3), 6))
    )
    compact_target: Counter[LabelEdge] = Counter(coefficients)
    for first, second, third in positive_triangles:
        for edge in (
            _label_edge(first, second),
            _label_edge(first, third),
            _label_edge(second, third),
        ):
            compact_target[edge] -= 1
            if not compact_target[edge]:
                del compact_target[edge]

    compact_atoms = signed_compact_lattice_decomposition(compact_target)
    compact_record = tuple(sorted(compact_atoms.items()))
    degrees = [0] * P
    for (first, second), coefficient in coefficients.items():
        degrees[first] += coefficient
        degrees[second] += coefficient
    edge_sum = sum(coefficients.values())
    expected_edge_sum = (
        -required_compact if hard else 18 - required_compact
    )
    if (
        edge_sum != expected_edge_sum
        or sum(compact_atoms.values()) != required_compact
        or len(positive_triangles) != (0 if hard else 6)
        or any(degree % 2 for degree in degrees)
    ):
        raise ArithmeticError("a row lattice coordinate changed")
    return (
        {
            "direction_index": direction,
            "direction_type": "hard" if hard else "opposite",
            "parallel_quota": quota,
            "edge_sum": edge_sum,
            "required_compact_count": required_compact,
            "positive_triangle_count": len(positive_triangles),
            "all_signed_degrees_even": True,
            "signed_compact_support": len(compact_atoms),
            "signed_compact_l1": sum(
                abs(value) for value in compact_atoms.values()
            ),
            "minimum_compact_coefficient": min(compact_atoms.values()),
            "maximum_compact_coefficient": max(compact_atoms.values()),
            "signed_compact_decomposition_sha256": _canonical_json_hash(
                compact_record
            ),
        },
        compact_record,
    )


def coupled_atom_congruence_certificate() -> dict[str, object]:
    """Replay the shared graph and every signed integer atom decomposition."""
    edges, centers, counts_before_fill = construct_coupled_ledger_graph()
    target_boundary = _target_vertex_boundary(centers)
    rows: list[dict[str, object]] = []
    atom_records: list[tuple[tuple[CompactAtom, int], ...]] = []
    for direction in range(P + 1):
        expected_boundary = (
            set(range(P)) - {centers[direction]}
            if direction in HARD
            else set()
        )
        if _projected_boundary(edges, direction) != expected_boundary:
            raise ArithmeticError("a projected graph boundary changed")
        row, atom_record = _row_lattice_replay(
            normalized_row_coefficients(edges, direction, centers),
            direction,
        )
        rows.append(row)
        atom_records.append(atom_record)

    fixed_edges = tuple(edge for edge in edges if _is_fixed(edge))
    graph_hash = _canonical_json_hash(edges)
    atom_hash = _canonical_json_hash(atom_records)
    minimum_coefficient = min(
        int(row["minimum_compact_coefficient"]) for row in rows
    )
    proved = bool(
        len(edges) == sum(TARGET_PROFILE) == 479
        and len(target_boundary) == 452
        and len(fixed_edges) == 1
        and _spatial_direction_index(fixed_edges[0]) == FIXED_DIRECTION
        and len({_orbit_key(edge) for edge in edges}) == len(edges)
        and all(row["all_signed_degrees_even"] for row in rows)
        and minimum_coefficient < 0
    )
    if not proved:
        raise ArithmeticError("the coupled congruence certificate changed")
    return {
        "p": P,
        "endpoint": "top t=177 hard-fixed j=0, f=1, d=0 ledger",
        "parallel_profile": TARGET_PROFILE,
        "graph_edge_count": len(edges),
        "fixed_direction_index": FIXED_DIRECTION,
        "fixed_edge": fixed_edges[0],
        "fixed_edge_count": len(fixed_edges),
        "doubled_nonfixed_orbit_count": 0,
        "hard_centers": tuple(sorted(centers.items())),
        "vertex_boundary_formula": "XOR_{N hard} 1[N(v)=j_N]",
        "vertex_boundary_weight": len(target_boundary),
        "direction_counts_before_even_fill": counts_before_fill,
        "graph_sha256": graph_hash,
        "row_count": len(rows),
        "rows": tuple(rows),
        "all_rows_in_required_signed_integer_atom_lattice": True,
        "signed_atom_decompositions_sha256": atom_hash,
        "total_signed_compact_support": sum(
            int(row["signed_compact_support"]) for row in rows
        ),
        "total_signed_compact_l1": sum(
            int(row["signed_compact_l1"]) for row in rows
        ),
        "minimum_compact_coefficient": minimum_coefficient,
        "maximum_compact_coefficient": max(
            int(row["maximum_compact_coefficient"]) for row in rows
        ),
        "compact_lattice_theorem": (
            "Z-span{K(a,b;c)} = {integer edge vectors with even signed degrees}"
        ),
        "mod2_audit": (
            "survives: triangle supports span the cycle space, and the common "
            "graph boundary gives exactly the hard-star/opposite row boundaries"
        ),
        "mod3_audit": (
            "survives: K(a,c;b)+K(b,c;a)=-2e_ab=e_ab mod 3; "
            "the atom-count coordinate follows from the exact edge sum"
        ),
        "pure_lattice_congruence_cut_found": False,
        "nonnegative_atom_decomposition_constructed": False,
        "nonnegative_chvatal_gomory_cut_excluded": False,
        "reason_nonnegative_problem_remains": (
            "the exact signed replay has negative compact multiplicities"
        ),
        "residual_ii_witness": False,
        "proved": proved,
    }


if __name__ == "__main__":
    print(
        json.dumps(
            coupled_atom_congruence_certificate(),
            indent=2,
            sort_keys=True,
        )
    )
