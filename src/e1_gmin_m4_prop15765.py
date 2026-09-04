#!/usr/bin/env python3
r"""Proposition 15.765 -- a nonaffine first-defect-shell point at ``p=11``.

Kiss--Somlai give a 33-point subset ``E0`` of ``F_11^2`` having exactly
four special (non-equidistributed) directions: the two axes and slopes
``+1,-1``.  Their displayed matrix has an empty horizontal row.  Add that
row, then apply ``T(x,y)=(x,2y)``.  The resulting 44-point set ``E`` has
exactly the four special directions

    infinity, 0, 2, -2.

In ``F_121=F_11[a]/(a^2-2)``, all four are square spatial directions.  Thus
``1_E`` is constant plus a ``+11`` Paley eigenfunction.  Direct exact
verification gives

    Q 1_E = 11 1_E - 4 1,

where ``Q_uv=chi(v-u)``.  For ``D=F_11^2\E`` this is

    Q 1_D = 11 1_D - 7 1,
    3 + 2 sum_(v in D) chi(v-u) = 11(2 1_D(u)-1).          (1)

Consequently, for the normalized Paley conference matrix

    C = [[0, 1^T], [1, Q]],

the integral vector ``y_infinity=3``, ``y_u=2 1_D(u)-1`` satisfies
``Cy=11y``.  Replacing its unique 3 by 1 gives a Boolean point ``x`` with

    Phi - q_C(x) = 22 = 2p.

This point is genuinely nonaffine: ``E`` has four special directions,
whereas a nontrivial union of parallel affine lines has exactly one.  In
particular ``D`` is not a union of seven parallel lines.  This disproves the
suggested classification of the full first-defect shell by affine aliases.

It is not a residual-(ii) counterexample.  Proposition 15.765 constructs one
Paley shell point, not a common switching set ``H`` satisfying the
all-deletions hypotheses, and it changes no global predicate.
"""
from __future__ import annotations

import hashlib
import json
import struct
from pathlib import Path
from typing import Iterable


P = 11
NONRESIDUE = 2
ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15765.json"

# Figure 4 (right-hand 0/1 matrix) in Kiss--Somlai, Section 6.  Rows are
# indexed by y=0,...,10 and columns by x=0,...,10.
KISS_SOMLAI_ROWS: tuple[tuple[int, ...], ...] = (
    (0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    (0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1),
    (1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0),
    (0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    (0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0),
)

Point = tuple[int, int]
Direction = int | None  # slope, with None denoting the vertical direction


def points() -> tuple[Point, ...]:
    """Return finite vertices in the canonical lexicographic ``(x,y)`` order."""
    return tuple((x, y) for x in range(P) for y in range(P))


def legendre(a: int) -> int:
    """Quadratic character of ``F_11``, extended by zero."""
    a %= P
    if a == 0:
        return 0
    return 1 if pow(a, (P - 1) // 2, P) == 1 else -1


def paley_character(z: Point) -> int:
    r"""Quadratic character of ``x+y*a`` in ``F_121``, ``a^2=2``.

    An element of a quadratic extension is a square iff its norm is a
    square in the base field, and ``N(x+y*a)=x^2-2y^2``.
    """
    x, y = z
    if x % P == 0 and y % P == 0:
        return 0
    return legendre(x * x - NONRESIDUE * y * y)


def kiss_somlai_base() -> frozenset[Point]:
    """The 33-point set read directly from the cited matrix."""
    if len(KISS_SOMLAI_ROWS) != P or any(len(row) != P for row in KISS_SOMLAI_ROWS):
        raise ArithmeticError("the cited matrix is not 11 by 11")
    if any(bit not in (0, 1) for row in KISS_SOMLAI_ROWS for bit in row):
        raise ArithmeticError("the cited matrix is not Boolean")
    return frozenset(
        (x, y)
        for y, row in enumerate(KISS_SOMLAI_ROWS)
        for x, bit in enumerate(row)
        if bit
    )


def exceptional_complement() -> frozenset[Point]:
    r"""Return the 44-point set ``E=T(E0 union {y=1})``."""
    base = kiss_somlai_base()
    added_line = frozenset((x, 1) for x in range(P))
    if base & added_line:
        raise ArithmeticError("the cited empty row is not empty")
    preimage = base | added_line
    return frozenset((x, 2 * y % P) for x, y in preimage)


def exceptional_set() -> frozenset[Point]:
    r"""Return ``D=F_11^2\E``, of cardinality 77."""
    return frozenset(points()) - exceptional_complement()


def line_profile(subset: frozenset[Point], direction: Direction) -> tuple[int, ...]:
    r"""Intersection sizes with the eleven lines in one direction.

    For finite slope ``m`` the lines are ``y=m*x+b``.  For ``None`` they
    are ``x=b``.
    """
    if direction is not None and not 0 <= direction < P:
        raise ValueError("a finite direction must be a residue modulo 11")
    out = []
    for b in range(P):
        if direction is None:
            out.append(sum((b, y) in subset for y in range(P)))
        else:
            out.append(
                sum((x, (direction * x + b) % P) in subset for x in range(P))
            )
    return tuple(out)


def special_direction_profiles(
    subset: frozenset[Point],
) -> dict[str, tuple[int, ...]]:
    """Return exactly the nonconstant parallel-class profiles."""
    answer: dict[str, tuple[int, ...]] = {}
    for direction in (*range(P), None):
        profile = line_profile(subset, direction)
        if len(set(profile)) > 1:
            answer["infinity" if direction is None else str(direction)] = profile
    return answer


def paley_convolution(subset: frozenset[Point], u: Point) -> int:
    """Return ``sum_(v in subset) chi(v-u)`` exactly."""
    ux, uy = u
    return sum(
        paley_character(((vx - ux) % P, (vy - uy) % P)) for vx, vy in subset
    )


def indicator_convolution_identity() -> bool:
    """Check both exact forms of (1) at all 121 finite vertices."""
    e_set = exceptional_complement()
    d_set = exceptional_set()
    for u in points():
        e = int(u in e_set)
        d = int(u in d_set)
        if paley_convolution(e_set, u) != P * e - 4:
            return False
        conv_d = paley_convolution(d_set, u)
        if conv_d != P * d - 7:
            return False
        if 3 + 2 * conv_d != P * (2 * d - 1):
            return False
    return True


def integral_eigenvector() -> tuple[int, ...]:
    """Return ``(3, 2*1_D-1)`` with infinity first."""
    d_set = exceptional_set()
    return (3, *(2 * int(u in d_set) - 1 for u in points()))


def boolean_shadow() -> tuple[int, ...]:
    """Replace the unique value 3 of the integral eigenvector by 1."""
    y = integral_eigenvector()
    return (1, *y[1:])


def conference_action(vector: tuple[int, ...]) -> tuple[int, ...]:
    """Apply the normalized order-122 Paley conference matrix exactly."""
    if len(vector) != P * P + 1:
        raise ValueError("the conference vector must have length 122")
    finite = vector[1:]
    out = [sum(finite)]
    finite_points = points()
    for i, u in enumerate(finite_points):
        value = vector[0]
        ux, uy = u
        for j, (vx, vy) in enumerate(finite_points):
            value += paley_character(((vx - ux) % P, (vy - uy) % P)) * finite[j]
        out.append(value)
    return tuple(out)


def conference_quadratic(vector: tuple[int, ...]) -> int:
    """Return ``vector^T C vector / 2`` as an exact integer."""
    acted = conference_action(vector)
    numerator = sum(a * b for a, b in zip(vector, acted, strict=True))
    if numerator % 2:
        raise ArithmeticError("conference quadratic numerator is odd")
    return numerator // 2


def paley_neighbor_counts(subset: frozenset[Point]) -> dict[str, tuple[int, ...]]:
    """Return sorted internal counts inside and outside ``subset``."""
    inside = []
    outside = []
    for u in points():
        ux, uy = u
        degree = sum(
            paley_character(((vx - ux) % P, (vy - uy) % P)) == 1
            for vx, vy in subset
        )
        (inside if u in subset else outside).append(degree)
    return {"inside": tuple(sorted(set(inside))), "outside": tuple(sorted(set(outside)))}


def coordinate_hash(subset: Iterable[Point]) -> str:
    """SHA-256 of sorted coordinates encoded as consecutive unsigned bytes."""
    payload = bytes(c for point in sorted(subset) for c in point)
    return hashlib.sha256(payload).hexdigest()


def vector_hash(vector: Iterable[int]) -> str:
    """SHA-256 of signed integers encoded as little-endian int16 values."""
    values = tuple(vector)
    payload = struct.pack(f"<{len(values)}h", *values)
    return hashlib.sha256(payload).hexdigest()


def theorem_record() -> dict[str, object]:
    """Return the complete exact Proposition 15.765 evidence record."""
    base = kiss_somlai_base()
    e_set = exceptional_complement()
    d_set = exceptional_set()
    profiles = special_direction_profiles(e_set)
    square_directions = {
        "infinity": paley_character((0, 1)),
        "0": paley_character((1, 0)),
        "2": paley_character((1, 2)),
        "9": paley_character((1, 9)),
    }
    y = integral_eigenvector()
    x = boolean_shadow()
    cy = conference_action(y)
    phi = P * (P * P + 1) // 2
    qx = conference_quadratic(x)
    neighbor_counts = paley_neighbor_counts(d_set)
    record = {
        "prop": "15.765",
        "title": "Nonaffine first-defect-shell Paley point at p=11",
        "status": "PROVED EXACT COUNTEREXAMPLE TO AFFINE-SHELL CLASSIFICATION",
        "primary_source": {
            "authors": "Gergely Kiss and Gabor Somlai",
            "title": "Special directions on the finite affine plane",
            "journal": "Designs, Codes and Cryptography 92 (2024), 2587-2597",
            "doi": "10.1007/s10623-024-01404-y",
            "arxiv": "2109.13992v3",
            "location": "Section 6, Figure 4, right-hand 11 by 11 matrix",
            "arxiv_pdf_sha256": "e35a184f64f4e7af03744a10a8cd3eed3ecf05308f4d2dda286136810e85df0b",
            "arxiv_source_sha256": "ebbf527ea4e16e42416de3239e8be5fb7635a1ca2458117bcfe2fd7f3a77e525",
        },
        "construction": {
            "field": "F_121=F_11[a]/(a^2-2)",
            "character": "chi(x+y*a)=Legendre(x^2-2*y^2)",
            "base_size": len(base),
            "added_disjoint_line": "{(x,1): x in F_11}",
            "linear_map": "T(x,y)=(x,2y)",
            "E_size": len(e_set),
            "D_size": len(d_set),
            "E_coordinate_encoding": "sorted (x,y), consecutive uint8",
            "E_coordinate_sha256": coordinate_hash(e_set),
            "D_coordinate_sha256": coordinate_hash(d_set),
        },
        "nonaffineness": {
            "special_direction_profiles": {key: list(value) for key, value in profiles.items()},
            "special_directions": list(profiles),
            "all_special_directions_are_square": all(v == 1 for v in square_directions.values()),
            "direction_characters": square_directions,
            "E_is_union_of_four_parallel_lines": False,
            "D_is_union_of_seven_parallel_lines": False,
            "certificate": "four special directions; a nontrivial parallel-line union has exactly one",
        },
        "exact_identities": {
            "Q_1E": "11*1_E-4*1",
            "Q_1D": "11*1_D-7*1",
            "pointwise_convolution_checked": indicator_convolution_identity(),
            "paley_quotient_D_E": [[40, 20], [35, 25]],
            "D_neighbor_counts": {
                key: list(value) for key, value in neighbor_counts.items()
            },
            "integral_vector_norm_squared": sum(value * value for value in y),
            "Cy_equals_11y": cy == tuple(P * value for value in y),
            "magnitude_three_coordinates": [i for i, value in enumerate(y) if abs(value) == 3],
            "integral_vector_int16le_sha256": vector_hash(y),
            "boolean_shadow_int16le_sha256": vector_hash(x),
            "Phi": phi,
            "q_C_boolean_shadow": qx,
            "first_defect": phi - qx,
            "first_defect_equals_2p": phi - qx == 2 * P,
        },
        "literature_boundary": {
            "bailey_cameron_gavrilyuk_goryainov_2019": (
                "Section 8 poses classification of (n-m)-perfect sets for MOLS graphs; "
                "its Latin-square classification is only m=3"
            ),
            "high_dimensional_boolean_degree_one_theorem_inapplicable": (
                "the OA(6,11) point-facet complex is not proper: its top up-map "
                "has 726-dimensional domain and 121-dimensional codomain"
            ),
        },
        "proved": {
            "nonaffine_first_defect_shell_point_exists": True,
            "all_first_shell_points_are_affine": False,
            "common_all_deletions_H_constructed": False,
            "residual_ii_closed": False,
            "e1_closed_general": False,
            "L": False,
        },
        "L_status": "OPEN",
    }
    if set(profiles) != {"0", "2", "9", "infinity"}:
        raise ArithmeticError("unexpected set of special directions")
    if record["exact_identities"]["D_neighbor_counts"] != {
        "inside": [40],
        "outside": [35],
    }:
        raise ArithmeticError("the Paley equitable-partition parameters failed")
    return record


def main() -> dict[str, object]:
    out = theorem_record()
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print("Prop. 15.765 nonaffine p=11 first shell: proved")
    print("  |D|=77, Cy=11y, Phi-q_C(x)=22")
    print("  affine-shell classification: false")
    print("  residual (ii): OPEN")
    print("wrote", EV)
    return out


if __name__ == "__main__":
    main()
