"""Exact constructor and regularizability witness for the sporadic Peisert graph.

This is the rank-3 graph ``G(23^2, 2)``.  Vertices are ``F_23^2`` and the
connection set is the orbit of ``(1, 0)`` under a concrete copy of
``Z(GL(2,23)) o SL(2,3)``.  No external graph or finite-field package is
needed; all checks use integer arithmetic modulo 23.

The module is a deduplication artifact, not a residual-(ii) closure: the
scalar subgroup makes the connection set a union of projective directions,
so this graph is in the already-regularizable linear-OA/PN class.
"""

from __future__ import annotations

import hashlib
import json
import struct
from collections import deque
from typing import Iterable

P = 23
V = tuple((a, b) for a in range(P) for b in range(P))

Matrix2 = tuple[int, int, int, int]
Vector2 = tuple[int, int]

# A quaternion generator, a binary-tetrahedral order-three generator, and a
# scalar of order 11.  Together they generate a group of order 264.
GENERATORS: tuple[Matrix2, ...] = (
    (0, -1, 1, 0),
    (8, 16, 17, 14),
    (2, 0, 0, 2),
)

PUBLISHED_CSV_URL = (
    "https://www.math.mun.ca/distanceregular/graphdata/"
    "sporadicpeisert529.am.csv"
)
PUBLISHED_CSV_SHA256 = (
    "fcdd847709adc1527374781fb81857bf3ce741c538e35f59126b1c6256e3cda6"
)
CONSTRUCTED_CSV_SHA256 = (
    "eeba88c7441385065a235fecc62240de25474f22ad15cab3dcc307b2cbf0c3bc"
)
CANONICAL_CERTIFICATE_SHA256 = (
    "64c08ce8c6cacedbf5201441a5230637b1961b19f43d784985b5a1d57332a3b0"
)


def matmul(x: Matrix2, y: Matrix2) -> Matrix2:
    """Multiply two 2 by 2 matrices modulo 23."""

    a, b, c, d = x
    e, f, g, h = y
    return (
        (a * e + b * g) % P,
        (a * f + b * h) % P,
        (c * e + d * g) % P,
        (c * f + d * h) % P,
    )


def act(m: Matrix2, v: Vector2) -> Vector2:
    a, b, c, d = m
    x, y = v
    return ((a * x + b * y) % P, (c * x + d * y) % P)


def group_closure(generators: Iterable[Matrix2] = GENERATORS) -> frozenset[Matrix2]:
    identity = (1, 0, 0, 1)
    seen = {identity}
    todo = deque([identity])
    generators = tuple(tuple(x % P for x in g) for g in generators)
    while todo:
        x = todo.popleft()
        for g in generators:
            y = matmul(x, g)
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return frozenset(seen)


def connection_set(group: Iterable[Matrix2] | None = None) -> frozenset[Vector2]:
    if group is None:
        group = group_closure()
    return frozenset(act(g, (1, 0)) for g in group)


def direction_points(direction: Vector2) -> frozenset[Vector2]:
    a, b = direction
    return frozenset(((t * a) % P, (t * b) % P) for t in range(1, P))


def projective_directions() -> tuple[Vector2, ...]:
    return tuple((1, m) for m in range(P)) + ((0, 1),)


def is_adjacent(v: Vector2, w: Vector2, dset: frozenset[Vector2]) -> bool:
    return ((w[0] - v[0]) % P, (w[1] - v[1]) % P) in dset


def common_neighbor_histogram(dset: frozenset[Vector2]) -> dict[int, int]:
    """Check the SRG identity by translation, without dense matrix products."""

    hist: dict[int, int] = {}
    for t in V[1:]:
        count = sum(((x - t[0]) % P, (y - t[1]) % P) in dset for x, y in dset)
        expected = 131 if t in dset else 132
        assert count == expected
        hist[count] = hist.get(count, 0) + 1
    return hist


def line_value(v: Vector2, direction: Vector2) -> int:
    """A linear functional whose fibres are parallel to ``direction``."""

    x, y = v
    a, b = direction
    return (b * x - a * y) % P


def boolean_witness(lam: int) -> tuple[int, ...]:
    """Return x with x_infinity=1 and Cx=lam*x for lam in {+23,-23}."""

    if lam == P:
        direction, levels = (1, 1), frozenset(range(11))
    elif lam == -P:
        direction, levels = (1, 0), frozenset(range(12))
    else:
        raise ValueError("lambda must be +23 or -23")
    return (1,) + tuple(-1 if line_value(v, direction) in levels else 1 for v in V)


def verify_boolean_witness(
    x: tuple[int, ...], lam: int, dset: frozenset[Vector2]
) -> None:
    assert len(x) == P * P + 1
    assert x[0] == 1 and set(x) <= {-1, 1}
    finite = x[1:]
    assert sum(finite) == lam
    for i, v in enumerate(V):
        neighbor_sum = sum(
            finite[j] for j, w in enumerate(V) if is_adjacent(v, w, dset)
        )
        # C has first row/column +1 and finite block S=J-I-2A.
        image = 1 + sum(finite) - finite[i] - 2 * neighbor_sum
        assert image == lam * finite[i]


def int16_sha256(x: Iterable[int]) -> str:
    return hashlib.sha256(b"".join(struct.pack("<h", value) for value in x)).hexdigest()


def constructed_csv_sha256(dset: frozenset[Vector2]) -> str:
    digest = hashlib.sha256()
    for v in V:
        row = ",".join("1" if is_adjacent(v, w, dset) else "0" for w in V)
        digest.update(row.encode("ascii") + b"\n")
    return digest.hexdigest()


def exact_report(check_csv_hash: bool = True) -> dict[str, object]:
    group = group_closure()
    assert len(group) == 264
    dset = connection_set(group)
    assert len(dset) == 264 and (0, 0) not in dset
    assert all(((-x) % P, (-y) % P) in dset for x, y in dset)

    directions = projective_directions()
    selected = tuple(d for d in directions if direction_points(d) <= dset)
    assert len(selected) == 12
    assert all(
        direction_points(d) <= dset or direction_points(d).isdisjoint(dset)
        for d in directions
    )

    plus = boolean_witness(P)
    minus = boolean_witness(-P)
    verify_boolean_witness(plus, P, dset)
    verify_boolean_witness(minus, -P, dset)
    csv_hash = constructed_csv_sha256(dset) if check_csv_hash else CONSTRUCTED_CSV_SHA256
    assert csv_hash == CONSTRUCTED_CSV_SHA256

    return {
        "status": "regularizable_linear_OA_PN_not_residual_ii_closure",
        "p": P,
        "graph_parameters": [529, 264, 131, 132],
        "conference_order": 530,
        "conference_identity": "C^2=529I",
        "group_order": len(group),
        "connection_size": len(dset),
        "selected_projective_directions": [list(d) for d in selected],
        "common_neighbor_histogram": {
            str(k): v for k, v in sorted(common_neighbor_histogram(dset).items())
        },
        "witnesses": {
            "+23": {
                "minus_count": sum(value == -1 for value in plus[1:]),
                "parallel_direction": [1, 1],
                "line_levels": list(range(11)),
                "int16_sha256": int16_sha256(plus),
            },
            "-23": {
                "minus_count": sum(value == -1 for value in minus[1:]),
                "parallel_direction": [1, 0],
                "line_levels": list(range(12)),
                "int16_sha256": int16_sha256(minus),
            },
        },
        "published_csv": {
            "url": PUBLISHED_CSV_URL,
            "sha256": PUBLISHED_CSV_SHA256,
        },
        "constructed_csv_sha256": csv_hash,
        "canonical_certificate_sha256": CANONICAL_CERTIFICATE_SHA256,
    }


def main() -> None:
    print(json.dumps(exact_report(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
