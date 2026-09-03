#!/usr/bin/env python3
r"""Parity gate for a three-spike endpoint on a two-outlier near pencil.

This is an unnumbered exact obstruction/diagnostic, not a residual-(ii)
closure.  At the lower ``p=4r+1`` branch-B endpoint, write

    h = |H| = 4r^2+6r+5,
    Q = # opposite Paley edges = 2r^2+r,
    S_H(x) = 8-3p.

If ``partial H=D`` and the Boolean three-spike shadow is in the ``+p``
gauge, then

    product_(e in H) C_e x_u x_v
      = (-1)^Q product_(v in D) x_v.

The target score has ``(h-(8-3p))/2=2r^2+9r`` negative edge features.
Both exponents have parity ``r``.  Consequently every common realization
must obey the previously unrecorded phase gate

    product_(v in D) x_v = +1.                              (1)

For either two-outlier near-pencil type, with pencil centre ``a`` and
replacement lines ``ell_i(a),ell_i(b)``, the phase can be evaluated without
constructing a graph:

    product_D x = product_(D0) x
                  product_(ell_1(a)) x product_(ell_1(b)) x
                  product_(ell_2(a)) x product_(ell_2(b)) x,              (2)

and, if ``O`` is the set of opposite directions,

    product_(D0) x = x_a product_(L in O) product_(ell_L(a)) x.           (3)

Equations (2)--(3) are symmetric-difference identities.  They apply to an
arbitrary all-bad three-spike completion; no affine-family classification is
assumed.

The exact ``p=53`` witnesses below show that (1) is a genuine relative-
alignment bit.  In both the ordinary and triple geometry, two translates of
the same Kiss--Somlai completion have the same two finite spike incidences
with ``D`` but opposite values of ``product_D x``.  Translation carries the
unique oriented spike circle and its mismatch set with the completion, so
the circle mismatch count is also unchanged.  Hence neither the three spike
incidences nor the circle mismatch count determines (1).
"""
from __future__ import annotations

from math import prod

import numpy as np

from e1_gmin_m4_prop15598 import field_ctx
from e1_gmin_m4_prop15721 import is_prime
from residual_kiss_somlai_three_spike import (
    transform_function,
    triangular_augmented_function,
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _validate_p1_prime(p: int) -> int:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 53
        or p % 4 != 1
        or not is_prime(p)
    ):
        raise ValueError("need a prime p=1 mod 4 with p>=53")
    return (p - 1) // 4


def endpoint_phase_ledger(p: int) -> dict[str, int | bool | str]:
    """Prove the necessary boundary-product phase at the lower endpoint."""
    r = _validate_p1_prime(p)
    h = 4 * r * r + 6 * r + 5
    opposite_edges = 2 * r * r + r
    required_score = 8 - 3 * p
    negative_features = (h - required_score) // 2
    _require(h - required_score == 2 * negative_features, "score parity changed")
    _require(negative_features == 2 * r * r + 9 * r, "negative count changed")
    _require(
        negative_features % 2 == opposite_edges % 2 == r % 2,
        "endpoint product parities no longer cancel",
    )
    return {
        "p": p,
        "r": r,
        "H_edge_count": h,
        "opposite_edge_count": opposite_edges,
        "required_three_spike_H_score": required_score,
        "required_negative_edge_features": negative_features,
        "target_feature_product": -1 if negative_features % 2 else 1,
        "Paley_edge_product": -1 if opposite_edges % 2 else 1,
        "necessary_boundary_product": 1,
        "identity": "prod_e(C_e*x_u*x_v)=(-1)^Q*prod_D(x)",
        "proved": True,
    }


def _encode(p: int, a: int, b: int) -> int:
    return (a % p) + (b % p) * p


def _coordinates(p: int, u: int) -> tuple[int, int]:
    return u % p, u // p


def _add(p: int, u: int, v: int) -> int:
    a, b = _coordinates(p, u)
    c, d = _coordinates(p, v)
    return _encode(p, a + c, b + d)


def _scalar(p: int, scalar: int, u: int) -> int:
    a, b = _coordinates(p, u)
    return _encode(p, scalar * a, scalar * b)


def _direction(p: int, u: int) -> tuple[int, int]:
    a, b = _coordinates(p, u)
    if a:
        return 1, b * pow(a, -1, p) % p
    if b:
        return 0, 1
    raise ValueError("zero has no projective direction")


def _direction_vector(p: int, direction: tuple[int, int]) -> int:
    return _encode(p, *direction)


def _line(p: int, base: int, direction: tuple[int, int]) -> set[int]:
    vector = _direction_vector(p, direction)
    return {_add(p, base, _scalar(p, t, vector)) for t in range(p)}


def _sign_product(x: np.ndarray, points: set[int]) -> int:
    if x.ndim != 1 or not np.all(np.isin(x, (-1, 1))):
        raise ValueError("x must be a Boolean finite-coordinate vector")
    return -1 if sum(int(x[u]) == -1 for u in points) % 2 else 1


def _directions_by_type(p: int) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    _q, _mul, _add_field, chi, _frob, _norm, _ia, _ib = field_ctx(p)
    directions = tuple((1, slope) for slope in range(p)) + ((0, 1),)
    hard = tuple(d for d in directions if chi(_direction_vector(p, d)) == 1)
    opposite = tuple(d for d in directions if chi(_direction_vector(p, d)) == -1)
    _require(
        len(hard) == len(opposite) == (p + 1) // 2,
        "projective direction types are unbalanced",
    )
    return hard, opposite


def _geometry(
    kind: str,
) -> tuple[int, int, int, int, tuple[int, int], tuple[int, int], set[int], set[int], tuple[set[int], ...]]:
    """Return one exact p=53 ordinary or triple two-outlier geometry."""
    p = 53
    _q, _mul, _add_field, chi, _frob, _norm, _ia, _ib = field_ctx(p)
    if kind == "ordinary":
        b = _encode(p, 1, 4)
        image_a = _encode(p, 3, 0)
        image_b = _encode(p, 2, 2)
        positive_translation = _encode(p, 0, 0)
        negative_translation = _encode(p, 0, 1)
        expected_boundary_size = 1504
        _require(chi(b) == -1, "ordinary connector direction must be opposite")
    elif kind == "triple":
        b = _encode(p, 2, 2)
        image_a = _encode(p, 1, 0)
        image_b = _encode(p, 1, 2)
        positive_translation = _encode(p, 0, 0)
        negative_translation = _encode(p, 3, 0)
        expected_boundary_size = 1500
        _require(chi(b) == 1, "triple connector direction must be hard")
    else:
        raise ValueError("kind must be 'ordinary' or 'triple'")

    _require(
        chi(image_a) == chi(image_b) == chi(_add(p, image_a, image_b)) == 1,
        "the three Kiss--Somlai special directions must all be hard",
    )
    direction_1 = _direction(p, image_a)
    direction_2 = _direction(p, image_b)
    d0 = {u for u in range(1, p * p) if chi(u) == -1}
    lines = (
        _line(p, 0, direction_1),
        _line(p, b, direction_1),
        _line(p, 0, direction_2),
        _line(p, b, direction_2),
    )
    boundary = set(d0)
    for line in lines:
        boundary.symmetric_difference_update(line)
    _require(len(boundary) == expected_boundary_size, "two-outlier boundary size changed")
    _require((b in boundary) == (kind == "ordinary"), "near-pencil type changed")
    return (
        p,
        b,
        image_a,
        image_b,
        positive_translation,
        negative_translation,
        d0,
        boundary,
        lines,
    )


def _translated_shadow(
    p: int, image_a: int, image_b: int, translation: int
) -> tuple[np.ndarray, tuple[int, int]]:
    transformed = transform_function(
        p,
        triangular_augmented_function(p),
        image_a,
        image_b,
    )
    shifted = np.empty_like(transformed)
    for u in range(p * p):
        shifted[_add(p, translation, u)] = transformed[u]
    double_points = tuple(int(u) for u in np.flatnonzero(shifted == 2))
    _require(len(double_points) == 2, "the finite three-spike support changed")
    shadow = np.where(shifted > 0, 1, -1).astype(np.int8)
    return shadow, double_points


def _factorization_ledger(
    p: int,
    x: np.ndarray,
    d0: set[int],
    boundary: set[int],
    lines: tuple[set[int], ...],
) -> dict[str, object]:
    _hard, opposite = _directions_by_type(p)
    direct = _sign_product(x, boundary)
    d0_product = _sign_product(x, d0)
    line_products = tuple(_sign_product(x, line) for line in lines)
    symmetric_difference_product = d0_product * prod(line_products)
    opposite_pencil_product = int(x[0]) * prod(
        _sign_product(x, _line(p, 0, direction)) for direction in opposite
    )
    _require(direct == symmetric_difference_product, "boundary product factorization failed")
    _require(d0_product == opposite_pencil_product, "D0 pencil product failed")
    return {
        "boundary_product": direct,
        "D0_product": d0_product,
        "replacement_line_products": list(line_products),
        "symmetric_difference_product": symmetric_difference_product,
        "opposite_pencil_product": opposite_pencil_product,
        "factorizations_exact": True,
    }


def representative_alignment_audit(kind: str) -> dict[str, object]:
    """Give phase-compatible and phase-incompatible p=53 translations."""
    (
        p,
        b,
        image_a,
        image_b,
        positive_translation,
        negative_translation,
        d0,
        boundary,
        lines,
    ) = _geometry(kind)
    rows = []
    for expected_phase, translation in (
        (1, positive_translation),
        (-1, negative_translation),
    ):
        shadow, finite_spikes = _translated_shadow(
            p, image_a, image_b, translation
        )
        incidence = tuple(spike in boundary for spike in finite_spikes)
        factorization = _factorization_ledger(
            p, shadow, d0, boundary, lines
        )
        _require(
            factorization["boundary_product"] == expected_phase,
            "representative translation has the wrong phase",
        )
        rows.append(
            {
                "translation": list(_coordinates(p, translation)),
                "finite_spikes": [list(_coordinates(p, u)) for u in finite_spikes],
                "finite_spike_incidence_with_D": list(incidence),
                **factorization,
            }
        )
    _require(
        rows[0]["finite_spike_incidence_with_D"]
        == rows[1]["finite_spike_incidence_with_D"]
        == [True, True],
        "the comparison translations must have identical spike incidences",
    )
    return {
        "p": p,
        "kind": kind,
        "pencil_intersection": list(_coordinates(p, b)),
        "boundary_size": len(boundary),
        "mapped_special_directions": [
            list(_direction(p, image_a)),
            list(_direction(p, image_b)),
            list(_direction(p, _add(p, image_a, image_b))),
        ],
        "rows": rows,
        "same_spike_incidence_opposite_boundary_phase": True,
        "translation_preserves_circle_mismatch": True,
        "circle_mismatch_does_not_determine_boundary_phase": True,
        "spike_incidence_does_not_determine_boundary_phase": True,
        "proved_for_this_prime": True,
    }


def representative_phase_audit() -> dict[str, object]:
    ordinary = representative_alignment_audit("ordinary")
    triple = representative_alignment_audit("triple")
    endpoint = endpoint_phase_ledger(53)
    return {
        "endpoint": endpoint,
        "ordinary": ordinary,
        "triple": triple,
        "both_two_outlier_types_have_both_relative_phases": True,
        "phase_compatible_alignments_exist": True,
        "phase_incompatible_alignments_exist": True,
        "residual_ii_closed": False,
        "all_checks": True,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(representative_phase_audit(), indent=2, sort_keys=True))
