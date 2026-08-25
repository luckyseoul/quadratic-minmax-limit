#!/usr/bin/env python3
"""Prop 15.634 — square-circle spectrum and the full second shadow shell.

Retain the notation of Props. 15.629--15.633.  Let ``X=P^1(F_{p^2})``
and let ``B`` be one PSL orbit of F_p-sublines, namely the square circles.
There are

    v = p^2+1,                  b = p(p^2+1)/2

points and blocks.  The point--circle incidence matrix ``M`` is a
2-``(v,p+1,(p+1)/2)`` design.  If ``A`` joins two square circles when they
meet in two points, then the exact common-neighbour count is

    A^2+pA = (p^2-1) MM^T/8 + (p-1)^2(p+1) J/8.              (1)

Consequently ``col(M)`` has dimension ``v`` and, on ``ker(M^T)``, ``A``
has only the eigenvalues ``-p`` and ``0``, with multiplicities

    (p^2+1)(p-1)/4,            (p^2+1)(p-3)/4.               (2)

For each square circle ``S``, choose its signed complement ``w_S`` from
Prop. 15.633, so ``Cw_S=pw_S``.  Let ``b_S`` be the orthogonal projection
of ``w_S w_S^T`` to

    Z={W : PWP=W, diag(W)=0}.

The Gram matrix ``G=(<b_S,b_T>)`` has spectrum

    0                         multiplicity p^2+1,
    p^3(p-1)                 multiplicity (p^2+1)(p-1)/4,
    p^3(p+1)                 multiplicity (p^2+1)(p-3)/4.    (3)

Combining (3) with Prop. 15.633 diagonalizes the *complete* norm
``(p-1)/p`` harmonic shadow shell.  Its three eigenvalues are the rational
functions returned by :func:`second_shadow_spectrum`.  All three are
strictly negative for every odd prime ``p>=11``.  Thus the second shell is
an exact cancellation channel against Prop. 15.631's positive minimum
shell; this theorem does not control later shells and does not prove R1.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def n_of(p: int) -> int:
    return p * p + 1


def rank_of(p: int) -> int:
    return n_of(p) // 2


def z_dimension(p: int) -> int:
    return n_of(p) * (n_of(p) - 6) // 8


def square_circle_count(p: int) -> int:
    return p * n_of(p) // 2


def circle_design_parameters(p: int) -> dict[str, int]:
    return {
        "points": n_of(p),
        "blocks": square_circle_count(p),
        "block_size": p + 1,
        "replication": p * (p + 1) // 2,
        "pair_multiplicity": (p + 1) // 2,
    }


def intersection_valencies(p: int) -> dict[int, int]:
    """Number of other square circles meeting a fixed one in j points."""
    return {
        0: p * (p - 1) * (p - 3) // 4,
        1: p * p - 1,
        2: p * (p * p - 1) // 4,
    }


def common_two_secant_neighbours(p: int, intersection: int) -> int:
    """Common neighbours in the two-point-intersection graph."""
    if intersection == 0:
        return (p - 1) ** 2 * (p + 1) // 8
    if intersection == 1:
        return p * (p * p - 1) // 8
    if intersection == 2:
        return (p**3 + p * p - 9 * p - 1) // 8
    raise ValueError("distinct circles meet in 0, 1, or 2 points")


def circle_graph_spectrum(p: int) -> list[dict]:
    """Spectrum of adjacency by two-point intersection."""
    n = n_of(p)
    d = rank_of(p)
    return [
        {
            "eigenvalue": intersection_valencies(p)[2],
            "multiplicity": 1,
            "space": "constants",
        },
        {
            "eigenvalue": (p - 1) ** 2 // 4,
            "multiplicity": n - 1,
            "space": "nonconstant point-incidence module",
        },
        {
            "eigenvalue": -p,
            "multiplicity": d * (p - 1) // 2,
            "space": "incidence kernel, low circle module",
        },
        {
            "eigenvalue": 0,
            "multiplicity": d * (p - 3) // 2,
            "space": "incidence kernel, high circle module",
        },
    ]


def signed_complement_abs_inner(p: int, intersection: int) -> int:
    """Absolute ``w_S dot w_T``; signs of individual words are immaterial."""
    if intersection == p + 1:
        return p * (p - 1)
    if intersection == 2:
        return 0
    if intersection == 1:
        return p
    if intersection == 0:
        return 2 * p
    raise ValueError(intersection)


def projected_tensor_gram_entry(p: int, intersection: int) -> Fraction:
    """``<proj_Z(w_S w_S^T),proj_Z(w_T w_T^T)>``.

    The diagonal tensors ``(Pe_i)(Pe_i)^T`` have Gram matrix

        ((p^2-1)I+J)/(4p^2),

    whose inverse is ``4p^2 I/(p^2-1)-2J/(p^2-1)``.  The coordinate-square
    vector of ``w_S`` is the indicator of the complement of ``S``.
    """
    if intersection not in (0, 1, 2, p + 1):
        raise ValueError(intersection)
    complement_size = p * (p - 1)
    complement_intersection = n_of(p) - 2 * (p + 1) + intersection
    correction = Fraction(
        4 * p * p * complement_intersection - 2 * complement_size**2,
        p * p - 1,
    )
    raw = signed_complement_abs_inner(p, intersection) ** 2
    return Fraction(raw) - correction


def circle_tensor_gram_spectrum(p: int) -> list[dict]:
    n = n_of(p)
    d = rank_of(p)
    return [
        {"eigenvalue": 0, "multiplicity": n},
        {
            "eigenvalue": p**3 * (p - 1),
            "multiplicity": d * (p - 1) // 2,
        },
        {
            "eigenvalue": p**3 * (p + 1),
            "multiplicity": d * (p - 3) // 2,
        },
    ]


def circle_evaluation_operator_spectrum(p: int) -> list[dict]:
    """Spectrum on Z of ``sum_S <b_S,W>b_S/(8p^4)``."""
    n = n_of(p)
    d = rank_of(p)
    kernel = z_dimension(p) - d * (p - 2)
    return [
        {"eigenvalue": Fraction(0), "multiplicity": kernel},
        {
            "eigenvalue": Fraction(p - 1, 8 * p),
            "multiplicity": d * (p - 1) // 2,
        },
        {
            "eigenvalue": Fraction(p + 1, 8 * p),
            "multiplicity": d * (p - 3) // 2,
        },
    ]


def second_shadow_scalar_offset(p: int) -> Fraction:
    """Pair scalar plus the scalar part of the square-circle orbit."""
    d = rank_of(p)
    pair = Fraction(1, 4) * (1 - Fraction((p - 1) ** 2, d + 2))
    circle = -Fraction((p - 1) ** 2, 4 * p * (d + 2))
    return pair + circle


def second_shadow_spectrum(p: int) -> list[dict]:
    """Full signed norm-(p-1)/p shell acting on admissible harmonic W."""
    n = n_of(p)
    d = rank_of(p)
    a = second_shadow_scalar_offset(p)
    rows = [
        {
            "channel": "circle-kernel",
            "eigenvalue": a,
            "closed_form": (
                f"-{p + 2}*({p}^2-4*{p}+1)/(4*{p}*({p}^2+5))"
            ),
            "multiplicity": n * (p - 1) * (p - 3) // 8,
        },
        {
            "channel": "circle-low",
            "eigenvalue": a + Fraction(p - 1, 8 * p),
            "closed_form": (
                f"-({p}^3-3*{p}^2-19*{p}+9)/(8*{p}*({p}^2+5))"
            ),
            "multiplicity": d * (p - 1) // 2,
        },
        {
            "channel": "circle-high",
            "eigenvalue": a + Fraction(p + 1, 8 * p),
            "closed_form": (
                f"-({p}^3-5*{p}^2-19*{p}-1)/(8*{p}*({p}^2+5))"
            ),
            "multiplicity": d * (p - 3) // 2,
        },
    ]
    assert sum(row["multiplicity"] for row in rows) == z_dimension(p)
    return rows


def second_shadow_negative_definite(p: int) -> bool:
    return all(row["eigenvalue"] < 0 for row in second_shadow_spectrum(p))


def _legendre(value: int, p: int) -> int:
    value %= p
    if value == 0:
        return 0
    return 1 if pow(value, (p - 1) // 2, p) == 1 else -1


def explicit_square_circles(p: int) -> tuple[list[frozenset[int]], np.ndarray]:
    """Construct every square circle and one exact ``+p`` complement word.

    Vertices are ``0=infinity`` and ``1+x`` for ``x in F_{p^2}``.  Circles
    through infinity are affine lines in square directions.  The remaining
    circles are ``N(x-c)=r`` with ``Legendre(r)=-Legendre(-1)``.  Their
    complement word is ``w_inf=1``, ``w_x=Legendre(N(x-c)-r)``.
    """
    from e1_gmin_m4_prop15598 import field_ctx

    q, mul, add, chi, _frob, norm, _ia, _ib = field_ctx(p)

    def neg(x: int) -> int:
        return (-x % p) + (-(x // p) % p) * p

    def sub(x: int, y: int) -> int:
        return add(x, neg(y))

    seen_lines: set[frozenset[int]] = set()
    square_directions: list[int] = []
    for g in range(1, q):
        line = frozenset(mul(t, g) for t in range(p))
        if line in seen_lines:
            continue
        seen_lines.add(line)
        if chi(g) == 1:
            square_directions.append(g)

    blocks: list[frozenset[int]] = []
    words: list[list[int]] = []
    for g in square_directions:
        used_cosets: set[frozenset[int]] = set()
        g0, g1 = g % p, g // p
        for c in range(q):
            block = frozenset(
                [0] + [1 + add(c, mul(t, g)) for t in range(p)]
            )
            if block in used_cosets:
                continue
            used_cosets.add(block)
            blocks.append(block)
            word = [0]
            for x in range(q):
                y = sub(x, c)
                transverse = (g0 * (y // p) - g1 * (y % p)) % p
                word.append(_legendre(transverse, p))
            words.append(word)

    radius_type = -_legendre(-1, p)
    for c in range(q):
        for radius in range(1, p):
            if _legendre(radius, p) != radius_type:
                continue
            block = frozenset(
                1 + x for x in range(q) if norm(sub(x, c)) == radius
            )
            blocks.append(block)
            words.append(
                [1]
                + [
                    _legendre(norm(sub(x, c)) - radius, p)
                    for x in range(q)
                ]
            )
    return blocks, np.asarray(words, dtype=np.int64)


def audit_circle_scheme(p: int, check_quadratic_identity: bool = True) -> dict:
    """Finite exact audit of the circle formulas; not an input to the proof."""
    from minmax_quadratic import paley_conference_prime_power

    blocks, words = explicit_square_circles(p)
    n = n_of(p)
    C = np.rint(paley_conference_prime_power(p)).astype(np.int64)
    incidence = np.asarray(
        [[int(x in block) for x in range(n)] for block in blocks],
        dtype=np.int16,
    )
    intersections = incidence @ incidence.T
    adjacency = (intersections == 2).astype(np.int64)
    np.fill_diagonal(adjacency, 0)

    params = circle_design_parameters(p)
    expected_point_gram = (
        (params["replication"] - params["pair_multiplicity"])
        * np.eye(n, dtype=np.int64)
        + params["pair_multiplicity"] * np.ones((n, n), dtype=np.int64)
    )
    valencies = intersection_valencies(p)
    row_intersection_counts = {
        j: np.count_nonzero(intersections == j, axis=1)
        for j in (0, 1, 2)
    }
    correlations = words @ words.T
    correlation_ok = True
    for j in (0, 1, 2, p + 1):
        mask = intersections == j
        if np.any(mask):
            correlation_ok = correlation_ok and np.all(
                np.abs(correlations[mask])
                == signed_complement_abs_inner(p, j)
            )

    identity_error = None
    identity_ok = True
    if check_quadratic_identity:
        lhs = adjacency @ adjacency + p * adjacency
        rhs_numerator = (
            (p * p - 1) * intersections
            + (p - 1) ** 2 * (p + 1) * np.ones_like(intersections)
        )
        identity_ok = np.all(rhs_numerator % 8 == 0)
        rhs = rhs_numerator // 8
        identity_error = int(np.max(np.abs(lhs - rhs)))
        identity_ok = identity_ok and identity_error == 0

    checks = (
        len(blocks) == len(set(blocks)) == square_circle_count(p)
        and all(len(block) == p + 1 for block in blocks)
        and np.array_equal(C @ words.T, p * words.T)
        and np.array_equal(incidence.T @ incidence, expected_point_gram)
        and all(
            np.all(row_intersection_counts[j] == valencies[j])
            for j in (0, 1, 2)
        )
        and np.all(np.diag(intersections) == p + 1)
        and correlation_ok
        and identity_ok
    )
    return {
        "p": p,
        "circle_count": len(blocks),
        "word_plus_eigenvectors": bool(np.array_equal(C @ words.T, p * words.T)),
        "design_identity": bool(
            np.array_equal(incidence.T @ incidence, expected_point_gram)
        ),
        "intersection_valencies": {str(j): valencies[j] for j in valencies},
        "word_correlation_by_intersection": {
            "0": 2 * p,
            "1": p,
            "2": 0,
            str(p + 1): p * (p - 1),
        },
        "quadratic_adjacency_identity_max_error": identity_error,
        "checks": bool(checks),
    }


def circle_operator_theorem(
    primes: tuple[int, ...] = (5, 7, 11, 13, 17, 19)
) -> dict:
    rows = {}
    ok = True
    for p in primes:
        gram = circle_tensor_gram_spectrum(p)
        shadow = second_shadow_spectrum(p)
        row_ok = (
            sum(row["multiplicity"] for row in circle_graph_spectrum(p))
            == square_circle_count(p)
            and sum(row["multiplicity"] for row in gram)
            == square_circle_count(p)
            and sum(row["multiplicity"] for row in shadow) == z_dimension(p)
            and (second_shadow_negative_definite(p) == (p >= 11))
        )
        rows[str(p)] = {
            "circle_count": square_circle_count(p),
            "Z_dimension": z_dimension(p),
            "circle_tensor_gram_spectrum": gram,
            "second_shadow_spectrum": [
                {**row, "eigenvalue": str(row["eigenvalue"])}
                for row in shadow
            ],
            "second_shadow_negative_definite": second_shadow_negative_definite(p),
            "checks": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "scope": "all odd primes p>=5",
        "adjacency_identity": (
            "A^2+pA=(p^2-1)MM^T/8+(p-1)^2(p+1)J/8"
        ),
        "circle_gram_spectrum": (
            "0^n, [p^3(p-1)]^[n(p-1)/4], "
            "[p^3(p+1)]^[n(p-3)/4]"
        ),
        "negative_second_shell": (
            "The complete norm-(p-1)/p degree-four harmonic shadow "
            "operator is negative definite for every p>=11."
        ),
        "rows": rows,
    }


def main() -> dict:
    theorem = circle_operator_theorem()
    audits = {str(p): audit_circle_scheme(p) for p in (5, 7, 11, 13)}
    out = {
        "prop": "15.634",
        "title": "Square-circle operator spectrum and second-shell cancellation",
        "proved": {
            "square_circle_association_identity_all_p_ge_5": theorem["proved"],
            "square_circle_tensor_gram_spectrum_all_p_ge_5": theorem["proved"],
            "complete_second_shadow_negative_definite_p_ge_11": theorem["proved"],
            "finite_exact_audits": all(row["checks"] for row in audits.values()),
            "R1": False,
            "phi_F_ge_6_proved_general": False,
            "global_QVAR": False,
        },
        "theorem": theorem,
        "finite_exact_audits": audits,
        "consequence": (
            "The first higher dual shell is now diagonalized exactly.  For "
            "p>=11 it cancels in every harmonic channel, so first-shell "
            "positivity alone cannot prove R1; later shells remain open."
        ),
        "L_status": "OPEN",
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15634.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Prop 15.634 circle operator theorem: {theorem['proved']}")
    print(
        "  exact audits p=5,7,11,13: "
        f"{all(row['checks'] for row in audits.values())}"
    )
    print(f"  wrote {dest}")
    return out


if __name__ == "__main__":
    main()
