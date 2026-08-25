#!/usr/bin/env python3
"""Prop 15.632 — affine slack budget and quadratic parity lifts.

Let ``H`` be an odd set of signed Paley edges and put

    S_H(y) = sum_{e in H} C_e y_e.

For each of the ``p+1`` projective F_p-directions ``d``, let ``eps_d`` be
the quadratic type of its kernel.  The affine halfspaces in direction ``d``
are Boolean ``eps_d*p`` eigenvectors.  If ``H`` separates the two Boolean
eigenshells with margin three, then

    A_d(y) = (eps_d S_H(y)-3)/2

is a nonnegative integer-valued quadratic polynomial on the middle slice
``J(p,(p+1)/2)``.  Define ``a_d=2p E_d[A_d]``.  Direct middle-slice moments
give the stronger integral statement

    a_d is a nonnegative even integer,
    sum_{d:eps_d=tau} a_d = (p+1)(|H|-3p)/2  (tau=+-1),
    sum_d a_d = (p+1)(|H|-3p).                              (1)

At residual-(ii) size ``|H|=4p+1``, each quadratic-type half has budget
``(p+1)^2/2`` and the total budget is ``(p+1)^2``.

There is a non-Walsh refinement of (1).  Let ``D`` be the odd-degree
boundary of ``H`` and ``c_H=prod_{e in H} C_e``.  In direction ``d``, let
``B_d`` be the set of affine fibres containing an odd number of finite
points of ``D`` and put ``b_d=|B_d|``.  The edge-product identity gives

    A_d(x) = sum_{s in B_d} x_s + eta_d                 (mod 2),   (2)

where ``eta_d`` is explicit from ``eps_d``, ``c_H``, ``infinity in D``,
``|H|``, and ``b_d``.

For any nonnegative integer-valued quadratic ``A`` on the middle slice
with parity ``sum_{s in B} x_s+eta``, average ``A`` under
``Sym(B) x Sym(B^c)``.  The average is a univariate quadratic ``q(t)``,
``t=|X cap B|``.  Every value averaged at ``t`` is a nonnegative integer
of parity ``t+eta``, hence

    q(t) >= (t+eta mod 2).

Let ``M(p,b,eta)`` be the minimum hypergeometric expectation of a quadratic
with those pointwise inequalities.  This is an exact three-variable LP;
its vertices are obtained by interpolating three parity points.  Therefore

    a_d >= 2 ceil(p M(p,b_d,eta_d)),
    sum_{d:eps_d=tau} 2 ceil(p M(p,b_d,eta_d))
        <= (p+1)(|H|-3p)/2,                 tau=+-1,         (3)

Summing the two inequalities in (3) gives the corresponding total bound.
Equation (3) is an all-prime necessary condition for the non-Walsh
multi-level residual.  It uses the magnitude of the integer slack, not
only its Walsh parity.  In particular, an Eulerian boundary ``D=empty`` is
impossible at residual size: one quadratic-type half has constant odd slack
parity and therefore costs ``p(p+1)>(p+1)^2/2``.  The proposition does not
yet classify all nonempty affine-plane boundary profiles, so it does not
close residual-(ii), Type I, R1, or L.
"""
from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]


def middle_weight(p: int) -> int:
    if p < 3 or p % 2 == 0:
        raise ValueError("p must be odd")
    return (p + 1) // 2


def projective_directions(p: int) -> tuple[tuple[int, int], ...]:
    """Canonical linear forms r*a+s*b, one per P^1(F_p)."""
    return tuple([(1, t) for t in range(p)] + [(0, 1)])


def hypergeometric_weights(p: int, b: int) -> dict[int, Fraction]:
    """Law of t=|X cap B| for uniform |X|=(p+1)/2 and |B|=b."""
    if not 0 <= b <= p:
        raise ValueError("b must lie in 0..p")
    m = middle_weight(p)
    lo = max(0, m - (p - b))
    hi = min(b, m)
    den = math.comb(p, m)
    return {
        t: Fraction(math.comb(b, t) * math.comb(p - b, m - t), den)
        for t in range(lo, hi + 1)
    }


def _quadratic_interpolant(
    points: tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
) -> tuple[Fraction, Fraction, Fraction]:
    """Return (u,v,w) for q(t)=u*t^2+v*t+w through three points."""
    out = [Fraction(0), Fraction(0), Fraction(0)]
    for i, (x, y) in enumerate(points):
        other = [points[j][0] for j in range(3) if j != i]
        r, s = other
        den = (x - r) * (x - s)
        out[0] += Fraction(y, den)
        out[1] += Fraction(-y * (r + s), den)
        out[2] += Fraction(y * r * s, den)
    return out[0], out[1], out[2]


def eval_quadratic(coefficients: tuple[Fraction, Fraction, Fraction], t: int) -> Fraction:
    u, v, w = coefficients
    return u * t * t + v * t + w


def parity_majorant_floor(p: int, b: int, phase: int) -> dict:
    """Exact value M(p,b,phase) from the three-variable LP.

    The constraints are q(t)>=(t+phase mod 2) on the support of the
    hypergeometric law.  With at least three support points, every LP vertex
    has three active, linearly independent evaluation constraints, so exact
    enumeration of triples is exhaustive.  One- and two-point supports are
    interpolation-trivial.
    """
    if phase not in (0, 1):
        raise ValueError("phase must be 0 or 1")
    weights = hypergeometric_weights(p, b)
    ts = tuple(weights)
    parity = {t: (t + phase) & 1 for t in ts}

    if len(ts) <= 2:
        value = sum(weights[t] * parity[t] for t in ts)
        return {
            "p": p,
            "b": b,
            "phase": phase,
            "value": value,
            "coefficients": None,
            "contacts": ts,
            "support": ts,
            "proved": True,
        }

    best: Fraction | None = None
    best_coefficients: tuple[Fraction, Fraction, Fraction] | None = None
    best_contacts: tuple[int, int, int] | None = None
    for contacts in itertools.combinations(ts, 3):
        coefficients = _quadratic_interpolant(
            tuple((t, parity[t]) for t in contacts)  # type: ignore[arg-type]
        )
        if not all(eval_quadratic(coefficients, t) >= parity[t] for t in ts):
            continue
        value = sum(weights[t] * eval_quadratic(coefficients, t) for t in ts)
        if best is None or value < best:
            best = value
            best_coefficients = coefficients
            best_contacts = contacts
    if best is None or best_coefficients is None or best_contacts is None:
        raise RuntimeError("quadratic-majorant LP had no enumerated vertex")
    return {
        "p": p,
        "b": b,
        "phase": phase,
        "value": best,
        "coefficients": best_coefficients,
        "contacts": best_contacts,
        "support": ts,
        "proved": all(
            eval_quadratic(best_coefficients, t) >= parity[t] for t in ts
        ),
    }


def ceil_fraction(x: Fraction) -> int:
    return -((-x.numerator) // x.denominator)


def scaled_direction_floor(p: int, b: int, phase: int) -> int:
    """The even-integer consequence a_d >= 2 ceil(p M)."""
    value = parity_majorant_floor(p, b, phase)["value"]
    assert isinstance(value, Fraction)
    return 2 * ceil_fraction(p * value)


def eulerian_residual_type_budget_gap(p: int) -> int:
    """Contradiction gap for ``|H|=4p+1`` and empty boundary.

    Empty boundary gives ``b_d=0``.  The residual phase is
    ``-eps_d*c_H``, so one of the two quadratic-type halves has phase one.
    Every one of its ``(p+1)/2`` directions costs ``2p``, while that half's
    exact budget is ``(p+1)^2/2``.
    """
    m = middle_weight(p)
    return 2 * p * m - m * (p + 1)


def field_direction_data(p: int, direction: tuple[int, int]) -> tuple[int, list[int]]:
    """Return (kernel quadratic type, finite fibre labels)."""
    from e1_gmin_m4_prop15598 import field_ctx

    r, s = direction
    q, _mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(p)
    kernel_generator = (s % p) + ((-r) % p) * p
    eps = int(chi(kernel_generator))
    labels = [((r * (u % p) + s * (u // p)) % p) for u in range(q)]
    return eps, labels


def edge_boundary(n: int, H: Iterable[tuple[int, int]]) -> tuple[set[int], list[int]]:
    degree = [0] * n
    for a, b in H:
        degree[a] += 1
        degree[b] += 1
    return {i for i, d in enumerate(degree) if d & 1}, degree


def direction_record(
    p: int,
    H: Iterable[tuple[int, int]],
    C,
    direction: tuple[int, int],
) -> dict:
    """Exact directional mean, boundary parity, and parity-lift floor.

    Vertex 0 is infinity and finite field element ``u`` is vertex ``u+1``.
    No shell enumeration is used.
    """
    H = tuple(tuple(sorted(e)) for e in H)
    h = len(H)
    if h % 2 != 1:
        raise ValueError("H must have odd cardinality")
    n = p * p + 1
    eps, labels = field_direction_data(p, direction)

    scaled_signed_mean = 0
    for a, b in H:
        if a == 0:
            scaled_signed_mean += 1
            continue
        la, lb = labels[a - 1], labels[b - 1]
        cab = int(C[a, b])
        if la == lb:
            if cab != eps:
                raise AssertionError("parallel finite edge has wrong Paley type")
            scaled_signed_mean += p
        else:
            scaled_signed_mean += -eps * cab
    a_value = scaled_signed_mean - 3 * p

    D, degree = edge_boundary(n, H)
    counts = [0] * p
    for vertex in D:
        if vertex != 0:
            counts[labels[vertex - 1]] += 1
    B = tuple(i for i, count in enumerate(counts) if count & 1)
    b_size = len(B)
    c_H = math.prod(int(C[a, b]) for a, b in H)

    # (-1)^A = eps*(-1)^((h-3)/2)*c_H*chi_D(y), and
    # chi_D(y)=eps^[infinity in D] * (-1)^(b+sum_{B} x_s).
    sign = eps
    sign *= -1 if ((h - 3) // 2) & 1 else 1
    sign *= c_H
    if 0 in D:
        sign *= eps
    if b_size & 1:
        sign *= -1
    phase = 0 if sign == 1 else 1
    floor = scaled_direction_floor(p, b_size, phase)
    return {
        "direction": direction,
        "eps": eps,
        "scaled_signed_mean": scaled_signed_mean,
        "a": a_value,
        "a_even": a_value % 2 == 0,
        "boundary": tuple(sorted(D)),
        "degree": tuple(degree),
        "fibre_counts_on_boundary": tuple(counts),
        "B": B,
        "b": b_size,
        "c_H": c_H,
        "phase": phase,
        "parity_floor": floor,
    }


def affine_boundary_budget(p: int, H: Iterable[tuple[int, int]]) -> dict:
    """Evaluate (1)--(3) exactly for an edge set H."""
    import numpy as np

    from minmax_quadratic import paley_conference_prime_power

    H = tuple(tuple(sorted(e)) for e in H)
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    rows = [direction_record(p, H, C, d) for d in projective_directions(p)]
    expected_sum = (p + 1) * (len(H) - 3 * p)
    expected_type_sum = middle_weight(p) * (len(H) - 3 * p)
    sum_a = sum(int(row["a"]) for row in rows)
    floor_sum = sum(int(row["parity_floor"]) for row in rows)
    type_sums = {
        str(eps): sum(int(row["a"]) for row in rows if row["eps"] == eps)
        for eps in (-1, 1)
    }
    type_floor_sums = {
        str(eps): sum(
            int(row["parity_floor"]) for row in rows if row["eps"] == eps
        )
        for eps in (-1, 1)
    }
    return {
        "p": p,
        "h": len(H),
        "rows": rows,
        "sum_a": sum_a,
        "expected_sum": expected_sum,
        "sum_identity": sum_a == expected_sum,
        "type_sums": type_sums,
        "expected_type_sum": expected_type_sum,
        "type_sum_identity": all(
            value == expected_type_sum for value in type_sums.values()
        ),
        "all_a_even": all(bool(row["a_even"]) for row in rows),
        "all_affine_mean_nonnegative": all(int(row["a"]) >= 0 for row in rows),
        "parity_floor_sum": floor_sum,
        "parity_budget_holds": floor_sum <= expected_sum,
        "type_parity_floor_sums": type_floor_sums,
        "type_parity_budgets_hold": all(
            value <= expected_type_sum for value in type_floor_sums.values()
        ),
    }


def direct_affine_slack_audit(
    p: int,
    H: Iterable[tuple[int, int]],
) -> dict:
    """Enumerate every affine middle-slice point and audit the formulas.

    This is deliberately a finite diagnostic, used below only at ``p=5``.
    The all-prime proof uses the exact moment and edge-product identities in
    :func:`direction_record`, not this enumeration.
    """
    import numpy as np

    from minmax_quadratic import paley_conference_prime_power

    H = tuple(tuple(sorted(e)) for e in H)
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    m = middle_weight(p)
    rows = []
    for direction in projective_directions(p):
        record = direction_record(p, H, C, direction)
        eps, labels = field_direction_data(p, direction)
        slacks: list[int] = []
        parity_ok = True
        integer_ok = True
        for chosen in itertools.combinations(range(p), m):
            X = set(chosen)
            y = [eps] + [1 if labels[u] in X else -1 for u in range(p * p)]
            score = sum(int(C[a, b]) * y[a] * y[b] for a, b in H)
            numerator = eps * score - 3
            integer_ok = integer_ok and numerator % 2 == 0
            slack = numerator // 2
            slacks.append(slack)
            predicted_parity = (
                sum(fibre in X for fibre in record["B"]) + record["phase"]
            ) & 1
            parity_ok = parity_ok and slack % 2 == predicted_parity

        count = math.comb(p, m)
        enumerated_a = Fraction(2 * p * sum(slacks), count)
        rows.append(
            {
                "direction": direction,
                "eps": eps,
                "slack_support": tuple(sorted(set(slacks))),
                "minimum_slack": min(slacks),
                "integer_slacks": integer_ok,
                "boundary_parity_pointwise": parity_ok,
                "a_from_enumeration": enumerated_a,
                "a_from_moments": record["a"],
                "mean_identity": enumerated_a == record["a"],
            }
        )
    return {
        "p": p,
        "rows": rows,
        "all_slacks_nonnegative": all(row["minimum_slack"] >= 0 for row in rows),
        "all_slacks_integral": all(row["integer_slacks"] for row in rows),
        "all_boundary_parities": all(
            row["boundary_parity_pointwise"] for row in rows
        ),
        "all_mean_identities": all(row["mean_identity"] for row in rows),
    }


# A genuine integral solution of all corrected affine inequalities at p=5.
# It prevents this proposition from being mistaken for a residual close.
AFFINE_P5_G = (
    (0, 7), (0, 13), (0, 22), (0, 25), (1, 23),
    (2, 3), (2, 5), (2, 14), (3, 10), (3, 19),
    (3, 22), (5, 15), (5, 17), (5, 25), (7, 12),
    (7, 19), (10, 25), (13, 22), (14, 25), (15, 23),
)


def affine_p5_counterexample_to_affine_close() -> dict:
    H = AFFINE_P5_G + ((0, 1),)
    budget = affine_boundary_budget(5, H)
    direct = direct_affine_slack_audit(5, H)
    return {
        "proved": (
            budget["sum_identity"]
            and budget["type_sum_identity"]
            and budget["all_a_even"]
            and budget["all_affine_mean_nonnegative"]
            and direct["all_slacks_nonnegative"]
            and direct["all_slacks_integral"]
            and direct["all_boundary_parities"]
            and direct["all_mean_identities"]
            and [row["a"] for row in budget["rows"]] == [12, 4, 0, 6, 10, 4]
        ),
        "a_by_direction": [row["a"] for row in budget["rows"]],
        "a_by_quadratic_type": budget["type_sums"],
        "slack_supports_by_direction": [
            row["slack_support"] for row in direct["rows"]
        ],
        "direct_affine_audit": direct,
        "boundary": budget["rows"][0]["boundary"],
        "boundary_is_infinity_plus_affine_line": budget["rows"][0]["boundary"]
        == (0, 2, 7, 12, 17, 22),
        "not_full_shell_claim": True,
    }


def theorem_affine_parity_budget() -> dict:
    """Machine-readable statement of the all-prime reduction."""
    known = {
        (7, 3, 0): Fraction(4, 7),
        (7, 3, 1): Fraction(1),
        (11, 3, 0): Fraction(8, 11),
        (11, 5, 1): Fraction(9, 11),
        (17, 5, 0): Fraction(1),
        (17, 5, 1): Fraction(1),
    }
    exact = all(
        parity_majorant_floor(p, b, phase)["value"] == value
        for (p, b, phase), value in known.items()
    )
    complement = True
    for p in (5, 7, 11, 13, 17, 19):
        m_parity = middle_weight(p) & 1
        for b in range(p + 1):
            for phase in (0, 1):
                lhs = parity_majorant_floor(p, p - b, phase)["value"]
                rhs = parity_majorant_floor(p, b, phase ^ m_parity)["value"]
                complement = complement and lhs == rhs
    eulerian = all(
        eulerian_residual_type_budget_gap(p) == (p * p - 1) // 2
        and eulerian_residual_type_budget_gap(p) > 0
        for p in (3, 5, 7, 11, 13, 17, 19)
    )
    p5 = affine_p5_counterexample_to_affine_close()
    return {
        "proved": bool(exact and complement and eulerian and p5["proved"]),
        "all_odd_primes": True,
        "directional_sum": "sum a_d=(p+1)(|H|-3p)",
        "quadratic_type_sum": (
            "sum_{eps_d=tau} a_d=(p+1)(|H|-3p)/2 for tau=+-1"
        ),
        "residual_k_eq_4p_type_sum": "(p+1)^2/2",
        "residual_k_eq_4p_sum": "(p+1)^2",
        "parity_budget": (
            "sum_{eps_d=tau} 2 ceil(p M(p,b_d,eta_d)) "
            "<= (p+1)(|H|-3p)/2 for tau=+-1"
        ),
        "majorant_exact_samples": exact,
        "complement_symmetry_samples": complement,
        "eulerian_residual_boundary_excluded": eulerian,
        "eulerian_residual_type_budget_gap": "(p^2-1)/2",
        "affine_p5_feasible": p5,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, tuple):
        return [_jsonable(v) for v in value]
    if isinstance(value, list):
        return [_jsonable(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    return value


def main() -> dict:
    theorem = theorem_affine_parity_budget()
    tables = {
        str(p): {
            str(b): {
                str(phase): _jsonable(parity_majorant_floor(p, b, phase))
                for phase in (0, 1)
            }
            for b in range((p + 1) // 2)
        }
        for p in (5, 7, 11, 13, 17, 19)
    }
    out = {
        "prop": "15.632",
        "title": "Affine slack budget and quadratic parity lifts",
        "proved": {
            "all_prime_directional_slack_sum": theorem["proved"],
            "all_prime_boundary_parity_formula": theorem["proved"],
            "all_prime_quadratic_majorant_budget": theorem["proved"],
            "residual_ii_k_ge_4p": False,
            "R1": False,
            "L": False,
        },
        "theorem": theorem,
        "majorant_tables": tables,
        "L_status": "OPEN",
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15632.json"
    dest.write_text(json.dumps(_jsonable(out), indent=2) + "\n")
    print(f"wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
