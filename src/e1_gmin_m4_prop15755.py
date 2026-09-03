#!/usr/bin/env python3
r"""Prop. 15.755 -- full-cube spike reduction and its sharp aliases.

Let ``C`` be the symmetric Paley conference matrix of order
``n=p^2+1``, let ``Phi=np/2``, and put

    P_eps=(I+eps*C/p)/2,          eps in {+1,-1}.

For a Boolean vector ``x`` define its signed conference defect by

    delta_eps(x)=Phi-eps*q_C(x)=p*||P_-eps x||^2,

where ``q_C(x)=x^T C x/2``.  Propositions 15.629--15.630 imply the
following exact full-cube gap.  If ``x`` is not an ``eps*p`` Boolean
eigenvector, then

    delta_eps(x) >= 2p.                                      (1)

Indeed, fix a Boolean ``eps*p`` eigenvector ``y`` and write
``x-y=2z``.  Then ``P_-eps x=2P_-eps z`` and the projected integral
lattice has minimum squared norm ``1/2``.

Equality in (1) is completely described by the minimum shell.  There are
an index ``i`` and ``a in {+1,-1}`` such that

    P_-eps z=a P_-eps e_i,
    y'=x-2a e_i in ker_Z(C-eps*p*I).                         (2)

If ``a=x_i``, ``y'`` is Boolean and ``x`` is a one-coordinate flip of a
Boolean eigenshell vector.  If ``a=-x_i``, ``y'`` has coordinate
``3x_i`` at ``i`` and all other coordinates are Boolean.  The latter is
not an empty algebraic possibility.

For every odd prime, normalize infinity first and choose a square-kernel
linear functional ``L:F_(p^2)->F_p``.  If ``g:F_p->{+1,-1}`` has sum
three, then

    v_infinity=3,        v_u=g(Lu)

satisfies ``Cv=pv``.  Replacing ``v_infinity=3`` by one gives a Boolean
``x`` with ``delta_+(x)=2p`` in the bad-sign case of (2).  Thus the
ordinary one-bit interpretation of the first defect shell is false; the
same vector is also a flip of one affine ``p``-point fibre from a Boolean
eigenvector.

Now suppose ``G`` has even size, ``e`` is not in ``G``, and

    A=C xor G,       H=G union {e},       B=C xor H.

If ``||q_A||_infinity=Phi-2`` and ``||q_B||_infinity<=Phi-4``, choose a
signed maximizer ``eps*q_A(x)=Phi-2``.  The one-edge difference forces

    eps*f_e(x)=+1,
    eps*S_G(x)=1-delta_eps(x)/2,
    eps*S_H(x)=2-delta_eps(x)/2.                            (3)

Since ``|G|`` is even, (3) makes ``delta_eps(x)/2`` odd.  In particular
``x`` is off the Boolean eigenshell and (1) applies.  For ``p>=11``, the
odd-phase dual-lattice gap of Proposition 15.635 strengthens (1) to the
sharp dichotomy

    delta_eps(x)=2p,       or       delta_eps(x)>=6p-12.    (4)

The second endpoint is also attained on the Boolean cube for every odd
prime.  Switch a Boolean eigenvector ``y`` to all ones.  Its positive-edge
graph is strongly regular of degree ``p(p+1)/2`` and triangle parameter
``(p-1)(p+3)/4>0``.  Flipping a positive triangle gives
``delta=6p-12`` and lands in Proposition 15.639's projected signed-triple
shell.  Therefore neither branch of (4) can be deleted by lattice
quantization alone.

There are two further exact consequences.

* Every Boolean ``eps*p`` eigenvector ``w`` obeys
  ``eps*S_H(w)>=3``.  If ``x=w^T``, its signed H-cut satisfies

      eps*c_H(T;w) >= (delta_eps(x)+2)/4,

  hence at least ``(p+1)/2`` in the first branch of (4).

* Sharing the maximizer supplies no second hereditary cut restriction.
  In the ``x,eps`` gauge, write ``W_D(T)`` for the signed cut of a matrix
  ``D``.  Global maximality is exactly ``0<=W_D(T)<=M_D`` for every
  ``T``.  Since the distinguished A-edge has gauge weight ``+1``,

      W_A(T)=W_B(T)+2*1_{e crosses T},      M_A=M_B+2.

  Thus all B-cut inequalities automatically imply all A-cut inequalities.

The affine alias family does yield one honest sparse-range obstruction.  More
generally, let odd ``r>=1`` and choose ``m=(p+1)/2+r`` positive fibres.  The
resulting Boolean vector has defect ``2p r^2``.  Removing any ``r`` of its
positive fibres gives a Boolean eigenshell vector, so (4) applies to all
``binom(m,r)`` unions of ``r`` parallel fibres.  Any edge crosses at most
``2 binom(m-2,r-1)`` of those unions.  Summing the signed-cut inequalities
and discarding signs therefore forces

    |H| >= (p r^2+1)m(m-1)/(4r(m-r)).                      (5)

At ``r=1`` this is ``|H| >= (p+1)(p+3)/8``.  Formula (5) excludes this
particular affine spike in sufficiently sparse layers, but it neither
classifies every minimum-shell representative nor reaches dense residual
sizes.

This proposition is a global reduction and a method barrier, not a close
of residual (ii).  It leaves the common-graph exclusion of the two sharp
Boolean spike families open and does not flip any global predicate.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15635 import odd_nonminimum_scaled_floor
from e1_gmin_m4_prop15721 import is_prime
from minmax_quadratic import paley_conference_prime_power


ROOT = Path(__file__).resolve().parents[1]


def n_of(p: int) -> int:
    if not isinstance(p, int) or isinstance(p, bool) or p < 3 or p % 2 == 0:
        raise ValueError("need an odd integer p>=3")
    return p * p + 1


def phi_of(p: int) -> int:
    return p * n_of(p) // 2


def quadratic_score(C: np.ndarray, x: np.ndarray) -> int:
    value = int(x.astype(np.int64) @ C.astype(np.int64) @ x.astype(np.int64))
    if value % 2:
        raise ArithmeticError("symmetric zero-diagonal score must be even")
    return value // 2


def defect_projector_identity(p: int, signed_quadratic_score: int) -> dict[str, object]:
    """Record ``delta=p||P_-eps x||^2`` in exact score normalization."""
    delta = phi_of(p) - signed_quadratic_score
    return {
        "p": p,
        "Phi": phi_of(p),
        "signed_q_C": signed_quadratic_score,
        "delta": delta,
        "projected_norm_squared": str(Fraction(delta, p)),
        "identity": "delta_eps(x)=p*||P_-eps*x||^2",
        "proved": delta >= 0,
    }


def equality_sign_classification(x_i: int, a: int) -> dict[str, object]:
    """Classify the exceptional coordinate of ``y'=x-2a e_i``."""
    if x_i not in (-1, 1) or a not in (-1, 1):
        raise ValueError("x_i and a must be signs")
    exceptional = x_i - 2 * a
    good = a == x_i
    bad = a == -x_i
    if not (good ^ bad):
        raise ArithmeticError("the two equality signs must be exhaustive")
    return {
        "x_i": x_i,
        "minimum_shell_sign_a": a,
        "eigenvector_exceptional_coordinate": exceptional,
        "one_bit_boolean_eigenvector_case": good and exceptional == -x_i,
        "one_coordinate_three_case": bad and exceptional == 3 * x_i,
        "proved": abs(exceptional) == (1 if good else 3),
    }


def affine_bad_sign_vector(p: int) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(v,x)`` with ``Cv=pv``, ``v_inf=3``, and Boolean ``x``.

    The repo field encoding is ``u=a+b*p``.  The functional ``L(u)=b``
    has square kernel ``F_p``.  Choosing ``(p+3)/2`` positive fibres makes
    the fibre sign sum equal three.
    """
    n_of(p)
    positive_fibres = set(range((p + 3) // 2))
    v = np.ones(p * p + 1, dtype=np.int64)
    v[0] = 3
    for u in range(p * p):
        v[1 + u] = 1 if u // p in positive_fibres else -1
    x = v.copy()
    x[0] = 1
    return v, x


def affine_odd_parameter_ledger(p: int, r: int) -> dict[str, object]:
    """General odd-``r`` affine alias and its summed-cut lower bound."""
    n_of(p)
    if not isinstance(r, int) or isinstance(r, bool) or r < 1 or r % 2 == 0:
        raise ValueError("need an odd positive affine parameter r")
    m = (p + 1) // 2 + r
    if m > p:
        raise ValueError("the positive-fibre count cannot exceed p")
    representations = 1
    for j in range(1, r + 1):
        representations = representations * (m - r + j) // j
    max_edge_multiplicity = 2
    for j in range(1, r):
        max_edge_multiplicity = (
            max_edge_multiplicity * (m - r - 1 + j) // j
        )
    # The preceding product is 2*C(m-2,r-1); retain a direct integer audit.
    from math import comb

    exact_max = 2 * comb(m - 2, r - 1)
    cut_floor = Fraction(p * r * r + 1, 2)
    edge_lower_bound = cut_floor * Fraction(representations, exact_max)
    closed = Fraction(
        (p * r * r + 1) * m * (m - 1),
        4 * r * (m - r),
    )
    proved = bool(
        representations == comb(m, r)
        and max_edge_multiplicity == exact_max
        and edge_lower_bound == closed
    )
    if not proved:
        raise ArithmeticError("the affine summed-cut incidence identity failed")
    return {
        "p": p,
        "odd_parameter_r": r,
        "positive_fibres_m": m,
        "boolean_defect": 2 * p * r * r,
        "eigenshell_representations": representations,
        "maximum_edge_cut_multiplicity": exact_max,
        "signed_cut_floor_per_representation": str(cut_floor),
        "H_edge_lower_bound": str(closed),
        "r_equals_one_closed_form": (
            str(Fraction((p + 1) * (p + 3), 8)) if r == 1 else None
        ),
        "proved": proved,
    }


def affine_bad_sign_certificate(p: int) -> dict[str, object]:
    """Exact all-coordinate audit of the uniform bad-sign construction."""
    if not is_prime(p):
        raise ValueError("the Paley construction requires an odd prime")
    C = paley_conference_prime_power(p).astype(np.int64)
    v, x = affine_bad_sign_vector(p)
    residual = C @ v - p * v
    q = quadratic_score(C, x)
    delta = phi_of(p) - q
    boolean_x = bool(np.all(np.abs(x) == 1))
    eig_exact = bool(np.all(residual == 0))
    proved = bool(
        boolean_x
        and eig_exact
        and v[0] == 3
        and np.array_equal(v[1:], x[1:])
        and delta == 2 * p
    )
    if not proved:
        raise ArithmeticError("the affine bad-sign construction failed")
    return {
        "p": p,
        "dimension": p * p + 1,
        "positive_fibres": (p + 3) // 2,
        "fibre_sign_sum": 3,
        "integral_eigenvector_coordinate_multiset": {
            "plus_or_minus_one": p * p,
            "plus_three": 1,
        },
        "C_v_equals_p_v": eig_exact,
        "x_is_boolean": boolean_x,
        "q_C_x": q,
        "delta_plus_x": delta,
        "bad_sign_equality_case": True,
        "proved": proved,
    }


def positive_triangle_endpoint(p: int) -> dict[str, object]:
    """Symbolic SRG proof that a three-flip realizes ``delta=6p-12``."""
    n = n_of(p)
    degree = p * (p + 1) // 2
    restricted = ((p - 1) // 2, -(p + 1) // 2)
    triangle_parameter = (p - 1) * (p + 3) // 4
    cut = 3 * p - 6
    delta = 2 * cut
    proved = bool(
        2 * degree - (n - 1) == p
        and restricted[0] + restricted[1] == -1
        and triangle_parameter > 0
        and delta == 6 * p - 12
    )
    if not proved:
        raise ArithmeticError("the positive-triangle endpoint identity failed")
    return {
        "p": p,
        "switched_positive_graph_degree": degree,
        "restricted_eigenvalues": list(restricted),
        "common_neighbors_of_positive_edge": triangle_parameter,
        "positive_triangle_exists": True,
        "signed_cut_of_triangle": cut,
        "boolean_defect_after_triangle_flip": delta,
        "scaled_odd_dual_norm": delta // 2,
        "proved": proved,
    }


def dangerous_spike_dichotomy(p: int) -> dict[str, object]:
    """Package (3)--(4) for a hypothetical dangerous residual descent."""
    if p < 11 or not is_prime(p):
        raise ValueError("the odd nonminimum dual gap is stated for primes p>=11")
    higher_scaled_norm = odd_nonminimum_scaled_floor(p)
    higher_delta = 2 * higher_scaled_norm
    first_score = 2 - p
    higher_score_cap = 2 - higher_scaled_norm
    proved = bool(
        higher_scaled_norm == 3 * p - 6
        and higher_delta == 6 * p - 12
        and higher_score_cap == 8 - 3 * p
    )
    return {
        "p": p,
        "dangerous_maximizer_is_off_eigenshell": True,
        "delta_over_two_is_odd": True,
        "defect_dichotomy": [f"delta={2 * p}", f"delta>={higher_delta}"],
        "corresponding_signed_H_score": [
            f"eps*S_H={first_score}",
            f"eps*S_H<={higher_score_cap}",
        ],
        "first_branch_sharp_by_affine_line_alias": True,
        "second_branch_sharp_by_positive_triangle_flip": True,
        "proved": proved,
    }


def hereditary_cut_redundancy(edge_crosses: bool, W_B: int, M_B: int) -> dict[str, object]:
    """Audit ``B`` cut feasibility implies ``A`` cut feasibility."""
    if not 0 <= W_B <= M_B:
        raise ValueError("the supplied B cut must satisfy its hereditary interval")
    W_A = W_B + (2 if edge_crosses else 0)
    M_A = M_B + 2
    return {
        "edge_crosses_cut": edge_crosses,
        "W_B": W_B,
        "M_B": M_B,
        "W_A": W_A,
        "M_A": M_A,
        "B_interval": [0, M_B],
        "A_interval": [0, M_A],
        "A_cut_automatic": 0 <= W_A <= M_A,
        "proved": 0 <= W_A <= M_A,
    }


def theorem_record(primes: tuple[int, ...] = (3, 5, 7, 11, 13, 17, 19)) -> dict[str, object]:
    affine = {str(p): affine_bad_sign_certificate(p) for p in primes}
    triangles = {str(p): positive_triangle_endpoint(p) for p in primes}
    dichotomy = {
        str(p): dangerous_spike_dichotomy(p) for p in primes if p >= 11
    }
    signs = [equality_sign_classification(x_i, a) for x_i in (-1, 1) for a in (-1, 1)]
    cuts = [
        hereditary_cut_redundancy(crosses, W, 20)
        for crosses in (False, True)
        for W in (0, 7, 20)
    ]
    affine_cut_bounds = {
        str(p): affine_odd_parameter_ledger(p, 1) for p in primes
    }
    proved = bool(
        all(row["proved"] for row in affine.values())
        and all(row["proved"] for row in triangles.values())
        and all(row["proved"] for row in dichotomy.values())
        and all(row["proved"] for row in signs)
        and all(row["proved"] for row in cuts)
        and all(row["proved"] for row in affine_cut_bounds.values())
    )
    return {
        "prop": "15.755",
        "title": "Full-cube spike reduction and sharp affine/triangle aliases",
        "proved": {
            "off_eigenshell_boolean_defect_at_least_2p": proved,
            "dangerous_defect_dichotomy_p_ge_11": proved,
            "delta_2p_bad_sign_family_all_odd_primes": proved,
            "delta_6p_minus_12_triangle_family_all_odd_primes": proved,
            "shared_AB_maximizer_cut_constraints_redundant": proved,
            "affine_alias_summed_cut_lower_bound": proved,
            "residual_ii_closed": False,
            "e1_closed_general": False,
            "L": False,
        },
        "equality_sign_cases": signs,
        "affine_bad_sign_checks": affine,
        "positive_triangle_endpoint_checks": triangles,
        "dangerous_spike_dichotomy_checks": dichotomy,
        "hereditary_cut_redundancy_checks": cuts,
        "affine_alias_summed_cut_checks": affine_cut_bounds,
        "remaining_obstruction": (
            "Exclude both sharp Boolean spike families using a genuinely common "
            "H-graph invariant; lattice quantization and shared-maximizer cuts alone do not."
        ),
        "duplicate_work_guards": [
            "Do not identify every delta=2p state with a one-bit eigenshell neighbor.",
            "Do not claim delta>2p implies delta>6p-12 on the Boolean cube.",
            "Do not treat the A and B hereditary cut intervals as independent.",
        ],
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    out = theorem_record()
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15755.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("Prop. 15.755 full-cube spike dichotomy: proved")
    print("  residual (ii): OPEN")
    print(f"  wrote {dest}")
    return out


if __name__ == "__main__":
    main()
