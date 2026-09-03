#!/usr/bin/env python3
r"""Proposition 15.762 -- integral cube gap for symmetric conferences.

Let ``p>=3`` be odd, ``n=p^2+1``, and let ``C`` be a symmetric conference
matrix of order ``n``.  With the repository normalization

    Q_C(x)=(1/2)x^t Cx,       Phi_sph=pn/2,

put ``delta=Phi_sph-Q_C(x)`` and ``z=(C-pI)x/2`` for a Boolean vector
``x``.  Then ``z`` is integral, ``Cz=-pz``, and

    ||z||^2=p delta.                                      (1)

After switching by ``x``, write ``D=diag(x) C diag(x)`` and
``w=diag(x)z``.  The graph with Seidel matrix ``D`` has all degrees of one
parity: the off-diagonal entries of ``D^2=p^2I`` give

    0=n-2-2(d_i+d_j)+4A_ij+4(A^2)_ij,

and ``n-2=p^2-1`` is divisible by eight.  Since
``d_i=(p^2-p)/2-w_i``, all coordinates of ``w`` have one parity.

This parity, (1), and the elementary integral-eigenvector bound

    ||v||_1 >= (p+1)||v||_infinity,       Dv=-pv,          (2)

exclude ``delta=2,4`` for every odd ``p>=3`` and ``delta=6`` for every
``p>=5``.  Hence for ``p>=5`` either ``Cx=px`` or ``delta>=8``.  Applying
the same theorem to ``-C`` proves that a symmetric conference matrix with
neither Boolean eigenshell satisfies

    Phi(C) <= pn/2-8.                                     (3)

At the first arithmetically possible gap ``delta=8``, necessarily
``w=2v`` where ``Dv=-pv``, ``v`` has exactly ``2p`` unit coordinates,
and ``sum(v)=-4``.  Its signs therefore split ``p-2`` positive versus
``p+2`` negative entries.

This is a universal conference/cube theorem and a counterexample-search
criterion.  It does not construct a nonregularizable conference matrix,
does not show that the ``delta=8`` shell is nonempty, and does not close
residual (ii), E1, or the MathOverflow limit.
"""
from __future__ import annotations

from typing import Sequence


def _validate_p(p: int) -> None:
    if p < 3 or p % 2 == 0:
        raise ValueError("p must be odd and at least three")


def order(p: int) -> int:
    _validate_p(p)
    return p * p + 1


def spherical_ceiling(p: int) -> int:
    """Return ``pn/2`` in the repository's half-quadratic normalization."""
    return p * order(p) // 2


def projected_norm_squared(p: int, q_value: int) -> int:
    """Return the forced value ``||(C-pI)x/2||^2=p(Phi_sph-Q)``."""
    delta = spherical_ceiling(p) - q_value
    if delta < 0:
        raise ValueError("q_value exceeds the spherical ceiling")
    return p * delta


def integral_eigenvector_norm_floor(p: int) -> int:
    """The consequence ``||v||^2>=p+1`` of the coordinate/l1 bound."""
    _validate_p(p)
    return p + 1


def low_gap_exclusion_certificate(p: int) -> dict[str, object]:
    """Return the exact arithmetic excluding gaps 2, 4, and (for p>=5) 6.

    This records the uniform proof branches; it is not a prime census.
    """
    n = order(p)
    excluded: dict[int, str] = {
        2: (
            "even w contradicts ||w||^2=2p mod 4; odd w has "
            "||w||^2>=n>2p"
        ),
        4: (
            "odd w has norm 2 mod 4, not 4p; even w=2v would have "
            "||v||^2=p<p+1"
        ),
    }
    if p >= 5:
        excluded[6] = (
            "even w contradicts ||w||^2=6p mod 4; odd w has n>6p for "
            "p>=7, while p=5 fails modulo 8"
        )
    checks = {
        "n_minus_two_divisible_by_eight": (n - 2) % 8 == 0,
        "ceiling_is_odd": spherical_ceiling(p) % 2 == 1,
        "gap_two_odd_branch_too_large": n > 2 * p,
        "gap_four_even_branch_below_eigenvector_floor": p < p + 1,
        "gap_four_odd_branch_mod_four_mismatch": n % 4 != (4 * p) % 4,
        "gap_six_uniform_for_p_ge_five": (
            p < 5
            or (p == 5 and n % 8 != (6 * p) % 8)
            or (p >= 7 and n > 6 * p)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError("conference cube-gap arithmetic changed")
    return {
        "p": p,
        "n": n,
        "all_cube_gaps_are_even": True,
        "all_switched_projection_coordinates_have_one_parity": True,
        "excluded_positive_gaps": sorted(excluded),
        "reasons": {str(gap): reason for gap, reason in excluded.items()},
        "checks": checks,
        "proved": True,
    }


def first_possible_shell_certificate(p: int) -> dict[str, object]:
    """State the exact necessary form at the first possible gap eight."""
    _validate_p(p)
    if p < 5:
        raise ValueError("the gap-eight first-shell statement requires p>=5")
    return {
        "gap": 8,
        "projection_z_norm_squared": 8 * p,
        "switched_projection": "w=2v",
        "v_eigen_equation": "Dv=-pv",
        "v_norm_squared": 2 * p,
        "v_coordinate_alphabet": [-1, 0, 1],
        "v_support": 2 * p,
        "v_sum": -4,
        "v_positive_count": p - 2,
        "v_negative_count": p + 2,
        "existence_claimed": False,
        "proved_necessary_form": True,
    }


def exact_projection_audit(
    matrix: Sequence[Sequence[int]], x: Sequence[int], p: int
) -> dict[str, object]:
    """Audit (1) and the switched parity statement for one exact input.

    The routine is a fail-when-wrong verifier only.  The theorem is the
    symbolic argument recorded above, not an enumeration using this helper.
    """
    n = order(p)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("matrix has the wrong order")
    if len(x) != n or any(value not in (-1, 1) for value in x):
        raise ValueError("x must be a Boolean vector of the conference order")
    if any(matrix[i][i] != 0 for i in range(n)):
        raise ValueError("matrix diagonal must vanish")
    if any(
        matrix[i][j] != matrix[j][i] or matrix[i][j] not in (-1, 1)
        for i in range(n)
        for j in range(i + 1, n)
    ):
        raise ValueError("matrix must be symmetric with off-diagonal signs")

    cx = [sum(matrix[i][j] * x[j] for j in range(n)) for i in range(n)]
    if any((cx[i] - p * x[i]) % 2 for i in range(n)):
        raise ArithmeticError("the projected numerator is not even")
    z = [(cx[i] - p * x[i]) // 2 for i in range(n)]
    cz = [sum(matrix[i][j] * z[j] for j in range(n)) for i in range(n)]
    q_value = sum(x[i] * cx[i] for i in range(n)) // 2
    delta = spherical_ceiling(p) - q_value
    w = [x[i] * z[i] for i in range(n)]
    return {
        "q_value": q_value,
        "gap": delta,
        "z": z,
        "w": w,
        "z_is_minus_p_eigenvector": all(cz[i] == -p * z[i] for i in range(n)),
        "sum_w_identity": sum(w) == -delta,
        "norm_identity": sum(value * value for value in z) == p * delta,
        "one_coordinate_parity": len({value % 2 for value in w}) == 1,
    }


def theorem_record() -> dict[str, object]:
    return {
        "prop": "15.762",
        "title": "Integral cube gap for square-order symmetric conferences",
        "status": "PROVED UNIVERSAL GAP AND FIRST-POSSIBLE-SHELL NECESSARY FORM",
        "normalization": {
            "Q": "(1/2)x^T Cx",
            "spherical_ceiling": "p(p^2+1)/2",
            "z": "(C-pI)x/2",
            "z_norm_squared": "p(Phi_sph-Q)",
            "norm_2p_means_gap": 2,
            "gap_4_means_z_norm_squared": "4p",
        },
        "proved": {
            "switched_graph_degrees_have_one_parity": True,
            "cube_gaps_2_and_4_absent_for_every_odd_p_at_least_3": True,
            "cube_gap_6_absent_for_every_odd_p_at_least_5": True,
            "non_eigen_boolean_gap_at_least_8_for_p_at_least_5": True,
            "first_possible_gap_8_sparse_form_is_necessary": True,
            "first_possible_gap_8_is_attained": False,
            "nonregularizable_conference_constructed": False,
            "residual_ii_closed": False,
            "e1_closed_general": False,
            "L": False,
        },
        "counterexample_criterion": (
            "If a symmetric conference matrix of order p^2+1 has neither a "
            "Boolean +p nor -p eigenvector, then Phi(C)<=p(p^2+1)/2-8."
        ),
        "first_possible_shell": first_possible_shell_certificate(5),
        "not_claimed": [
            "existence of a nonregularizable symmetric conference matrix of square order",
            "attainment of the gap-eight shell",
            "a residual-(ii) exclusion or counterexample",
            "E1 or the MathOverflow limit",
        ],
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    out = theorem_record()
    print("Prop. 15.762 conference cube gap: proved")
    print("  non-eigen Boolean gap >=8 for every odd p>=5")
    print("  residual (ii): OPEN")
    return out


if __name__ == "__main__":
    main()
