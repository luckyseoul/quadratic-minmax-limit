#!/usr/bin/env python3
r"""Proposition 15.763 -- signed affine-alias incidence bound.

Continue the dangerous shared-maximizer setup of Proposition 15.755.  Thus
``G`` has even size, ``H=G union {e}``, ``B=C xor H`` has norm at most
``Phi-4``, and a signed common maximizer ``x`` has

    delta_eps(x)=2*p*r^2,          eps*S_H(x)=2-p*r^2,          (1)

where ``r`` is positive and odd.  Assume specifically that ``x`` is in the
affine alias family: there are

    m=(p+1)/2+r

positive parallel fibres, and flipping the union ``T_J`` of any ``r`` of
them gives a Boolean ``eps*p`` eigenvector ``w_J=x^T_J``.

Put ``a_h=eps*C_h*x_u*x_v`` for ``h={u,v} in H`` and let ``mu_h`` be the
number of alias cuts crossed by ``h``.  For every alias, the gap of ``B`` and
the odd cardinality of ``H`` give

    eps*S_H(w_J) >= 3.

Writing ``L_J=sum_(h crossing T_J) a_h``, (1) therefore gives

    -L_J >= (p*r^2+1)/2.                                  (2)

There are ``N=binom(m,r)`` aliases, and

    mu_h <= M=2*binom(m-2,r-1).                            (3)

Unlike Proposition 15.755, retain the signs in the sum of (2).  Since
``sum_h a_h=2-p*r^2``, the number of negative ``a_h`` is exactly

    N_minus=(|H|+p*r^2-2)/2.

Positive edges only decrease ``sum_h mu_h*(-a_h)``.  Consequently

    N*(p*r^2+1)/2 <= M*N_minus,

and hence

    |H| >= (p*r^2+1)m(m-1)/(2r(m-r)) - p*r^2 + 2.         (4)

The cardinality ``|H|`` is odd, so the executable bound is the least odd
integer at least the rational right side.  At ``r=1`` the right side is the
already odd integer

    (p^2+11)/4,

strictly stronger than the parity-adjusted Proposition 15.755 bound from
``(p+1)(p+3)/8``.

There is also a scoped critical-alias alternative.  The distinguished edge
has ``a_e=+1``.  If an alias cut does not cross ``e`` and has
``eps*S_H(w_J)=3``, then ``eps*S_G(w_J)=2``.  If no such alias exists, every
noncrossing alias has H-score at least five.  Summing this extra unit and
retaining the negative contribution of ``e`` cancels the unknown number of
alias cuts crossing ``e`` and yields (4) with ``p*r^2+1`` replaced by
``p*r^2+3``.

This is a proved conditional theorem for the affine-alias subfamily.  It
does not classify all minimum-shell representatives (the three-coordinate
integral-eigenvector branch of Proposition 15.755 is larger), does not align
the affine coordinates belonging to different deletions, and does not close
residual (ii), E1, or the MathOverflow limit.
"""
from __future__ import annotations

import json
from fractions import Fraction
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15721 import is_prime


ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15763.json"


def _validate_parameters(p: int, r: int) -> None:
    if not isinstance(p, int) or isinstance(p, bool) or p < 5 or not is_prime(p):
        raise ValueError("p must be an odd prime at least five")
    if not isinstance(r, int) or isinstance(r, bool) or r < 1 or r % 2 == 0:
        raise ValueError("r must be a positive odd integer")
    if (p + 1) // 2 + r > p:
        raise ValueError("the affine positive-fibre count m cannot exceed p")


def ceil_fraction(value: Fraction) -> int:
    """Return the exact integer ceiling of ``value``."""
    return -(-value.numerator // value.denominator)


def least_odd_integer_at_least(value: Fraction) -> int:
    """Return the least odd integer greater than or equal to ``value``."""
    answer = ceil_fraction(value)
    return answer if answer % 2 else answer + 1


def affine_alias_parameters(p: int, r: int) -> dict[str, int]:
    """Return the exact alias-family incidence parameters ``m,N,M``."""
    _validate_parameters(p, r)
    m = (p + 1) // 2 + r
    aliases = comb(m, r)
    between_positive = 2 * comb(m - 2, r - 1)
    positive_to_outside = comb(m - 1, r - 1)
    if positive_to_outside > between_positive:
        raise ArithmeticError("the asserted maximum alias-cut multiplicity failed")
    return {
        "p": p,
        "r": r,
        "m": m,
        "aliases": aliases,
        "max_edge_multiplicity": between_positive,
        "positive_to_outside_multiplicity": positive_to_outside,
    }


def prop15755_unsigned_bound(p: int, r: int) -> Fraction:
    """Return Proposition 15.755's sign-discarded rational lower bound."""
    row = affine_alias_parameters(p, r)
    m = row["m"]
    return Fraction((p * r * r + 1) * m * (m - 1), 4 * r * (m - r))


def signed_affine_alias_bound(p: int, r: int) -> Fraction:
    """Return the signed rational right side of (4)."""
    row = affine_alias_parameters(p, r)
    m = row["m"]
    return (
        Fraction((p * r * r + 1) * m * (m - 1), 2 * r * (m - r))
        - p * r * r
        + 2
    )


def critical_alias_avoidance_bound(p: int, r: int) -> Fraction:
    """Bound forced if no noncrossing H-score-three alias exists."""
    row = affine_alias_parameters(p, r)
    m = row["m"]
    return (
        Fraction((p * r * r + 3) * m * (m - 1), 2 * r * (m - r))
        - p * r * r
        + 2
    )


def signed_incidence_feasibility(p: int, r: int, h: int) -> dict[str, object]:
    """Audit the exact signed-count implication for a proposed odd ``|H|``.

    This is not a graph search.  It checks the numerical inequality obtained
    after (2)--(3) and the forced signed row sum in (1).
    """
    row = affine_alias_parameters(p, r)
    if not isinstance(h, int) or isinstance(h, bool) or h < 1 or h % 2 == 0:
        raise ValueError("h=|H| must be a positive odd integer")
    signed_sum = 2 - p * r * r
    positive = (h + signed_sum) // 2
    negative = (h - signed_sum) // 2
    integral_counts = (
        h + signed_sum == 2 * positive and h - signed_sum == 2 * negative
    )
    nonnegative_counts = positive >= 0 and negative >= 0
    cut_floor = Fraction(p * r * r + 1, 2)
    lower_total = row["aliases"] * cut_floor
    upper_total = row["max_edge_multiplicity"] * negative
    incidence_compatible = bool(
        integral_counts and nonnegative_counts and lower_total <= upper_total
    )
    bound_satisfied = Fraction(h) >= signed_affine_alias_bound(p, r)
    if incidence_compatible != (nonnegative_counts and bound_satisfied):
        raise ArithmeticError("signed incidence and closed-form bound disagree")
    return {
        **row,
        "H_size": h,
        "signed_H_score_at_x": signed_sum,
        "positive_edge_count": positive,
        "negative_edge_count": negative,
        "cut_floor": str(cut_floor),
        "summed_cut_lower": str(lower_total),
        "negative_edge_upper": str(upper_total),
        "incidence_compatible": incidence_compatible,
        "bound_satisfied": bound_satisfied,
        "proved_equivalent": True,
    }


def affine_alias_bound_ledger(p: int, r: int) -> dict[str, object]:
    """Return the exact Proposition 15.763 arithmetic ledger."""
    row = affine_alias_parameters(p, r)
    old = prop15755_unsigned_bound(p, r)
    signed = signed_affine_alias_bound(p, r)
    critical = critical_alias_avoidance_bound(p, r)
    old_odd = least_odd_integer_at_least(old)
    signed_odd = least_odd_integer_at_least(signed)
    critical_odd = least_odd_integer_at_least(critical)
    row_sum_bound = p * r * r - 2
    effective = max(old_odd, signed_odd, row_sum_bound)
    critical_effective = max(old_odd, critical_odd, row_sum_bound)
    identity = signed == 2 * old - p * r * r + 2
    if not identity:
        raise ArithmeticError("signed/unsigned comparison identity failed")
    r_one_formula = Fraction(p * p + 11, 4) if r == 1 else None
    if r == 1 and not (signed == r_one_formula and signed.denominator == 1):
        raise ArithmeticError("the r=1 signed closed form failed")
    return {
        **row,
        "boolean_defect": 2 * p * r * r,
        "signed_H_score_at_active_x": 2 - p * r * r,
        "cut_floor_per_alias": str(Fraction(p * r * r + 1, 2)),
        "forced_negative_edge_count": "(|H|+p*r^2-2)/2",
        "prop15755_unsigned_rational_bound": str(old),
        "prop15755_unsigned_odd_integer_bound": old_odd,
        "signed_rational_bound": str(signed),
        "signed_odd_integer_bound": signed_odd,
        "absolute_row_sum_bound": row_sum_bound,
        "effective_odd_integer_bound": effective,
        "signed_strictly_improves_parity_adjusted_15755": signed_odd > old_odd,
        "signed_comparison_identity": "B_signed=2*B_15755-p*r^2+2",
        "critical_alias_alternative": {
            "avoidance_rational_bound": str(critical),
            "avoidance_odd_integer_bound": critical_odd,
            "effective_avoidance_odd_integer_bound": critical_effective,
            "conclusion_below_bound": (
                "some affine alias has eps*S_H=3 and does not cross e, "
                "hence eps*S_(H\\{e})=2"
            ),
        },
        "r_equals_one_closed_form": str(r_one_formula) if r == 1 else None,
        "proved": True,
    }


def theorem_record() -> dict[str, object]:
    return {
        "prop": "15.763",
        "title": "Signed affine-alias incidence bound",
        "status": "PROVED CONDITIONAL AFFINE-ALIAS THEOREM",
        "formula": (
            "|H| >= oddceil((p*r^2+1)m(m-1)/(2r(m-r))-p*r^2+2), "
            "m=(p+1)/2+r"
        ),
        "r_equals_one": {
            "signed_bound": "(p^2+11)/4",
            "already_an_odd_integer": True,
            "prop15755_bound": "(p+1)(p+3)/8",
        },
        "proved": {
            "signed_affine_alias_incidence_bound": True,
            "all_admissible_odd_r": True,
            "critical_internal_alias_or_stronger_size_bound": True,
            "all_minimum_shell_representatives_are_affine": False,
            "deletion_alias_coordinates_are_common": False,
            "residual_ii_closed": False,
            "e1_closed_general": False,
            "L": False,
        },
        "exact_audits": [
            affine_alias_bound_ledger(5, 1),
            affine_alias_bound_ledger(11, 1),
            affine_alias_bound_ledger(11, 3),
            affine_alias_bound_ledger(11, 5),
        ],
        "not_claimed": [
            "classification of the full defect-2p shell",
            "a common affine coordinate system across all deletions",
            "exclusion of the nonaffine three-coordinate eigenvector branch",
            "closure of residual (ii), E1, or the MathOverflow limit",
        ],
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    out = theorem_record()
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print("Prop. 15.763 signed affine-alias incidence: proved")
    print("  r=1: |H| >= (p^2+11)/4")
    print("  residual (ii): OPEN")
    print("wrote", EV)
    return out


if __name__ == "__main__":
    main()
