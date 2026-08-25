#!/usr/bin/env python3
"""Prop 15.636 — complete third Paley-dual shell for every ``p>=11``.

Proposition 15.635 reduced equality at scaled norm ``2(p+1)`` to the
projected positive point-pair orbit except for one zero-profile case.  In
that case, with ``m=(p-1)/2``, the positive profile has one repeated root
and ``m-2`` other roots, while the negative profile has ``m`` distinct
roots.  Their power sums agree through degree ``m-1``.

Let ``A,B`` be the corresponding monic degree-``m`` root polynomials.
Newton identities give ``A-B=constant``.  Their distinct roots cover all
of ``F_p`` except two points ``u,v``, and the repeated root is ``alpha``.
After the affine normalization ``u=0, v=1, alpha=lambda``, one obtains

    ((A+B)/2)^2
      = x^(p-1) + a(x+...+x^(p-2)) + b,       a=1-lambda,

where ``a`` and ``c=a-1`` are both nonzero.

Reverse the monic square root and compare it with the formal series

    R(y)=sqrt((1+c*y)/(1-y)).

Uniqueness of the square root at ``y=0`` forces the coefficients of ``R``
in degrees ``m+1,...,p-2`` to vanish.  Below degree ``p`` those coefficients
are the coefficients of

    K(y)=(1+c*y)^(m+1)(1-y)^m.

Thus ``K=U+q*y^(p-1)+r*y^p``, with ``deg(U)<=m`` and ``q!=0``.  For
``1<=j<=m``, its Hasse derivative of order ``m+j`` is the monomial

    D^[m+j] K = q*(-1)^(m+j)*y^(m-j).

Evaluating it at the two roots ``1`` and ``-1/c`` and comparing with the
product formula for ``K`` gives

    c^(m-j) = -1/(2j).                              (1)

Taking ``j=m-1`` and ``j=m-2`` forces ``c=1/3`` and ``c^2=1/5``.
Hence ``1/9=1/5``, or ``p|4``, impossible for ``p>=11``.  The exceptional
profile cannot occur.

Consequently, for every odd prime ``p>=11``, the complete third shell is

    { +/- P(e_i + C_ij e_j) : i != j },

with signed size ``p^2(p^2+1)``.  Its degree-four harmonic operator is the
negative scalar already computed in Proposition 15.635.  This completes a
dual shell, not the full theta tail, so R1 remains open.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15635 import (
    third_norm,
    third_pair_harmonic_coefficient,
    third_pair_signed_count,
    third_scaled_norm,
)

ROOT = Path(__file__).resolve().parents[1]


def mod_inverse(value: int, p: int) -> int:
    return pow(value % p, -1, p)


def hasse_forced_power(p: int, j: int) -> tuple[int, int]:
    """Return ``(m-j, -1/(2j))`` from equation (1), modulo ``p``."""
    m = (p - 1) // 2
    if not 1 <= j <= m:
        raise ValueError("need 1 <= j <= (p-1)/2")
    return m - j, (-mod_inverse(2 * j, p)) % p


def exceptional_profile_contradiction(p: int) -> dict:
    """Exact finite-field certificate for the uniform Hasse contradiction.

    The proof above is symbolic in ``p``.  These residues expose its final
    two equations and make denominator/sign regressions testable.
    """
    if p < 11 or p % 2 == 0:
        raise ValueError("the theorem is stated for odd p>=11")
    m = (p - 1) // 2
    exponent_one, forced_c = hasse_forced_power(p, m - 1)
    exponent_two, forced_c_sq = hasse_forced_power(p, m - 2)
    one_third = mod_inverse(3, p)
    one_fifth = mod_inverse(5, p)
    contradiction = (forced_c * forced_c - forced_c_sq) % p
    return {
        "p": p,
        "m": m,
        "j_m_minus_1": m - 1,
        "first_exponent": exponent_one,
        "forced_c": forced_c,
        "one_third": one_third,
        "j_m_minus_2": m - 2,
        "second_exponent": exponent_two,
        "forced_c_squared": forced_c_sq,
        "one_fifth": one_fifth,
        "contradiction_residue": contradiction,
        "equivalent_nonzero_numerator": 4 % p,
        "exceptional_profile_impossible": bool(
            exponent_one == 1
            and exponent_two == 2
            and forced_c == one_third
            and forced_c_sq == one_fifth
            and contradiction != 0
        ),
    }


def complete_third_shell_theorem(
    primes: tuple[int, ...] = (11, 13, 17, 19, 23, 29, 31),
) -> dict:
    rows = {}
    ok = True
    for p in primes:
        obstruction = exceptional_profile_contradiction(p)
        row_ok = (
            obstruction["exceptional_profile_impossible"]
            and third_scaled_norm(p) == 2 * (p + 1)
            and third_norm(p) == Fraction(p + 1, p)
            and third_pair_signed_count(p) == p * p * (p * p + 1)
            and third_pair_harmonic_coefficient(p)
            == -Fraction(p * p + 4 * p - 3, 4 * (p * p + 5))
            and third_pair_harmonic_coefficient(p) < 0
        )
        rows[str(p)] = {
            "third_scaled_norm": third_scaled_norm(p),
            "third_norm": str(third_norm(p)),
            "complete_signed_count": third_pair_signed_count(p),
            "harmonic_scalar": str(third_pair_harmonic_coefficient(p)),
            "hasse_certificate": obstruction,
            "checks": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "scope": "complete third dual shell for every odd prime p>=11",
        "classification": "+/-P(e_i+C_ij e_j), i!=j",
        "uniform_obstruction": (
            "The sole extra equality profile would force c=1/3 and "
            "c^2=1/5, hence p divides 4."
        ),
        "rows": rows,
    }


def main() -> dict:
    theorem = complete_third_shell_theorem()
    out = {
        "prop": "15.636",
        "title": "Complete third Paley-dual shell",
        "proved": {
            "complete_third_shell_all_odd_p_ge_11": theorem["proved"],
            "third_shell_harmonic_scalar_all_odd_p_ge_11": theorem["proved"],
            "R1": False,
            "global_QVAR": False,
        },
        "theorem": theorem,
        "remaining_obstruction": (
            "The fourth and later norm-parity-twisted harmonic shells are "
            "not controlled uniformly; R1 is not implied by three shells."
        ),
        "L_status": "OPEN",
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15636.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Prop 15.636 complete third shell: {theorem['proved']}")
    print(f"  wrote {dest}")
    return out


if __name__ == "__main__":
    main()
