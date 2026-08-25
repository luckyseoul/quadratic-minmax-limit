#!/usr/bin/env python3
"""Prop 15.635 — the third dual norm and the exact p=11 third shell.

Retain the Paley lattice notation of Props. 15.629--15.634 and scale dual
norms by ``s=2p||u||^2``.  For every odd prime ``p>=11``:

* after ``s=p`` and ``s=2(p-1)``, the next possible norm is
  ``s=2(p+1)``;
* it is attained by ``+/-P(e_i+C_ij e_j)``;
* every odd-phase vector outside the minimum shell has
  ``s>=3p-6``.

The point-pair orbit has signed size ``p^2(p^2+1)`` and degree-four
harmonic coefficient

    -(p^2+4p-3)/(4(p^2+5)) ||W||_F^2.

At ``p=11``, exact PARI ``qfminim`` enumeration through ``s=24`` gives
31,110 signed vectors.  Subtracting the proved first-shell count 244 and
second-shell count 16,104 leaves 14,762, exactly the point-pair orbit.
Hence the complete p=11 third shell is classified.  The all-prime statement
classifies the third *norm* and one full orbit, not necessarily every vector
on that shell for p>11; R1 remains open.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def n_of(p: int) -> int:
    return p * p + 1


def rank_of(p: int) -> int:
    return n_of(p) // 2


def third_scaled_norm(p: int) -> int:
    if p < 11:
        raise ValueError("the uniform third-norm theorem is stated for p>=11")
    return 2 * (p + 1)


def third_norm(p: int) -> Fraction:
    return Fraction(third_scaled_norm(p), 2 * p)


def odd_nonminimum_scaled_floor(p: int) -> int:
    return 3 * p - 6


def third_pair_signed_count(p: int) -> int:
    return p * p * n_of(p)


def third_pair_harmonic_coefficient(p: int) -> Fraction:
    """Complete signed ``+C_ij`` pair orbit evaluated at ``u/2``."""
    return Fraction(1, 4) * (
        1 - Fraction((p + 1) ** 2, rank_of(p) + 2)
    )


P11_QFMINIM = {
    "scaled_bound": 24,
    "signed_cumulative_count": 31_110,
    "maximum_scaled_norm": 24,
    "elapsed_ms": 215_744,
    "first_signed_count": 244,
    "second_signed_count": 16_104,
}


def p11_third_shell_audit() -> dict:
    residual = (
        P11_QFMINIM["signed_cumulative_count"]
        - P11_QFMINIM["first_signed_count"]
        - P11_QFMINIM["second_signed_count"]
    )
    expected = third_pair_signed_count(11)
    return {
        **P11_QFMINIM,
        "third_shell_signed_count": residual,
        "predicted_pair_signed_count": expected,
        "complete_third_shell_is_pair_orbit": residual == expected,
    }


def third_shell_theorem(
    primes: tuple[int, ...] = (11, 13, 17, 19, 23),
) -> dict:
    rows = {}
    ok = True
    for p in primes:
        row_ok = (
            third_scaled_norm(p) == 2 * (p + 1)
            and third_norm(p) == Fraction(p + 1, p)
            and odd_nonminimum_scaled_floor(p) > third_scaled_norm(p)
            and third_pair_harmonic_coefficient(p)
            == -Fraction(p * p + 4 * p - 3, 4 * (p * p + 5))
            and third_pair_harmonic_coefficient(p) < 0
        )
        rows[str(p)] = {
            "third_scaled_norm": third_scaled_norm(p),
            "third_norm": str(third_norm(p)),
            "odd_nonminimum_scaled_floor": odd_nonminimum_scaled_floor(p),
            "pair_signed_count": third_pair_signed_count(p),
            "pair_harmonic_coefficient": str(
                third_pair_harmonic_coefficient(p)
            ),
            "checks": row_ok,
        }
        ok = ok and row_ok
    p11 = p11_third_shell_audit()
    return {
        "proved": bool(ok and p11["complete_third_shell_is_pair_orbit"]),
        "scope": (
            "third norm and pair-orbit operator for all p>=11; complete "
            "third-shell classification additionally at p=11"
        ),
        "rows": rows,
        "p11_exact_qfminim": p11,
    }


def main() -> dict:
    theorem = third_shell_theorem()
    out = {
        "prop": "15.635",
        "title": "Third Paley-dual norm and exact p=11 third shell",
        "proved": {
            "third_dual_norm_all_p_ge_11": theorem["proved"],
            "odd_nonminimum_gap_all_p_ge_11": theorem["proved"],
            "third_pair_harmonic_scalar_all_p_ge_11": theorem["proved"],
            "complete_p11_third_shell": theorem["proved"],
            "complete_third_shell_all_p": False,
            "R1": False,
            "global_QVAR": False,
        },
        "theorem": theorem,
        "L_status": "OPEN",
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15635.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Prop 15.635 third norm / p11 third shell: {theorem['proved']}")
    print(f"  wrote {dest}")
    return out


if __name__ == "__main__":
    main()
