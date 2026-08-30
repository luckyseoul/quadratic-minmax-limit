#!/usr/bin/env python3
"""Prop 15.720 — a degree congruence excludes the required bi-tight levels.

This is the replacement for the invalid 15.55/15.167 spectral implication.
It uses the already-proved 15.272 -> 15.207 identity

    ker(G_+ + G_-) = scheme + cross

and no new finite-prime search.

If H is bi-tight of level s, |H|=sp, and d_i is its degree sequence, then

    d_i + d_j = 2ps  (mod (p^2-1)/2)       for every i != j.

Consequently all degrees have one common residue modulo (p^2-1)/2.  This
is incompatible with the handshake identity for s=2,3,4 at every prime
p>=5 (the small size comparisons are two-line residue calculations).

Only levels 2 and 3 are required by the E(1) bi-tight no-descent reductions.
Level 4 is a valid bi-tight corollary, but it does *not* exclude the one-sided
Max+/-tight level-4 covers relevant to residual (ii). Such covers exist
(15.402); its four-line family is outside residual by a 15.272/15.588 cylinder
witness. The open issue is joint residual compatibility. No QVAR, R1, or
spectral floor is used.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


def n_of(p: int) -> int:
    return p * p + 1


def degree_modulus(p: int) -> int:
    """M=(p^2-1)/2 in the pair-degree congruence."""
    return (p * p - 1) // 2


def centered_density(p: int, s: int) -> Fraction:
    """Density |H|/C(n,2) for a level-s tight cover, |H|=sp."""
    return Fraction(2 * s, n_of(p) * p)


def scheme_coordinate_from_degree(p: int, s: int, degree: int) -> Fraction:
    """The scheme coordinate g_i in B=D_g C+C D_g+X."""
    return (Fraction(degree, 1) - Fraction(s, p)) / (p * p - 1)


def degree_pair_congruence_holds(p: int, s: int, d_i: int, d_j: int) -> bool:
    """Check the necessary congruence for one pair of degrees."""
    return (d_i + d_j - 2 * p * s) % degree_modulus(p) == 0


def theorem_degree_congruence() -> dict:
    return {
        "proved": True,
        "statement": (
            "If a centered size-sp indicator lies in scheme+cross, then "
            "d_i+d_j == 2ps mod (p^2-1)/2 for every i!=j; hence all "
            "degrees have one common residue modulo (p^2-1)/2."
        ),
        "derivation": (
            "For B=C*h, Comm(B)=D_g C+C D_g with "
            "g_i=(d_i-s/p)/(p^2-1).  Off diagonal, "
            "h_ij+C_ij(CBC)_ij/p^2=2(g_i+g_j).  The integer "
            "C_ij(CBC)_ij forces p^2-1 to divide "
            "2(d_i+d_j)-4ps, which is the displayed congruence."
        ),
    }


def level_2_arithmetic_obstruction_general() -> dict:
    """Universal arithmetic after the degree congruence, s=2."""
    # M>2p iff p^2-4p-1>0.  It is 4 at p=5 and strictly increasing thereon.
    base = 5 * 5 - 4 * 5 - 1
    increasing_from_5 = 2 * 5 - 4 > 0
    return {
        "proved": base > 0 and increasing_from_5,
        "size_bound": "(p^2-1)/2 > 2p for p>=5",
        "consequence": (
            "All degrees are equal; handshake gives d=4p/(p^2+1), "
            "strictly between 0 and 1."
        ),
        "base_polynomial_p5": base,
    }


def level_3_arithmetic_obstruction_general() -> dict:
    """Universal arithmetic after the degree congruence, s=3."""
    # For p>=7, M>3p iff p^2-6p-1>0; base value at p=7 is 6.
    base = 7 * 7 - 6 * 7 - 1
    increasing_from_7 = 2 * 7 - 6 > 0
    # At p=5: M=12, n=26, sum d_i=30.  If d_i=r+12m_i,
    # then 30=26r+12 sum m_i.  Necessarily r in {0,1}; both fail mod 12.
    p5_residues = {
        0: (30 - 26 * 0) % 12,
        1: (30 - 26 * 1) % 12,
    }
    p5_empty = all(rem != 0 for rem in p5_residues.values())
    return {
        "proved": base > 0 and increasing_from_7 and p5_empty,
        "tail_size_bound": "(p^2-1)/2 > 3p for p>=7",
        "tail_consequence": (
            "All degrees are equal; handshake gives d=6p/(p^2+1), "
            "strictly between 0 and 1."
        ),
        "p5": {
            "modulus": 12,
            "vertices": 26,
            "degree_sum": 30,
            "possible_common_residues": [0, 1],
            "remainders": p5_residues,
            "empty": p5_empty,
        },
    }


def level_4_arithmetic_obstruction_general() -> dict:
    """Universal arithmetic after the degree congruence, s=4."""
    # For p>=11, M>4p iff p^2-8p-1>0; base value at p=11 is 32.
    base = 11 * 11 - 8 * 11 - 1
    increasing_from_11 = 2 * 11 - 8 > 0

    # p=5: M=12,n=26,sum d=40; p=7: M=24,n=50,sum d=56.
    # In either case n*r <= sum d forces r in {0,1}.
    small = {}
    for p in (5, 7):
        modulus = degree_modulus(p)
        n = n_of(p)
        total = 8 * p
        remainders = {r: (total - n * r) % modulus for r in (0, 1)}
        small[p] = {
            "modulus": modulus,
            "vertices": n,
            "degree_sum": total,
            "possible_common_residues": [0, 1],
            "remainders": remainders,
            "empty": all(rem != 0 for rem in remainders.values()),
        }
    return {
        "proved": base > 0 and increasing_from_11 and all(
            row["empty"] for row in small.values()
        ),
        "tail_size_bound": "(p^2-1)/2 > 4p for p>=11",
        "tail_consequence": (
            "All degrees are equal; handshake gives d=8p/(p^2+1), "
            "strictly between 0 and 1."
        ),
        "small": small,
    }


@lru_cache(maxsize=1)
def ker_gsum_eq_scheme_cross_proved_general() -> bool:
    """Import the existing 15.272 span through the 15.207 equivalence."""
    from e1_gmin_m4_prop15207 import theorem_ker_sc_reduction
    from e1_gmin_m4_prop15270 import gplus_pd_proved_general

    return bool(
        theorem_ker_sc_reduction()["proved"] and gplus_pd_proved_general()
    )


def bitight_level_obstruction(p: int, s: int) -> dict:
    """Per-prime bi-tight predicate for levels 2, 3, and the level-4 corollary."""
    valid = p >= 5 and is_prime(p) and s in (2, 3, 4)
    if s == 2:
        arithmetic = level_2_arithmetic_obstruction_general()["proved"]
    elif s == 3:
        arithmetic = level_3_arithmetic_obstruction_general()["proved"]
    elif s == 4:
        arithmetic = level_4_arithmetic_obstruction_general()["proved"]
    else:
        arithmetic = False
    ker_sc = ker_gsum_eq_scheme_cross_proved_general() if valid else False
    return {
        "p": p,
        "s": s,
        "prime_ge_5": p >= 5 and is_prime(p),
        "required_level": s in (2, 3),
        "level_4_corollary": s == 4,
        "degree_modulus": degree_modulus(p),
        "ker_gsum_eq_scheme_cross": ker_sc,
        "degree_congruence": theorem_degree_congruence()["proved"],
        "arithmetic_obstruction": arithmetic,
        "bi_tight_empty": bool(valid and ker_sc and arithmetic),
    }


def required_bitight_levels_empty_all_primes() -> bool:
    """The required bi-tight levels 2 and 3 are empty for every prime p>=5."""
    return bool(
        ker_gsum_eq_scheme_cross_proved_general()
        and theorem_degree_congruence()["proved"]
        and level_2_arithmetic_obstruction_general()["proved"]
        and level_3_arithmetic_obstruction_general()["proved"]
    )


def theorem_required_bitight_levels() -> dict:
    closed = required_bitight_levels_empty_all_primes()
    return {
        "proved": closed,
        "required_levels": [2, 3],
        "proved_bi_tight_levels": [2, 3, 4],
        "all_bi_tight_levels_claimed": False,
        "one_sided_tight_level_4_claimed": False,
        "spectral_floor_used": False,
        "qvar_or_r1_used": False,
        "ker_gsum_eq_scheme_cross": ker_gsum_eq_scheme_cross_proved_general(),
        "degree_congruence": theorem_degree_congruence(),
        "level_2": level_2_arithmetic_obstruction_general(),
        "level_3": level_3_arithmetic_obstruction_general(),
        "level_4": level_4_arithmetic_obstruction_general(),
        "consequence": (
            "The required level-2 and level-3 bi-tight alternatives in the Type-I "
            "and deep-tight no-descent reductions are empty for every prime p>=5. "
            "Bi-tight level 4 is also excluded. Generic one-sided level-4 "
            "covers exist; their joint compatibility with residual (ii) is not closed."
        ),
    }


def main() -> dict:
    return {
        "title": "Prop 15.720 degree-congruence bi-tight obstruction",
        "proved": theorem_required_bitight_levels(),
        "L_status": "OPEN",
        "remaining": [
            "multi-level Type-I bad case",
            "non-Walsh residual (ii) at even k>=4p",
        ],
    }


if __name__ == "__main__":
    import json

    print(json.dumps(main(), indent=2, default=str))
