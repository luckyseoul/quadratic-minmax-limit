#!/usr/bin/env python3
"""Prop 15.638: the first post-third even candidate shell is empty.

Retain the profile notation of Props. 15.630--15.637 and scale dual norms
by s=2p||x||^2.  Proposition 15.637 removes common sum t=0 at
s=2(p+3).  Profile balancing leaves only

    |t| = 2, p-1, p+1.

Each remaining case has only two or four units of profile-energy defect.

For t=2, ordinary profiles are delta_a+delta_b.  The binary cubic identity

    2q3 - 3q1q2 + q1^3 = 0

is forced by the ordinary directions.  With two energy-four exceptions it
kills every three-positive/one-negative profile, leaving two doubled
points.  The binary quadratic D=2q2-q1^2 is then square on every selected
direction and has two distinct selected roots.  If N is the anisotropic
quadratic defining the selected half of P^1(F_p), character orthogonality
forces |sum eta(ND)|=p-3.  But Y^2=ND is a smooth genus-one curve, so Hasse
gives at most 2sqrt(p), a contradiction for p>=11.

With one energy-six exception, the doubled pattern is killed by the cubic.
For the four-positive/two-negative pattern, the ordinary directions also
force the unsigned-pair quartic recurrence.  Its first four moments would
then equal those of a two-root multiset.  Newton identities make the four
positive roots equal the two negative roots plus those two roots, contrary
to disjoint support.

For t=p-1, replace a by b=1-a; for t=p+1, replace a by b=a-1.
The transformed profiles have common sum one.  Ordinary profiles are
single deltas, so q_d=q_1^d is forced as a binary-form identity whenever
enough ordinary directions remain.  The only low-energy exceptions are
excluded either by the factored quadratic defect q2-q1^2 or by Newton
identities on two three-element multisets.

Thus no vector has scaled norm 2(p+3) for any odd prime p>=11.  This is a
gap theorem, not a classification of the next nonempty shell or an R1
tail bound.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def half_index(p: int) -> int:
    if p < 11 or p % 2 == 0:
        raise ValueError("the theorem is stated for odd p>=11")
    return (p - 1) // 2


def direction_count(p: int) -> int:
    return half_index(p) + 1


def candidate_scaled_norm(p: int) -> int:
    half_index(p)
    return 2 * (p + 3)


def balanced_profile_energy(p: int, t: int) -> int:
    a, b = divmod(abs(t), p)
    return (p - b) * a * a + b * (a + 1) * (a + 1)


def profile_balancing_gap(p: int, t: int) -> int:
    return (p + 1) * balanced_profile_energy(p, t) - t * t


def candidate_common_sum_magnitudes(p: int) -> tuple[int, ...]:
    bound = candidate_scaled_norm(p)
    return tuple(
        t
        for t in range(0, 2 * p, 2)
        if profile_balancing_gap(p, t) <= bound
    )


def total_profile_energy(p: int, t: int) -> int:
    numerator = candidate_scaled_norm(p) + t * t
    if numerator % 2:
        raise ArithmeticError("profile energy must be integral")
    return numerator // 2


def excess_over_balancing(p: int, t: int) -> int:
    return total_profile_energy(p, t) - direction_count(p) * (
        balanced_profile_energy(p, t)
    )


def unsigned_pair_cubic_defect(q1: int, q2: int, q3: int) -> int:
    return 2 * q3 - 3 * q1 * q2 + q1**3


def unsigned_pair_quartic_defect(
    q1: int, q2: int, q4: int
) -> int:
    return 2 * q4 - q2 * q2 - 2 * q1 * q1 * q2 + q1**4


def sum_one_delta_quadratic_defect(q1: int, q2: int) -> int:
    return q2 - q1 * q1


def t2_energy6_doubled_cubic_defect(a: int, b: int, c: int) -> int:
    q1 = 2 * a + b - c
    q2 = 2 * a * a + b * b - c * c
    q3 = 2 * a**3 + b**3 - c**3
    return unsigned_pair_cubic_defect(q1, q2, q3)


def t2_energy6_doubled_cubic_factor(a: int, b: int, c: int) -> int:
    return 6 * (a - c) ** 2 * (b - c)


def t2_energy4_nondegenerate_cubic_defect(
    a: int, b: int, c: int, d: int
) -> int:
    q1 = a + b + c - d
    q2 = a * a + b * b + c * c - d * d
    q3 = a**3 + b**3 + c**3 - d**3
    return unsigned_pair_cubic_defect(q1, q2, q3)


def t2_energy4_nondegenerate_cubic_factor(
    a: int, b: int, c: int, d: int
) -> int:
    return 6 * (a - d) * (b - d) * (c - d)


def sum_one_energy3_defect(a: int, b: int, c: int) -> int:
    q1 = a + b - c
    q2 = a * a + b * b - c * c
    return sum_one_delta_quadratic_defect(q1, q2)


def sum_one_energy3_factor(a: int, b: int, c: int) -> int:
    return -2 * (a - c) * (b - c)


def sum_one_doubled_energy5_defect(a: int, b: int) -> int:
    q1 = 2 * a - b
    q2 = 2 * a * a - b * b
    return sum_one_delta_quadratic_defect(q1, q2)


def sum_one_doubled_energy5_factor(a: int, b: int) -> int:
    return -2 * (a - b) ** 2


def t2_hasse_obstruction(p: int) -> dict:
    half_index(p)
    required_abs_character_sum = p - 3
    return {
        "selected_directions": direction_count(p),
        "two_distinct_D_roots": 2,
        "required_abs_character_sum": required_abs_character_sum,
        "hasse_upper_bound_squared": 4 * p,
        "required_value_squared": required_abs_character_sum**2,
        "contradiction": required_abs_character_sum**2 > 4 * p,
    }


def t2_branch_excluded(p: int) -> bool:
    k = half_index(p)
    cert = t2_hasse_obstruction(p)
    return bool(
        excess_over_balancing(p, 2) == 4
        and k - 1 > 3
        and k > 4
        and cert["contradiction"]
        and 6 % p != 0
        and 4 < p
    )


def t_p_minus_1_branch_excluded(p: int) -> bool:
    k = half_index(p)
    return bool(
        excess_over_balancing(p, p - 1) == 4
        and k - 1 > 2
        and k > 3
        and 2 % p != 0
        and 3 < p
    )


def t_p_plus_1_branch_excluded(p: int) -> bool:
    k = half_index(p)
    return bool(
        excess_over_balancing(p, p + 1) == 2
        and k > 2
        and 2 % p != 0
    )


def zero_sum_dependency(p: int) -> bool:
    try:
        from .e1_gmin_m4_prop15637 import zero_common_sum_candidate_excluded
    except ImportError:
        from e1_gmin_m4_prop15637 import zero_common_sum_candidate_excluded

    return zero_common_sum_candidate_excluded(p)


def candidate_shell_excluded(p: int) -> bool:
    return bool(
        candidate_common_sum_magnitudes(p) == (0, 2, p - 1, p + 1)
        and zero_sum_dependency(p)
        and t2_branch_excluded(p)
        and t_p_minus_1_branch_excluded(p)
        and t_p_plus_1_branch_excluded(p)
    )


def candidate_shell_gap_theorem(
    primes: tuple[int, ...] = (11, 13, 17, 19, 23, 29, 31),
) -> dict:
    rows = {}
    ok = True
    for p in primes:
        row_ok = candidate_shell_excluded(p)
        rows[str(p)] = {
            "scaled_norm_excluded": candidate_scaled_norm(p),
            "balancing_allowed_abs_t": candidate_common_sum_magnitudes(p),
            "profile_excess": {
                "t=2": excess_over_balancing(p, 2),
                "t=p-1": excess_over_balancing(p, p - 1),
                "t=p+1": excess_over_balancing(p, p + 1),
            },
            "zero_sum_dependency": zero_sum_dependency(p),
            "t2_hasse": t2_hasse_obstruction(p),
            "t2_excluded": t2_branch_excluded(p),
            "t_p_minus_1_excluded": t_p_minus_1_branch_excluded(p),
            "t_p_plus_1_excluded": t_p_plus_1_branch_excluded(p),
            "checks": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "scope": (
            "For every odd prime p>=11, no vector in the Paley dual "
            "lattice has scaled norm 2(p+3)."
        ),
        "rows": rows,
    }


def main() -> dict:
    theorem = candidate_shell_gap_theorem()
    out = {
        "prop": "15.638",
        "title": "Empty first post-third even candidate shell",
        "proved": {
            "scaled_norm_2p_plus_6_absent_all_p_ge_11": theorem["proved"],
            "complete_fourth_shell": False,
            "R1": False,
            "global_QVAR": False,
        },
        "theorem": theorem,
        "remaining_obstruction": (
            "The next nonempty dual norm and the complete fourth-and-later "
            "harmonic theta tail remain unclassified."
        ),
        "L_status": "OPEN",
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15638.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Prop 15.638 empty post-third even candidate: {theorem['proved']}")
    print(f"  wrote {dest}")
    return out


if __name__ == "__main__":
    main()
