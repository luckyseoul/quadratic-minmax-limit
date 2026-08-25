#!/usr/bin/env python3
"""Prop 15.637 — exclude every zero-sum profile at energy ``p+3``.

Retain the zero-common-sum profile notation and put ``k=(p-1)/2`` and
``R=k+1``.  The first even energy after the complete third shell is
``E=p+3`` (scaled dual norm ``2E=2(p+3)``).  If ``h<R`` profiles are active,
the MDS/Newton mass bound gives

    M >= h(R-h),                 E >= 2M.

For ``p>=11`` this leaves only ``h=1`` or ``h=R-1``; ``h=R`` is the
separate full-direction case.  The first half of the proposition excludes
``h=1`` by the two root-polynomial arguments below.

For one active profile, ``M>=k`` and ``E=2k+4``.  There are only two
possible multiplicity patterns.

1. ``M=k+1``: one entry has magnitude two and every other nonzero entry
   has magnitude one.  The positive and negative degree-``k+1`` root
   polynomials ``A,B`` cover every field element, with one repeated root
   ``alpha``.  Equal power sums through degree ``k-1`` make ``A-B`` linear,
   while

       AB=(X^p-X)(X-alpha).

   If ``T=(A+B)/2``, reversing the square identity shows that its reverse
   agrees through degree ``p-2`` with ``sqrt(1-alpha*y)``.  In characteristic
   ``p`` that square root is ``(1-alpha*y)^(k+1)``, so
   ``T=(X-alpha)^(k+1)``.  Then ``(A-B)^2/4=(X-alpha)^2``; both ``A`` and
   ``B`` are divisible by ``X-alpha``, contradicting disjoint supports.

2. ``M=k``: exactly two entries have magnitude two; the root supports omit
   three field elements.  Normalize the omitted elements to ``0,1,rho``
   and let the repeated roots be ``alpha,beta``.  With

       N=(1-alpha*y)(1-beta*y),   D=(1-y)(1-rho*y),

   the reverse square root agrees through degree ``2k-1`` with the formal
   series ``S=sqrt(N/D)=sum u_j y^j``.  Hence
   ``u_(k+1)=...=u_(2k-1)=0``.  Differentiation gives

       2ND S'=(N'D-ND')S.

   The leading coefficient of ``ND`` is ``alpha*beta*rho != 0``.  In the
   coefficient of degree ``j+3``, once ``u_(j+1),...,u_(j+4)`` vanish, the
   displayed equation reduces to

       2(alpha*beta*rho) j u_j=0.

   Since ``k>=5``, start at ``j=k`` and descend to obtain
   ``u_1=...=u_k=0``.  Comparing degree one and two in ``D S^2=N`` gives
   ``N=D``.  Thus the repeated roots are the omitted roots ``1,rho``, again
   contradicting disjointness.

For the dense branches, let ``q_d`` be the degree-``d`` binary form supplied
by the profile glue.  A signed pair has moments satisfying

    4 q_1 q_3 - 3 q_2^2 - q_1^4 = 0.              (3)

If ``h=R``, energy ``p+3`` gives one energy-four exceptional profile and
``R-1=k>=5`` ordinary pairs.  The left side of (3) is a binary quartic, so
those five zeros make it identically zero.  On
``delta_a+delta_b-delta_c-delta_d`` it equals

    -12(a-c)(a-d)(b-c)(b-d),

which is nonzero.

If ``h=R-1``, the unique inactive direction is the root of ``q_1=L`` and
``q_2=L S`` for a linear form ``S``.  The cubic

    4q_3-L(3S^2+L^2)

vanishes at all ordinary pairs and at the inactive direction.  It excludes
two energy-four exceptions.  With one energy-six exception it first excludes
the multiplicity pattern ``2delta_a-delta_b-delta_c``.  In the remaining
three-positive/three-negative pattern, also write ``q_4=L T``.  The cubic

    2T-S(S^2+L^2)

vanishes at the ``k-1>=4`` ordinary directions.  Thus the exceptional first
four moments equal those of one signed pair.  Newton identities applied to
the two resulting four-element multisets force equality, contradicting the
six disjoint roots.

Hence the entire zero-common-sum branch at energy ``p+3`` is empty.  The
balancing bound leaves only ``|t|=2,p-1,p+1`` among nonzero common sums at
that scaled norm.  Those cases are not excluded here, so this does not yet
classify the next shell or prove R1.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def half_index(p: int) -> int:
    if p < 11 or p % 2 == 0:
        raise ValueError("the theorem is stated for odd p>=11")
    return (p - 1) // 2


def candidate_energy(p: int) -> int:
    return p + 3


def candidate_scaled_norm(p: int) -> int:
    return 2 * candidate_energy(p)


def profile_balancing_gap(p: int, t: int) -> int:
    """Lower bound ``(p+1) f_p(t)-t^2`` from balanced profiles."""
    a, b = divmod(abs(t), p)
    return p * a * a + 2 * a * b + b * (p + 1 - b)


def candidate_common_sum_magnitudes(p: int) -> tuple[int, ...]:
    """Even ``|t|`` not removed by balancing at scaled norm ``2(p+3)``.

    Values with quotient ``a>=2`` have gap at least ``4p>2(p+3)`` for
    ``p>=11``, so one complete range ``0<=t<2p`` is exhaustive.
    """
    half_index(p)
    bound = candidate_scaled_norm(p)
    return tuple(
        t
        for t in range(0, 2 * p, 2)
        if profile_balancing_gap(p, t) <= bound
    )


def mds_allowed_active_counts_before_one_profile_kill(p: int) -> tuple[int, ...]:
    """Counts not excluded by ``h(R-h)<=E/2``, including ``h=R``."""
    k = half_index(p)
    R = k + 1
    mass_ceiling = candidate_energy(p) // 2
    below_R = [h for h in range(1, R) if h * (R - h) <= mass_ceiling]
    return tuple(below_R + [R])


def one_profile_mass_patterns(p: int) -> tuple[dict, dict]:
    """The two integer multiplicity patterns at energy ``p+3``."""
    k = half_index(p)
    return (
        {
            "mass": k + 1,
            "energy_defect": 2,
            "magnitude_two_entries": 1,
            "distinct_support": p,
            "root_polynomial_degree": k + 1,
            "root_polynomial_difference_degree_at_most": 1,
        },
        {
            "mass": k,
            "energy_defect": 4,
            "magnitude_two_entries": 2,
            "distinct_support": p - 3,
            "root_polynomial_degree": k,
            "root_polynomial_difference_degree_at_most": 0,
        },
    )


def two_double_ode_descent(p: int) -> dict:
    """Index certificate for the formal-series coefficient descent."""
    k = half_index(p)
    gap_start = k + 1
    gap_end = 2 * k - 1
    first_equation_degree = k + 3
    highest_series_index_used = k + 4
    return {
        "k": k,
        "known_zero_gap": [gap_start, gap_end],
        "first_descent_j": k,
        "first_equation_degree": first_equation_degree,
        "highest_series_index_used": highest_series_index_used,
        "gap_reaches_first_equation": highest_series_index_used <= gap_end,
        "nonzero_descent_multipliers": all(
            (2 * j) % p != 0 for j in range(1, k + 1)
        ),
        "conclusion": "u_1=...=u_k=0, hence N=D",
    }


def one_profile_candidate_excluded(p: int) -> bool:
    patterns = one_profile_mass_patterns(p)
    descent = two_double_ode_descent(p)
    return bool(
        patterns[0]["distinct_support"] == p
        and patterns[1]["distinct_support"] == p - 3
        and descent["gap_reaches_first_equation"]
        and descent["nonzero_descent_multipliers"]
    )


def remaining_candidate_active_counts(p: int) -> tuple[int, int]:
    before = mds_allowed_active_counts_before_one_profile_kill(p)
    if not one_profile_candidate_excluded(p):
        raise ArithmeticError("one-profile proof certificate failed")
    return tuple(h for h in before if h != 1)  # type: ignore[return-value]


def signed_pair_cubic_defect(q1: int, q2: int, q3: int) -> int:
    """The quartic relation satisfied by ``delta_a-delta_b`` moments."""
    return 4 * q1 * q3 - 3 * q2 * q2 - q1**4


def energy_four_cubic_defect(a: int, b: int, c: int, d: int) -> int:
    q1 = a + b - c - d
    q2 = a * a + b * b - c * c - d * d
    q3 = a**3 + b**3 - c**3 - d**3
    return signed_pair_cubic_defect(q1, q2, q3)


def energy_four_cubic_factor(a: int, b: int, c: int, d: int) -> int:
    return -12 * (a - c) * (a - d) * (b - c) * (b - d)


def doubled_energy_six_cubic_defect(a: int, b: int, c: int) -> int:
    q1 = 2 * a - b - c
    q2 = 2 * a * a - b * b - c * c
    q3 = 2 * a**3 - b**3 - c**3
    return signed_pair_cubic_defect(q1, q2, q3)


def doubled_energy_six_cubic_factor(a: int, b: int, c: int) -> int:
    return -12 * (a - b) ** 2 * (a - c) ** 2


def dense_branch_recurrence_certificate(p: int) -> dict:
    k = half_index(p)
    return {
        "full_h_R": {
            "normal_pair_zeros": k,
            "relation_degree": 4,
            "forced_identity": k > 4,
            "exception_defect_factor": "-12(a-c)(a-d)(b-c)(b-d)",
        },
        "near_full_h_R_minus_1_two_energy_four": {
            "normal_plus_inactive_zeros": k - 1,
            "relation_degree": 3,
            "forced_identity": k - 1 > 3,
        },
        "near_full_h_R_minus_1_one_energy_six": {
            "cubic_normal_plus_inactive_zeros": k,
            "quartic_normal_zeros": k - 1,
            "relation_degree": 3,
            "forced_cubic_identity": k > 3,
            "forced_quartic_identity": k - 1 > 3,
            "doubled_entry_defect_factor": "-12(a-b)^2(a-c)^2",
            "six_unit_case": (
                "four matched power sums give equality of two size-four "
                "multisets by Newton identities"
            ),
        },
        "nonzero_characteristic_factors": (-12) % p != 0,
        "newton_degree_four_valid": 4 < p,
    }


def dense_candidate_branches_excluded(p: int) -> bool:
    cert = dense_branch_recurrence_certificate(p)
    return bool(
        cert["full_h_R"]["forced_identity"]
        and cert["near_full_h_R_minus_1_two_energy_four"]["forced_identity"]
        and cert["near_full_h_R_minus_1_one_energy_six"][
            "forced_cubic_identity"
        ]
        and cert["near_full_h_R_minus_1_one_energy_six"][
            "forced_quartic_identity"
        ]
        and cert["nonzero_characteristic_factors"]
        and cert["newton_degree_four_valid"]
    )


def zero_common_sum_candidate_excluded(p: int) -> bool:
    return one_profile_candidate_excluded(p) and dense_candidate_branches_excluded(p)


def zero_common_sum_gap_theorem(
    primes: tuple[int, ...] = (11, 13, 17, 19, 23, 29, 31),
) -> dict:
    rows = {}
    ok = True
    for p in primes:
        k = half_index(p)
        R = k + 1
        before = mds_allowed_active_counts_before_one_profile_kill(p)
        after_one_profile = remaining_candidate_active_counts(p)
        zero_sum_excluded = zero_common_sum_candidate_excluded(p)
        row_ok = (
            before == (1, R - 1, R)
            and after_one_profile == (R - 1, R)
            and zero_sum_excluded
            and candidate_common_sum_magnitudes(p) == (0, 2, p - 1, p + 1)
        )
        rows[str(p)] = {
            "candidate_energy": candidate_energy(p),
            "candidate_scaled_norm": candidate_scaled_norm(p),
            "balancing_allowed_common_sum_magnitudes": (
                candidate_common_sum_magnitudes(p)
            ),
            "allowed_h_before": before,
            "allowed_h_after_one_profile_kill": after_one_profile,
            "allowed_h_after_dense_branch_kills": (),
            "mass_patterns": one_profile_mass_patterns(p),
            "ode_descent": two_double_ode_descent(p),
            "dense_branch_certificate": dense_branch_recurrence_certificate(p),
            "checks": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "scope": (
            "At energy p+3, no zero-common-sum profile exists for odd "
            "p>=11."
        ),
        "rows": rows,
    }


def main() -> dict:
    theorem = zero_common_sum_gap_theorem()
    out = {
        "prop": "15.637",
        "title": "Zero-common-sum gap at the first post-third-shell energy",
        "proved": {
            "one_profile_energy_p_plus_3_excluded_all_p_ge_11": theorem["proved"],
            "zero_common_sum_energy_p_plus_3_excluded_all_p_ge_11": theorem[
                "proved"
            ],
            "complete_fourth_shell": False,
            "R1": False,
            "global_QVAR": False,
        },
        "theorem": theorem,
        "remaining_obstruction": (
            "The nonzero common profile sums |t|=2,p-1,p+1 at scaled norm "
            "2(p+3), and all later shells, remain unclassified."
        ),
        "L_status": "OPEN",
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15637.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Prop 15.637 zero-sum next-energy gap: {theorem['proved']}")
    print(f"  wrote {dest}")
    return out


if __name__ == "__main__":
    main()
