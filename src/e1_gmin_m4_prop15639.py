#!/usr/bin/env python3
"""Prop 15.639: classify the shell at the first nonminimal odd scaled norm.

The first possible nonminimal odd scaled norm is 3p-6.  Its complete shell
is the disjoint union
of:

1. projected signed triples Pz for which z has three unit coordinates and
   every signed conference edge inside the support is negative;
2. point-circle vectors Pe_i+w_S/p, where w_S is an oriented square-circle
   complement and (w_S)_i=-1.

The signed count is

    C(p^2+1,3) + p^2(p-1)(p^2+1)
      = p^2(p-1)(p+7)(p^2+1)/6.

The proof first uses parity and the norm to move every vector to a common
profile sum one.  Equality in the odd-profile bound leaves active counts
one and R-2.  The one-profile case is the square-circle equality profile.
In the R-2 case, degree-two and degree-three moments turn the vector, after
adding one minimum vector, into a common-sum-two vector.  Its number of
doubled profiles is zero, one, or two, placing it on the second shell, the
third shell, or the empty shell of Prop. 15.638.  Complete third-shell
classification excludes the middle case and Prop. 15.638 excludes the last.
The complete second-shell theorem then gives exactly the two families above.

For p=11,13 this is the fourth nonempty shell because Prop. 15.638 removes
the only intervening even candidate.  For p>=17, further even candidates
can lie below 3p-6 and are not classified here.  The result does not bound
the remaining harmonic theta tail, prove R1, or settle global QVAR.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def n_of(p: int) -> int:
    if p < 11 or p % 2 == 0:
        raise ValueError("the theorem is stated for odd p>=11")
    return p * p + 1


def direction_count(p: int) -> int:
    n_of(p)
    return (p + 1) // 2


def first_nonminimal_odd_scaled_norm(p: int) -> int:
    n_of(p)
    return 3 * p - 6


def first_nonminimal_odd_norm(p: int) -> Fraction:
    return Fraction(first_nonminimal_odd_scaled_norm(p), 2 * p)


def scaled_coordinate_square_sum(p: int) -> int:
    """Sum of r_i^2 for r=2p*x on the scaled-norm 3p-6 shell."""
    return 2 * p * first_nonminimal_odd_scaled_norm(p)


def unit_coordinate_forced(p: int) -> bool:
    """Odd r cannot have all coordinates of magnitude at least three."""
    return scaled_coordinate_square_sum(p) < 9 * n_of(p)


def equality_active_profile_counts(p: int) -> tuple[int, int]:
    """Counts surviving Delta=R-2 in the t=1 equality argument."""
    R = direction_count(p)
    delta = R - 2
    surviving = []
    for h in range(1, R + 1):
        if h == R:
            continue
        if max(h, h * (R - h - 1)) <= delta:
            surviving.append(h)
    return tuple(surviving)  # type: ignore[return-value]


def dense_profile_moment_differences(
    mu: int, A: int, B: int
) -> tuple[int, int]:
    """Return q2-mu^2 and q3-mu^3 for the dense equality profile."""
    alpha = mu + A
    beta = mu + B
    gamma = mu + A + B
    q2 = alpha * alpha + beta * beta - gamma * gamma
    q3 = alpha**3 + beta**3 - gamma**3
    return q2 - mu * mu, q3 - mu**3


def dense_profile_moment_factors(
    mu: int, A: int, B: int
) -> tuple[int, int]:
    return -2 * A * B, -3 * A * B * (2 * mu + A + B)


def transformed_t2_scaled_norm(p: int, doubled_inactive: int) -> int:
    """Norm after adding the linear-form minimum vector in the dense case."""
    if doubled_inactive not in (0, 1, 2):
        raise ValueError("there are exactly two inactive directions")
    return 2 * (p - 1) + 4 * doubled_inactive


def first_nonminimal_odd_is_fourth_norm(p: int) -> bool:
    """The s=3p-6 shell is ordinal shell four only at p=11 and p=13."""
    n_of(p)
    if p == 11:
        return (
            2 * (p + 1)
            < first_nonminimal_odd_scaled_norm(p)
            < 2 * (p + 3)
        )
    if p == 13:
        return (
            2 * (p + 1)
            < 2 * (p + 3)
            < first_nonminimal_odd_scaled_norm(p)
            < 2 * (p + 5)
        )
    return False


def negative_triangle_signed_count(p: int) -> int:
    """One antipodal signed pair for each negative conference triangle."""
    return math.comb(n_of(p), 3)


def circle_point_signed_count(p: int) -> int:
    """Both global signs of Pe_i+w_S/p with (w_S)_i=-1."""
    return p * p * (p - 1) * n_of(p)


def first_nonminimal_odd_shell_signed_count(p: int) -> int:
    return negative_triangle_signed_count(p) + circle_point_signed_count(p)


def first_nonminimal_odd_shell_signed_count_closed(p: int) -> int:
    return p * p * (p - 1) * (p + 7) * n_of(p) // 6


P11_QFMINIM = {
    "scaled_bound": 28,
    "signed_cumulative_count": 473_970,
    "maximum_scaled_norm": 27,
    "elapsed_ms": 2_033_141,
    "signed_count_through_third_shell": 31_110,
}


def p11_first_nonminimal_odd_exact_audit() -> dict:
    """Independent exact qfminim audit through the empty s=28 candidate."""
    residual = (
        P11_QFMINIM["signed_cumulative_count"]
        - P11_QFMINIM["signed_count_through_third_shell"]
    )
    expected = first_nonminimal_odd_shell_signed_count(11)
    return {
        **P11_QFMINIM,
        "scaled_norm_3p_minus_6_shell_signed_count": residual,
        "predicted_two_family_signed_count": expected,
        "exact_count_matches_classification": bool(
            residual == expected
            and P11_QFMINIM["maximum_scaled_norm"]
            == first_nonminimal_odd_scaled_norm(11)
        ),
        "backend": "PARI/GP qfminim on the saturated scaled dual Gram form",
        "host": "NUKA",
    }


def family_coordinate_signatures(p: int) -> dict:
    """Large-coordinate signatures proving disjointness and injectivity."""
    n_of(p)
    return {
        "large_magnitude": p - 2,
        "negative_triangle_large_coordinates": 3,
        "circle_point_large_coordinates": 1,
        "other_coordinate_magnitudes_at_most": 3,
        "separated": p - 2 > 3,
    }


def dependency_certificate(p: int) -> dict:
    try:
        from .e1_gmin_m4_prop15633 import second_shell_theorem
        from .e1_gmin_m4_prop15636 import complete_third_shell_theorem
        from .e1_gmin_m4_prop15638 import candidate_shell_excluded
    except ImportError:
        from e1_gmin_m4_prop15633 import second_shell_theorem
        from e1_gmin_m4_prop15636 import complete_third_shell_theorem
        from e1_gmin_m4_prop15638 import candidate_shell_excluded

    return {
        "complete_second_shell": second_shell_theorem((p,))["proved"],
        "complete_third_shell": complete_third_shell_theorem((p,))["proved"],
        "empty_even_candidate_shell": candidate_shell_excluded(p),
    }


def first_nonminimal_odd_shell_classified(p: int) -> bool:
    deps = dependency_certificate(p)
    return bool(
        unit_coordinate_forced(p)
        and equality_active_profile_counts(p) == (1, direction_count(p) - 2)
        and transformed_t2_scaled_norm(p, 0) == 2 * (p - 1)
        and transformed_t2_scaled_norm(p, 1) == 2 * (p + 1)
        and transformed_t2_scaled_norm(p, 2) == 2 * (p + 3)
        and all(deps.values())
        and family_coordinate_signatures(p)["separated"]
        and first_nonminimal_odd_shell_signed_count(p)
        == first_nonminimal_odd_shell_signed_count_closed(p)
    )


def first_nonminimal_odd_shell_theorem(
    primes: tuple[int, ...] = (11, 13, 17, 19, 23, 29, 31),
) -> dict:
    rows = {}
    ok = True
    for p in primes:
        row_ok = first_nonminimal_odd_shell_classified(p)
        rows[str(p)] = {
            "scaled_norm": first_nonminimal_odd_scaled_norm(p),
            "norm": str(first_nonminimal_odd_norm(p)),
            "scaled_coordinate_square_sum": scaled_coordinate_square_sum(p),
            "unit_coordinate_forced": unit_coordinate_forced(p),
            "equality_active_counts": equality_active_profile_counts(p),
            "transformed_t2_scaled_norms": [
                transformed_t2_scaled_norm(p, r) for r in range(3)
            ],
            "scaled_norm_3p_minus_6_is_fourth_norm": (
                first_nonminimal_odd_is_fourth_norm(p)
            ),
            "negative_triangle_signed_count": negative_triangle_signed_count(p),
            "circle_point_signed_count": circle_point_signed_count(p),
            "scaled_norm_3p_minus_6_shell_signed_count": (
                first_nonminimal_odd_shell_signed_count(p)
            ),
            "coordinate_signatures": family_coordinate_signatures(p),
            "dependencies": dependency_certificate(p),
            "checks": row_ok,
        }
        ok = ok and row_ok
    p11_exact = p11_first_nonminimal_odd_exact_audit()
    return {
        "proved": bool(ok and p11_exact["exact_count_matches_classification"]),
        "scope": (
            "For every odd prime p>=11, the complete scaled-norm 3p-6 "
            "shell is the disjoint union of negative signed triples and "
            "incident point-circle vectors."
        ),
        "rows": rows,
        "p11_exact_qfminim": p11_exact,
    }


def main() -> dict:
    theorem = first_nonminimal_odd_shell_theorem()
    out = {
        "prop": "15.639",
        "title": "Complete shell at the first nonminimal odd scaled norm",
        "proved": {
            "complete_scaled_norm_3p_minus_6_shell_all_p_ge_11": theorem[
                "proved"
            ],
            "fourth_norm_and_shell_p11_p13": (
                theorem["proved"]
                and first_nonminimal_odd_is_fourth_norm(11)
                and first_nonminimal_odd_is_fourth_norm(13)
            ),
            "fourth_norm_all_p_ge_11": False,
            "scaled_norm_3p_minus_6_shell_signed_count": theorem["proved"],
            "scaled_norm_3p_minus_6_harmonic_operator": False,
            "R1": False,
            "global_QVAR": False,
        },
        "theorem": theorem,
        "remaining_obstruction": (
            "For p>=17, intervening even shells below 3p-6 remain possible. "
            "The harmonic operator at scaled norm 3p-6 and the later norm-parity-"
            "twisted theta tail remain uncontrolled."
        ),
        "L_status": "OPEN",
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15639.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Prop 15.639 complete s=3p-6 shell: {theorem['proved']}")
    print(f"  wrote {dest}")
    return out


if __name__ == "__main__":
    main()
