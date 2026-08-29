#!/usr/bin/env python3
"""Prop. 15.692 -- binary affine-Radon isomorphism and p=19 reduction.

Let ``A`` be the line-point incidence matrix of ``AG(2,p)``, with rows
grouped into the ``p+1`` parallel classes.  Over ``F_2`` and for odd ``p``,

    A^T A = I + J.

Indeed, a point lies on ``p+1`` lines and two distinct points lie on one
common line.  Hence on the even-weight point space ``E_P`` one has
``A^T A=I``.  Every block of ``Ax`` has parity ``wt(x)``, so ``A`` maps
``E_P`` into the direct sum of the even-weight spaces of the parallel
classes.  Both spaces have dimension ``p^2-1``.  Therefore this map is an
isomorphism, with inverse ``x=A^T r``.

For the fourteen p=19 profiles left by Proposition 15.689, this turns
boundary-profile realizability into one exact nonlinear condition: choose
the prescribed even numbers of affine lines in each direction and require
the inverse Radon word ``A^T r`` to have weight sixteen.  There are no
additional linear compatibility equations.

The fixed first two moments of the number of chosen lines through a point
also do not obstruct the profiles.  For every survivor those moments admit
an exact probability distribution supported on the even multiplicities
``{4,6,8}``.  Thus pairwise independence/second moments cannot force a
positive odd-parity density.  The endpoint itself remains open.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15688 import p19_residue_zero_profiles


ROOT = Path(__file__).resolve().parents[1]
P = 19
BOUNDARY_SIZE = 16


def affine_binary_radon_isomorphism(p: int) -> dict[str, object]:
    """Return the exact dimension and inverse identities for odd ``p``."""
    if p < 3 or p % 2 == 0:
        raise ValueError("p must be odd and at least three")
    point_count = p * p
    direction_count = p + 1
    line_count = p * direction_count
    point_even_dimension = point_count - 1
    block_even_dimension = direction_count * (p - 1)
    if point_even_dimension != block_even_dimension:
        raise ArithmeticError("Radon source and target dimensions changed")
    return {
        "p": p,
        "point_count": point_count,
        "direction_count": direction_count,
        "affine_line_count": line_count,
        "point_degree": p + 1,
        "line_size": p,
        "common_lines_for_distinct_points": 1,
        "incidence_gram_over_F2": "A^T A = I + J",
        "restriction": "even-weight point space",
        "source_dimension": point_even_dimension,
        "target": "direct sum of even-weight line blocks by direction",
        "target_dimension": block_even_dimension,
        "injective": True,
        "surjective_by_equal_dimension": True,
        "inverse": "x = A^T r",
        "proved": True,
    }


def _even_moment_distribution(
    first_moment: Fraction, second_factorial_moment: Fraction
) -> dict[int, Fraction]:
    """Match two moments on support 4,6,8, using exact arithmetic."""
    # Put a=y_6+2y_8 from E[N]-4=2a.  Subtracting 18a from
    # E[N(N-1)]-12 leaves 8y_8.
    a = (first_moment - 4) / 2
    y8 = (second_factorial_moment - 12 - 18 * a) / 8
    y6 = a - 2 * y8
    y4 = 1 - y6 - y8
    result = {4: y4, 6: y6, 8: y8}
    if any(value < 0 for value in result.values()):
        raise ArithmeticError("even moment witness stopped being nonnegative")
    if sum(result.values()) != 1:
        raise ArithmeticError("even moment witness lost unit mass")
    if sum(n * value for n, value in result.items()) != first_moment:
        raise ArithmeticError("first moment witness changed")
    if (
        sum(n * (n - 1) * value for n, value in result.items())
        != second_factorial_moment
    ):
        raise ArithmeticError("second moment witness changed")
    return result


def p19_inverse_radon_profile_reduction() -> dict[str, object]:
    """Reduce the fourteen survivors to exact inverse-weight instances."""
    profiles = [
        row
        for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) >= 16
    ]
    if len(profiles) != 14:
        raise ArithmeticError("p=19 high-slack profile count changed")

    rows = []
    for index, profile in enumerate(profiles):
        sizes = []
        for b, count in profile["global_b_profile"].items():
            sizes.extend([int(b)] * int(count))
        if len(sizes) != P + 1 or any(b % 2 for b in sizes):
            raise ArithmeticError("p=19 Radon block sizes changed")

        selected_line_count = sum(sizes)
        cross_direction_intersections = sum(
            sizes[first] * sizes[second]
            for first in range(len(sizes))
            for second in range(first)
        )
        weight_mod_four = (
            P * selected_line_count - 2 * cross_direction_intersections
        ) % 4
        first_moment = Fraction(selected_line_count, P)
        second_factorial_moment = Fraction(
            selected_line_count * selected_line_count
            - sum(value * value for value in sizes),
            P * P,
        )
        moment_witness = _even_moment_distribution(
            first_moment, second_factorial_moment
        )
        rows.append(
            {
                "profile_index": index,
                "pair_slack": int(profile["pair_slack"]),
                "phase_profiles_b": profile["phase_profiles_b"],
                "global_b_profile": profile["global_b_profile"],
                "selected_affine_lines_in_radon_word": selected_line_count,
                "cross_direction_line_intersections": (
                    cross_direction_intersections
                ),
                "inverse_weight_mod_four": weight_mod_four,
                "target_weight_mod_four": BOUNDARY_SIZE % 4,
                "mod_four_compatible": weight_mod_four
                == BOUNDARY_SIZE % 4,
                "stripe_count_first_moment": first_moment,
                "stripe_count_second_factorial_moment": (
                    second_factorial_moment
                ),
                "all_even_second_moment_witness": moment_witness,
                "second_moments_force_positive_odd_density": False,
            }
        )

    if not all(bool(row["mod_four_compatible"]) for row in rows):
        raise ArithmeticError("a p=19 profile unexpectedly failed mod four")
    return {
        "proposition": "15.692",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "profile_count": len(rows),
        "radon_theorem": affine_binary_radon_isomorphism(P),
        "exact_remaining_condition": (
            "choose line-parity blocks with the prescribed phase weights "
            "and require wt(A^T r)=16"
        ),
        "additional_linear_compatibility_conditions": 0,
        "all_profiles_pass_inverse_weight_mod_four": True,
        "pairwise_independence_or_second_moments_can_close": False,
        "rows": rows,
        "p19_second_all_finite_endpoint_closed": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "closes_type_I": False,
        "L_status": "OPEN",
        "proved": True,
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    return value


def main() -> None:
    theorem = p19_inverse_radon_profile_reduction()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15692.json"
    target.write_text(json.dumps(_jsonable(theorem), indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.692: affine binary Radon isomorphism; "
        "p=19 remainder reduced to fourteen nonlinear inverse-weight instances"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
