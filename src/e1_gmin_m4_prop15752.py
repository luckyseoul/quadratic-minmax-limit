#!/usr/bin/env python3
r"""Prop. 15.752 -- influence rigidity closes a band of residual-II shells.

For every prime ``p>=23`` there is no nonzero nonnegative integer-valued
quadratic ``B`` on ``J(p,(p+1)/2)`` with

    4p E[B] = p+9.                                      (1)

The height-at-least-two case follows from paired-cube averaging, the sharp
stabilizer bound, and Proposition 15.751's dimension-free half-mean cube
theorem.  The height-one case is Boolean; the corrected transposition-
influence argument reduces it to the same fixed four-bit catalog used by
15.751, whose fourteen density profiles miss ``(p+9)/(4p)``.

Combined with the exact isolated-chart arithmetic of Propositions
15.734--15.735, (1) closes the contiguous higher-shell band

    p = 1 mod 4:  4 <= t <= (p-9)/2,
    p = 3 mod 4:  4 <= t <= (p-7)/2,

at ``k=4p+2t``.  Branches A/C force a phase-zero cell of mass ``p+9``;
branch B forces the already-excluded mass ``p+7`` from Proposition 15.751.
No boundary-size hypothesis or finite prime/configuration census is used.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    BRANCH_P1_LAST,
    BRANCH_P3_LAST,
    baseline_coefficient_rules,
    residual_even_floor_table,
)
from e1_gmin_m4_prop15751 import (
    atomic_write_json,
    cube_half_mean_height_certificate,
    density_profile_certificate,
    exact_four_cube_catalog,
    height_at_least_two_certificate,
    height_one_junta_certificate,
    profile_density,
)


ROOT = Path(__file__).resolve().parents[1]
FIRST_NEW_LAYER = 4
LOCAL_MASS_OFFSET = 9


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _check_local_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 23
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime p>=23")


def band_maximum_t(p: int) -> int:
    """Return the last layer covered by the uniform low-floor argument."""
    _check_local_prime(p)
    q = (p - 1) // 2
    return q - 4 if p % 4 == 1 else q - 3


def _check_band_parameters(p: int, t: int) -> None:
    _check_local_prime(p)
    if (
        not isinstance(t, int)
        or isinstance(t, bool)
        or not FIRST_NEW_LAYER <= t <= band_maximum_t(p)
    ):
        raise ValueError(
            f"need 4<=t<={band_maximum_t(p)} for the supplied prime"
        )


def p_plus_nine_height_certificate(p: int) -> dict[str, object]:
    """Exclude ``max(B)>=2`` when ``4p E[B]=p+9``."""
    _check_local_prime(p)
    scaled_mass = p + LOCAL_MASS_OFFSET
    lower_height = Fraction(p - 7, 4)
    if p % 4 == 3:
        stabilizer_upper = Fraction(scaled_mass, 4)
        paired_average_upper = Fraction(scaled_mass, 2 * (p + 1))
        stabilizer_branch = "4p E[B]>=4H"
    else:
        stabilizer_upper = Fraction(scaled_mass * (p + 3), 4 * (p - 1))
        paired_average_upper = Fraction(scaled_mass, 2 * (p - 1))
        stabilizer_branch = "4p E[B]>=4(p-1)H/(p+3)"

    half_mean = cube_half_mean_height_certificate()
    proved = bool(
        half_mean["proved"]
        and lower_height > 3
        and paired_average_upper < Fraction(3, 4)
    )
    _require(proved, "the p+9 height contradiction failed")
    return {
        "p": p,
        "slice": f"J({p},{(p + 1) // 2})",
        "scaled_mass_4p_E_B": scaled_mass,
        "assumed_maximum_at_least": 2,
        "every_paired_cube_mean_at_least": "1/2",
        "paired_cube_mean_lattice": "(1/4)Z",
        "height_lower_bound": str(lower_height),
        "stabilizer_branch": stabilizer_branch,
        "stabilizer_height_upper_bound": str(stabilizer_upper),
        "paired_cube_average_upper_bound": str(paired_average_upper),
        "some_paired_cube_has_mean_exactly": "1/2",
        "half_mean_cube_maximum_upper_bound": 3,
        "contradiction": True,
        "proved": proved,
    }


def _four_bit_densities(p: int) -> list[Fraction]:
    catalog = exact_four_cube_catalog()
    return sorted(
        {
            profile_density(tuple(row["layer_counts"]), p)
            for row in catalog["profiles"]  # type: ignore[index]
        }
    )


def p_plus_nine_boolean_certificate(p: int) -> dict[str, object]:
    """Exclude the Boolean ``max(B)=1`` branch at density ``(p+9)/(4p)``."""
    _check_local_prime(p)
    mu = Fraction(p + LOCAL_MASS_OFFSET, 4 * p)
    influence_floor = Fraction((p + 1) * (p - 3), 16 * p * (p - 2))
    total_influence_upper = (p - 1) * mu * (1 - mu)
    junta_bound = Fraction(
        2
        * (p - 1)
        * (p - 2)
        * (p + LOCAL_MASS_OFFSET)
        * (3 * p - LOCAL_MASS_OFFSET),
        p * p * (p + 1) * (p - 3),
    )
    cancelled_junta_bound = Fraction(
        6 * (p - 1) * (p - 2) * (p + 9), p * p * (p + 1)
    )
    seven_gap = p**3 - 29 * p**2 + 150 * p - 108
    seven_gap_at_23 = 168
    seven_gap_derivative_at_23 = 403
    seven_gap_second_derivative_at_23 = 80
    densities = _four_bit_densities(p)
    expected = sorted(
        {
            Fraction(0),
            Fraction(1),
            Fraction(p - 3, 4 * p),
            Fraction(p + 1, 4 * p),
            Fraction(p - 1, 2 * p),
            Fraction(p + 1, 2 * p),
            Fraction(3 * p - 1, 4 * p),
            Fraction(3 * (p + 1), 4 * p),
        }
    )
    proved = bool(
        junta_bound == cancelled_junta_bound
        and seven_gap > 0
        and seven_gap_at_23 == 23**3 - 29 * 23**2 + 150 * 23 - 108
        and seven_gap_derivative_at_23 == 3 * 23**2 - 58 * 23 + 150
        and seven_gap_second_derivative_at_23 == 6 * 23 - 58
        and seven_gap_second_derivative_at_23 > 0
        and junta_bound < 7
        and densities == expected
        and Fraction(p + 1, 4 * p) < mu < Fraction(p - 1, 2 * p)
        and mu not in densities
    )
    _require(proved, "the p+9 Boolean density exclusion failed")
    return {
        "p": p,
        "complementary_slice": f"J({p},{(p - 1) // 2})",
        "target_density": str(mu),
        "relevant_pair_influence_floor": str(influence_floor),
        "total_influence_upper_bound": str(total_influence_upper),
        "largest_zero_influence_class_complement_bound": str(junta_bound),
        "seven_gap_polynomial": "p^3-29p^2+150p-108",
        "seven_gap_at_p_23": seven_gap_at_23,
        "seven_gap_derivative_at_p_23": seven_gap_derivative_at_23,
        "seven_gap_second_derivative_at_p_23": (
            seven_gap_second_derivative_at_23
        ),
        "junta_coordinates_at_most": 6,
        "cube_coordinates_actually_needed_at_most": 4,
        "possible_four_bit_density_values": [str(value) for value in densities],
        "target_strict_bracket": [
            str(Fraction(p + 1, 4 * p)),
            str(Fraction(p - 1, 2 * p)),
        ],
        "target_absent": True,
        "proved": proved,
    }


def p_plus_nine_local_exclusion(p: int) -> dict[str, object]:
    """Package the all-height local theorem used by the residual band."""
    height = p_plus_nine_height_certificate(p)
    boolean = p_plus_nine_boolean_certificate(p)
    proved = bool(height["proved"] and boolean["proved"])
    return {
        "p": p,
        "statement": (
            "no nonzero nonnegative integer-valued quadratic B on "
            f"J({p},{(p + 1) // 2}) has 4p E[B]=p+9"
        ),
        "height_at_least_two": height,
        "height_one_boolean": boolean,
        "excluded": proved,
        "finite_slice_or_prime_census_used": False,
        "fixed_four_bit_catalog_reused": True,
        "proved": proved,
    }


def p19_sharp_mechanism_witness() -> dict[str, object]:
    """Show why the local p+9 lemma genuinely stops before p=19."""
    p = 19
    layer_values = [3 - 2 * r + r * (r - 1) // 2 for r in range(5)]
    mean = Fraction(p - 5, 2 * p)
    scaled_mass = 4 * p * mean
    proved = bool(
        layer_values == [3, 1, 0, 0, 1]
        and min(layer_values) == 0
        and scaled_mass == p + 9
    )
    _require(proved, "the p=19 sharp local witness changed")
    return {
        "p": p,
        "formula": "B=3-2r+binom(r,2), r=|X intersect R|, |R|=4",
        "layer_values": layer_values,
        "mean": str(mean),
        "scaled_mass_4p_E_B": int(scaled_mass),
        "equals_p_plus_9": True,
        "is_only_a_local_quadratic_not_a_residual_graph": True,
        "proved": proved,
    }


def band_arithmetic(p: int, t: int) -> dict[str, object]:
    """Audit the isolated chart and the exact admissible layer range."""
    _check_band_parameters(p, t)
    q = (p - 1) // 2
    m = q + 1
    edge_count = 4 * p + 2 * t + 1
    isolated_gap = p * p + 1 - 2 * edge_count
    maximum_low_mean = p + 1 + 2 * t
    if p % 4 == 1:
        next_nonbaseline_phase_one_floor = 2 * p - 6
        maximum_endpoint_lift_excess = 2 * t + 2
    else:
        next_nonbaseline_phase_one_floor = 2 * p
        maximum_endpoint_lift_excess = 2 * t + 2
    proved = bool(
        isolated_gap > 0
        and maximum_low_mean < next_nonbaseline_phase_one_floor
        and maximum_endpoint_lift_excess < p - 3
        and t + 1 < m
    )
    _require(proved, "the residual-band range arithmetic failed")
    return {
        "p": p,
        "p_mod_4": p % 4,
        "layer_index_t": t,
        "closed_t_interval": [FIRST_NEW_LAYER, band_maximum_t(p)],
        "original_k": 4 * p + 2 * t,
        "H_edge_count": edge_count,
        "ambient_vertex_count": p * p + 1,
        "guaranteed_isolated_vertices": isolated_gap,
        "q": q,
        "m": m,
        "type_budget": 2 * m * (m + t),
        "phase_one_mean_form": f"a_d=2u+{p + 1}k_d",
        "phase_one_quotient_sum": "sum k_d=m+t-u",
        "maximum_low_phase_one_mean": maximum_low_mean,
        "next_nonbaseline_phase_one_floor": next_nonbaseline_phase_one_floor,
        "maximum_endpoint_baseline_lift_excess": maximum_endpoint_lift_excess,
        "sharp_integral_lift_floor": p - 3,
        "transported_infinity_degree_I": 0,
        "every_directional_b_even": True,
        "boundary_size_hypothesis_used": False,
        "proved": proved,
    }


def band_hard_residue_certificate(p: int, t: int) -> dict[str, object]:
    """Reduce every common phase-one residue to branches A/B/C."""
    arithmetic = band_arithmetic(p, t)
    q = (p - 1) // 2
    m = q + 1
    phase_one = residual_even_floor_table(p)["phase_one_floors"]
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    low_rows: list[dict[str, object]] = []
    for u in range(t + 1):
        low_mean = p + 1 + 2 * u
        available = [
            (int(b), int(floor), low_mean - int(floor))
            for b, floor in phase_one.items()
            if int(floor) <= low_mean
        ]
        allowed_b = {2, p - 1}
        _require(
            all(b in allowed_b and 0 <= excess < lift_floor for b, _, excess in available),
            "a nonbaseline low phase-one cell entered the band",
        )
        exact = [(b, floor) for b, floor, excess in available if excess == 0]
        expected_exact = [(p - 1, p + 1)] if p % 4 == 1 and u == 0 else []
        _require(exact == expected_exact, "the low-residue survivor list changed")
        low_rows.append(
            {
                "u": u,
                "low_mean": low_mean,
                "available_endpoint_cells": [list(row) for row in available],
                "exact_survivors": [list(row) for row in exact],
                "surviving_branch": BRANCH_P1_LAST if exact else None,
            }
        )

    endpoint_candidates = [
        int(value)
        for value in residual_even_floor_table(p)[
            "phase_one_cells_at_mean_p_minus_one"
        ]
    ]
    if p % 4 == 1:
        branches = [BRANCH_B2, BRANCH_P1_LAST]
        _require(endpoint_candidates == [2], "the p=1 endpoint baseline changed")
    else:
        branches = [BRANCH_B2, BRANCH_P3_LAST]
        rules = baseline_coefficient_rules(p)
        _require(
            endpoint_candidates == [2, p - 1]
            and int(rules[BRANCH_B2]["offset"])
            - int(rules[BRANCH_P3_LAST]["offset"])
            == 1,
            "the p=3 endpoint-baseline dichotomy changed",
        )

    proved = bool(
        arithmetic["proved"]
        and low_rows
        and all(
            row["surviving_branch"] is None
            or (p % 4 == 1 and row["u"] == 0)
            for row in low_rows
        )
        and t + 1 < m
    )
    _require(proved, "the hard common-residue reduction failed")
    return {
        "p": p,
        "layer_index_t": t,
        "u_0_through_t_rows": low_rows,
        "t_less_than_u_at_most_m_minus_2": (
            "every k_d>=1 but sum k_d=m+t-u<m"
        ),
        "u_equals_m_minus_1_exact_cell_count_at_least": m - (t + 1),
        "u_equals_m_minus_1_endpoint_b_candidates": endpoint_candidates,
        "equal_mean_p3_endpoint_cells_cannot_mix": True,
        "possible_branches": branches,
        "proved": proved,
    }


def _branch_constants(branch: str) -> dict[str, int | str]:
    if branch == BRANCH_B2:
        return {
            "offset": 4,
            "P": 4,
            "opposite_delta": 0,
            "Q_min": 3,
            "minimum_mean": 8,
            "next_Q": 4,
            "next_mean_offset": 9,
            "target": "p+9",
        }
    if branch == BRANCH_P1_LAST:
        return {
            "offset": 5,
            "P": 5,
            "opposite_delta": 0,
            "Q_min": 2,
            "minimum_mean": 6,
            "next_Q": 3,
            "next_mean_offset": 7,
            "target": "p+7",
        }
    if branch == BRANCH_P3_LAST:
        return {
            "offset": 3,
            "P": 3,
            "opposite_delta": 1,
            "Q_min": 4,
            "minimum_mean": 8,
            "next_Q": 5,
            "next_mean_offset": 9,
            "target": "p+9",
        }
    raise ValueError("unknown residual branch")


def band_branch_exclusion(p: int, t: int, branch: str) -> dict[str, object]:
    """Force and exclude one ``p+7`` or ``p+9`` opposite cell."""
    _check_band_parameters(p, t)
    allowed = band_hard_residue_certificate(p, t)["possible_branches"]
    if not isinstance(branch, str) or branch not in allowed:
        raise ValueError(f"branch must be one of {allowed}")
    q = (p - 1) // 2
    m = q + 1
    data = _branch_constants(branch)
    offset = int(data["offset"])
    expected_P = int(data["P"])
    opposite_delta = t + int(data["opposite_delta"])

    parameter_rows = []
    for P in range(p + 1):
        numerator = P - offset
        if numerator % q:
            continue
        rho = numerator // q
        if rho < 0:
            continue
        s = P + rho
        opposite_edges = q * (8 - s) + opposite_delta
        if opposite_edges >= 0:
            parameter_rows.append((P, rho, s, opposite_edges))
    _require(
        len(parameter_rows) == 1
        and parameter_rows[0][:3] == (expected_P, 0, expected_P),
        "the coefficient congruence no longer fixes P",
    )

    opposite_edges = int(parameter_rows[0][3])
    Q_min = int(data["Q_min"])
    minimum_mean = int(data["minimum_mean"])
    surplus = opposite_edges - m * Q_min
    minimum_forbidden = bool(0 < minimum_mean < p - 3)
    after_forbidding_minimum = surplus - m
    next_Q_forced = 0 <= after_forbidding_minimum < m
    next_mean = p + int(data["next_mean_offset"])

    phase_zero = residual_even_floor_table(p)["phase_zero_floors"]
    nonzero_b_rows = [
        (int(b), int(floor), next_mean - int(floor))
        for b, floor in phase_zero.items()
        if int(b) != 0 and int(floor) <= next_mean
    ]
    _require(
        all(
            b in {2, p - 1} and 0 < excess < p - 3
            for b, _, excess in nonzero_b_rows
        ),
        "a nonzero-b next-mass cell escaped the lift floor",
    )

    if data["target"] == "p+9":
        local_dependency = p_plus_nine_local_exclusion(p)
    else:
        _require(p % 4 == 1 and p >= 29, "p+7 dependency used out of scope")
        local_dependency = {
            "height_at_least_two": height_at_least_two_certificate(p),
            "height_one_junta": height_one_junta_certificate(p),
            "density_exclusion": density_profile_certificate(p),
            "excluded": True,
            "source": "Proposition 15.751",
            "proved": True,
        }

    proved = bool(
        minimum_forbidden
        and next_Q_forced
        and not any(excess == 0 for _, _, excess in nonzero_b_rows)
        and local_dependency["proved"]
    )
    _require(proved, "a residual-band branch survived")
    return {
        "p": p,
        "layer_index_t": t,
        "branch": branch,
        "coefficient_offset": offset,
        "feasible_P_rho_s_opposite_edges": [list(row) for row in parameter_rows],
        "forced_P": expected_P,
        "forced_rho": 0,
        "forced_s": expected_P,
        "opposite_direction_count": m,
        "opposite_parallel_count_sum": opposite_edges,
        "minimum_Q": Q_min,
        "minimum_Q_mean": minimum_mean,
        "minimum_Q_forbidden_by_phase_zero_floor_and_15_688": True,
        "surplus_above_minimum_Q": surplus,
        "surplus_after_raising_every_Q_once": after_forbidding_minimum,
        "next_Q": int(data["next_Q"]),
        "a_next_Q_direction_is_forced": next_Q_forced,
        "forced_next_scaled_mean": next_mean,
        "nonzero_b_options_at_next_mean": [list(row) for row in nonzero_b_rows],
        "next_cell_forced_to_b_zero": True,
        "next_cell_is_A_equals_2B": True,
        "local_mass_exclusion": local_dependency,
        "branch_excluded": proved,
        "proved": proved,
    }


def residual_band_exclusion(p: int, t: int) -> dict[str, object]:
    """Exclude one complete layer in the new contiguous band."""
    arithmetic = band_arithmetic(p, t)
    residues = band_hard_residue_certificate(p, t)
    branches = {
        branch: band_branch_exclusion(p, t, branch)
        for branch in residues["possible_branches"]
    }
    proved = bool(
        arithmetic["proved"]
        and residues["proved"]
        and branches
        and all(row["proved"] for row in branches.values())
    )
    _require(proved, "the residual band layer did not close")
    return {
        "p": p,
        "layer_index_t": t,
        "original_k": 4 * p + 2 * t,
        "arithmetic": arithmetic,
        "hard_residue_reduction": residues,
        "branch_exclusions": branches,
        "all_boundary_sizes_excluded": True,
        "finite_prime_or_configuration_census_used": False,
        "residual_ii_layer_excluded": True,
        "proved": proved,
    }


def proposition_15752() -> dict[str, object]:
    """Package the local theorem and representative exact band replays."""
    sample_primes = (23, 29, 31, 37, 43)
    local = {str(p): p_plus_nine_local_exclusion(p) for p in sample_primes}
    endpoints = {
        str(p): residual_band_exclusion(p, band_maximum_t(p))
        for p in sample_primes
    }
    fifth_shell = {
        str(p): residual_band_exclusion(p, FIRST_NEW_LAYER)
        for p in sample_primes
    }
    proved = bool(
        all(row["proved"] for row in local.values())
        and all(row["proved"] for row in endpoints.values())
        and all(row["proved"] for row in fifth_shell.values())
    )
    return {
        "prop": "15.752",
        "title": "Influence-rigid p+9 exclusion and contiguous residual band",
        "result_status": "proved infinite-family theorem",
        "local_theorem": (
            "for every prime p>=23, no nonzero nonnegative integral quadratic "
            "on J(p,(p+1)/2) has 4p E[B]=p+9"
        ),
        "closed_band": {
            "p_1_mod_4": "4<=t<=(p-9)/2 at k=4p+2t (p>=29)",
            "p_3_mod_4": "4<=t<=(p-7)/2 at k=4p+2t (p>=23)",
        },
        "k_eq_4p_plus_8_closed_for_every_prime_p_ge_23": proved,
        "representative_local_replays": local,
        "representative_fifth_shell_replays": fifth_shell,
        "representative_band_endpoint_replays": endpoints,
        "p19_local_threshold_witness": p19_sharp_mechanism_witness(),
        "boundary_size_hypothesis_used": False,
        "finite_prime_or_configuration_census_used": False,
        "fixed_four_bit_certificate_reused": True,
        "residual_ii_k_ge_4p_closed": False,
        "E1_closed": False,
        "quadratic_minmax_limit_closed": False,
        "remaining_scope": (
            "critical p=5,7; p=11 at k>=50; p=13,k=60,u=6 and later; "
            "p=17,k=76 and later; p=19,k=84 and later; layers beyond the "
            "displayed band for p>=23; and positive p=7,z=7"
        ),
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    if path is None:
        path = ROOT / "evidence" / "e1_gmin_m4_prop15752.json"
    atomic_write_json(path, proposition_15752())
    return path


def main() -> None:
    path = write_evidence()
    print(json.dumps({"wrote": str(path), "prop": "15.752"}, sort_keys=True))


if __name__ == "__main__":
    main()
