#!/usr/bin/env python3
r"""Prop. 15.733 -- exclude the p=31 first-slack endpoint.

Continue the residual endpoint hypotheses of Proposition 15.728.  Thus
``p=31``, ``|H|=125``, the all-finite odd boundary ``D`` has size 32 and
outside pair slack ``R=10``, and one Paley type (call its sign ``eps``) has
fifteen directions of exact scaled mean 30 and one of mean 62.  At least
fourteen mean-30 directions have ``b=2``.

The mean-30 baseline contains substantially more coefficient information
than its averaged floor.  Fix a hard ``b=2`` direction, let its two odd
fibres be ``B={i,j}``, and put ``z_s=2x_s-1`` on the middle slice, so
``sum_s z_s=1``.  Exactness of the floor gives

    A=(1-x_i-x_j)^2,       eps S_H=4+z_i z_j.

Let ``I`` be the number of selected infinity edges, ``n_s`` their endpoint
counts in fibre ``s``, ``P_d`` the number of selected finite edges parallel
to the direction, and ``K_st`` the signed selected-edge sum between fibres.
Coefficient comparison modulo ``sum z_s=1`` gives an integer

    g_d=(I+P_d-4)/15

and

    K_st=eps(g_d-n_s-n_t+1_{st=B}),
    sum_{s<t}|g_d-n_s-n_t+1_{st=B}| <= 121-15g_d.       (1)

All fifteen mean-30 directions actually have ``b=2``.  Their equal means
and common Paley type force a common ``P_d=P``.  A mean-30 ``b=30`` baseline
would instead give ``15 | I+P-3``, incompatible with the ``b=2`` congruence
``15 | I+P-4`` already supplied by at least fourteen directions.

Write the common value as ``rho``.  The hard mean-62 direction has ``P+1``
parallel edges, so the number of finite selected edges of hard sign is
``16P+1``.  Hence the number of opposite-sign finite edges is

    E_opp=125-I-(16P+1)=15(8-s),       s=rho+P.          (2)

Here ``I=15rho+4-P`` is even, so ``s`` is even; nonnegativity in (2) leaves
``s in {0,2,4,6,8}``.  If an opposite-type direction has ``Q`` parallel
selected edges, its mean is

    a=30s-208+32Q,       sum_opp Q=15(8-s).             (3)

For ``s<8``, nonnegativity forces ``Q>=7-s``.  The total excess above this
minimum is only ``8+s<16``, so one of the sixteen directions has exactly
``Q=7-s`` and mean ``16-2s in {16,12,8,4}``.  A nonzero ``b`` costs at least
32 in phase zero.  At ``b=0`` the slack is ``A=2C`` with ``C`` a nonzero
nonnegative integral quadratic, and Proposition 15.688 gives
``4p E C>=28``.  Both alternatives contradict that mean.  Thus ``s=8``.

It follows that ``E_opp=0`` and

    I=124-16P,       0<=P<=7.                            (4)

The case ``P=0`` has ``I=124`` and only one finite edge, whereas
``D=U triangle partial(F)`` would give ``I<=|D|+2|F|=34``.  Hence
``1<=P<=7``.  Every selected finite edge now has hard sign, every
opposite-type direction has ``P_d=0``, and all sixteen opposite-type means
equal 32.

The phase-zero floor table permits only ``b=0,2,30`` at mean 32.  A phase-zero
``b=2`` equality has ``A=(x_i-x_j)^2`` and coefficient congruence
``15 | I-4``.  A phase-zero ``b=30`` equality has ``A=x_j`` and congruence
``15 | I-5``.  But (4) gives ``I=4-P (mod 15)`` with ``1<=P<=7``, so neither
congruence can hold.  Thus all sixteen opposite directions have ``b=0``.
The fifteen hard baselines contribute
only 30 to the global identity ``sum_d b_d=72``, so the remaining hard
direction would need ``b=42``.  This is impossible for 31 fibres (indeed
``b`` is even and at most 30).

Hence the ``p=31,R=10`` endpoint is excluded.  As an optional strengthening,
the entrywise norm in (1) forces ``rho`` even and reduces (4) further to
``(I,P)=(28,6),(60,4),(92,2)`` with explicit infinity-star fibre profiles;
the close does not need that refinement.  This is a symbolic coefficient
argument, not a configuration census.  Larger slack at ``p=31``, the endpoint
from ``p=37`` onward, residual (ii), Type I, and the limit remain open.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15727 import endpoint_block_row
from e1_gmin_m4_prop15728 import p31_even_floor_table


ROOT = Path(__file__).resolve().parents[1]
P = 31
R = 10
EDGE_COUNT = 4 * P + 1
TYPE_SIZE = (P + 1) // 2
TYPE_BUDGET = (P + 1) ** 2 // 2
GLOBAL_B_SUM = P + 1 + 4 * R


def p31_baseline_coefficient_rules() -> dict[str, object]:
    """Record the four exact baseline targets and congruences used."""
    return {
        "middle_slice_coordinates": "z_s=2x_s-1 and sum_s z_s=1",
        "hard_phase_one_b2": {
            "A": "(1-x_i-x_j)^2",
            "eps_S_H": "4+z_i z_j",
            "congruence": "15 divides I+P_d-4",
            "offset": 4,
        },
        "hard_phase_one_b30": {
            "A": "1-x_j after choosing the omitted fibre j",
            "eps_S_H": "4-z_j",
            "congruence": "15 divides I+P_d-3",
            "offset": 3,
        },
        "opposite_phase_zero_b2": {
            "A": "(x_i-x_j)^2",
            "eps_S_H": "4-z_i z_j",
            "congruence": "15 divides I+P_d-4",
            "offset": 4,
        },
        "opposite_phase_zero_b30": {
            "A": "x_j after choosing the omitted fibre j",
            "eps_S_H": "4+z_j",
            "congruence": "15 divides I+P_d-5",
            "offset": 5,
        },
        "proved": True,
    }


def hard_mean30_b2_upgrade() -> dict[str, object]:
    """Upgrade the fourteen-direction floor to all fifteen baselines."""
    rules = p31_baseline_coefficient_rules()
    return {
        "hard_type_sign": "eps_d=c_H",
        "mean_30_direction_count": 15,
        "previous_b2_floor": 14,
        "equal_means_force_equal_parallel_counts": True,
        "mean_formula": "a_d=I+32P_d-eps_d T-93",
        "b2_congruence_offset": rules["hard_phase_one_b2"]["offset"],
        "b30_congruence_offset": rules["hard_phase_one_b30"]["offset"],
        "coexistence_possible": False,
        "coexistence_contradiction": (
            "the same I+P cannot be congruent to both 4 and 3 modulo 15"
        ),
        "all_mean_30_directions_have_b2": True,
        "proved": True,
    }


def short_type_collapse_rows() -> list[dict[str, object]]:
    """Use opposite-type mass to force ``s=rho+P=8`` directly."""
    lift_floor = int(sharp_integral_quadratic_lift_floor(P)["sharp_scaled_floor"])
    rows: list[dict[str, object]] = []
    for s in (0, 2, 4, 6, 8):
        opposite_edge_count = 15 * (8 - s)
        if s < 8:
            minimum_parallel_count = 7 - s
            minimum_mean = 16 - 2 * s
            total_parallel_excess = opposite_edge_count - 16 * minimum_parallel_count
            base_direction_forced = total_parallel_excess < 16
            excluded = bool(
                base_direction_forced
                and 0 < minimum_mean < 32
                and minimum_mean < lift_floor
            )
            reason = (
                "one opposite direction has Q=7-s and mean 16-2s; nonzero "
                "b costs at least 32, while b=0 violates the lift floor 28"
            )
        else:
            minimum_parallel_count = 0
            minimum_mean = 32
            total_parallel_excess = 0
            base_direction_forced = True
            excluded = False
            reason = "all opposite-sign finite edges vanish"
        rows.append(
            {
                "s": s,
                "opposite_sign_finite_edge_count": opposite_edge_count,
                "opposite_parallel_count_sum": opposite_edge_count,
                "minimum_parallel_count_per_opposite_direction": minimum_parallel_count,
                "minimum_direction_mean": minimum_mean,
                "total_parallel_excess_above_minimum": total_parallel_excess,
                "a_minimum_direction_is_forced": base_direction_forced,
                "nonzero_integral_lift_floor": lift_floor,
                "excluded": excluded,
                "reason": reason,
            }
        )
    return rows


def short_s8_parallel_rows() -> list[dict[str, object]]:
    """List the ``s=8`` rows and apply boundary support plus congruences."""
    floors = p31_even_floor_table()["phase_zero_floors"]
    compatible_b = [int(b) for b, floor in floors.items() if int(floor) <= 32]
    rows: list[dict[str, object]] = []
    for parallel_count in range(8):
        infinity_degree = 124 - 16 * parallel_count
        finite_edges = EDGE_COUNT - infinity_degree
        boundary_support_upper = 32 + 2 * finite_edges
        boundary_support_excluded = infinity_degree > boundary_support_upper
        b2_congruence = (infinity_degree - 4) % 15 == 0
        b30_congruence = (infinity_degree - 5) % 15 == 0
        all_opposite_b0 = bool(
            not boundary_support_excluded
            and not b2_congruence
            and not b30_congruence
        )
        required_high_b = GLOBAL_B_SUM - 15 * 2
        contradiction = all_opposite_b0 and required_high_b > P - 1
        rows.append(
            {
                "P": parallel_count,
                "I": infinity_degree,
                "finite_edge_count": finite_edges,
                "boundary_support_upper_for_I": boundary_support_upper,
                "boundary_support_excluded": boundary_support_excluded,
                "opposite_direction_parallel_count": 0,
                "opposite_direction_mean": 32,
                "phase_zero_floor_compatible_b": compatible_b,
                "b2_congruence_I_minus_4": b2_congruence,
                "b30_congruence_I_minus_5": b30_congruence,
                "all_opposite_directions_have_b0": all_opposite_b0,
                "required_high_hard_direction_b": required_high_b,
                "maximum_possible_even_b": P - 1,
                "contradiction": contradiction,
                "excluded": boundary_support_excluded or contradiction,
            }
        )
    return rows


def hard_b2_coefficient_row(infinity_degree: int, parallel_count: int) -> dict[str, object]:
    """Return the exact coefficient/norm row for a hard b=2 baseline."""
    if (
        not isinstance(infinity_degree, int)
        or isinstance(infinity_degree, bool)
        or not isinstance(parallel_count, int)
        or isinstance(parallel_count, bool)
        or infinity_degree < 0
        or parallel_count < 0
        or infinity_degree + parallel_count > EDGE_COUNT
    ):
        raise ValueError("need nonnegative integer I,P with I+P<=125")
    numerator = infinity_degree + parallel_count - 4
    integral_g = numerator >= 0 and numerator % 15 == 0
    g = numerator // 15 if integral_g else None
    capacity = EDGE_COUNT - infinity_degree - parallel_count
    return {
        "I": infinity_degree,
        "P": parallel_count,
        "g": g,
        "g_integral_nonnegative": integral_g,
        "coefficient_identity": (
            "K_st=eps(g-n_s-n_t+1_{st=B_d}), with sum_s n_s=I"
        ),
        "entrywise_norm_bound": (
            "sum_{s<t}|g-n_s-n_t+1_{st=B_d}| <= E-P"
        ),
        "transverse_edge_capacity": capacity,
        "capacity_formula_121_minus_15g": (
            121 - 15 * g if g is not None else None
        ),
        "proved": bool(integral_g and capacity == 121 - 15 * int(g)),
    }


def odd_g_norm_obstruction() -> dict[str, object]:
    """Prove that the common coefficient parameter g cannot be odd."""
    pair_count = P * (P - 1) // 2
    maximum_zero_pairs = (P * P) // 4
    base_norm_floor = pair_count - maximum_zero_pairs
    corrected_norm_floor = base_norm_floor - 1
    odd_g_capacity_ceiling = 121 - 15
    return {
        "a_s": "2n_s-g, a nonzero odd integer when g is odd",
        "zero_cell_condition": "a_s=-a_t",
        "pair_count": pair_count,
        "maximum_zero_pairs": maximum_zero_pairs,
        "maximum_zero_pair_reason": (
            "positive_count times negative_count is at most floor(31^2/4)"
        ),
        "base_entrywise_norm_floor": base_norm_floor,
        "one_special_cell_correction": -1,
        "corrected_entrywise_norm_floor": corrected_norm_floor,
        "odd_g_is_at_least": 1,
        "transverse_capacity_ceiling": odd_g_capacity_ceiling,
        "contradiction": corrected_norm_floor > odd_g_capacity_ceiling,
        "g_must_be_even": True,
        "proved": (
            pair_count == 465
            and maximum_zero_pairs == 240
            and corrected_norm_floor == 224
            and corrected_norm_floor > odd_g_capacity_ceiling
        ),
    }


def pre_phase_zero_candidate_rows() -> list[dict[str, int]]:
    """Enumerate the elementary even-g rows after the aggregate norm bound."""
    rows: list[dict[str, int]] = []
    for infinity_degree in range(0, 95, 2):
        for parallel_count in range(0, EDGE_COUNT - infinity_degree + 1):
            numerator = infinity_degree + parallel_count - 4
            if numerator < 0 or numerator % 30:
                continue
            r = numerator // 30
            if r > 4:
                continue
            signed_aggregate_parameter = infinity_degree - P * r
            transverse_capacity = 121 - 30 * r
            signed_aggregate = 1 - 30 * signed_aggregate_parameter
            if transverse_capacity < 0 or abs(signed_aggregate) > transverse_capacity:
                continue
            rows.append(
                {
                    "r": r,
                    "g": 2 * r,
                    "I": infinity_degree,
                    "P": parallel_count,
                    "B0": signed_aggregate_parameter,
                    "signed_cell_aggregate": signed_aggregate,
                    "transverse_capacity": transverse_capacity,
                }
            )
    return rows


def phase_zero_filter_rows() -> list[dict[str, object]]:
    """Apply the opposite-type residue and the 15.688 lift floor."""
    lift_floor = int(sharp_integral_quadratic_lift_floor(P)["sharp_scaled_floor"])
    rows: list[dict[str, object]] = []
    for candidate in pre_phase_zero_candidate_rows():
        infinity_degree = int(candidate["I"])
        u0 = (infinity_degree + 4) % 16
        quotient_sum = 16 - u0
        if u0 == 0:
            excluded = False
            reason = "phase-zero residue is zero"
            forced_small_mean = None
        else:
            forced_small_mean = 2 * u0
            excluded = bool(quotient_sum < 16 and forced_small_mean < lift_floor)
            reason = (
                "some quotient k_d is zero; nonzero b costs at least 32, "
                "while b=0 gives A=2C and 4p E[C]>=28"
            )
        rows.append(
            {
                **candidate,
                "phase_zero_u0": u0,
                "phase_zero_quotient_sum": quotient_sum,
                "forced_zero_quotient_direction_mean": forced_small_mean,
                "nonzero_integral_lift_floor": lift_floor,
                "excluded": excluded,
                "reason": reason,
            }
        )
    return rows


def surviving_global_rows() -> list[dict[str, object]]:
    """Return the three rows and their forced sign/parallel consequences."""
    rows: list[dict[str, object]] = []
    for candidate in phase_zero_filter_rows():
        if candidate["excluded"]:
            continue
        r = int(candidate["r"])
        infinity_degree = int(candidate["I"])
        parallel_count = int(candidate["P"])
        finite_edges = EDGE_COUNT - infinity_degree
        hard_parallel_total = 15 * parallel_count + (parallel_count + 1)
        checks = {
            "row_formula_I": infinity_degree == 32 * r - 4,
            "row_formula_P": parallel_count == 8 - 2 * r,
            "hard_high_parallel_increment": hard_parallel_total == 16 * parallel_count + 1,
            "hard_parallel_counts_exhaust_finite_edges": hard_parallel_total == finite_edges,
            "opposite_mean": infinity_degree + finite_edges - 3 * P == 32,
        }
        rows.append(
            {
                "r": r,
                "I": infinity_degree,
                "baseline_hard_parallel_count_P": parallel_count,
                "high_hard_parallel_count": parallel_count + 1,
                "finite_edge_count": finite_edges,
                "hard_parallel_total": hard_parallel_total,
                "every_finite_selected_edge_has_hard_sign": True,
                "every_opposite_direction_parallel_count": 0,
                "every_opposite_direction_mean": 32,
                "checks": checks,
                "proved": all(checks.values()),
            }
        )
    return rows


def infinity_star_profile_options(r: int) -> dict[str, object]:
    """List the exact profile options implied by the unsigned cell law."""
    if not isinstance(r, int) or isinstance(r, bool) or r not in (1, 2, 3):
        raise ValueError("the surviving rows have r in {1,2,3}")
    deficit_sum = 4 - r
    histograms: list[dict[int, int]] = []

    def visit(deficit: int, remaining: int, counts: list[int]) -> None:
        if deficit > r:
            if remaining == 0:
                histogram = {r: P - sum(counts)}
                for value, multiplicity in enumerate(counts, start=1):
                    if multiplicity:
                        n_value = r - value
                        histogram[n_value] = histogram.get(n_value, 0) + multiplicity
                histograms.append(dict(sorted(histogram.items())))
            return
        for multiplicity in range(remaining // deficit + 1):
            visit(deficit + 1, remaining - deficit * multiplicity, counts + [multiplicity])

    visit(1, deficit_sum, [])
    histograms.sort(key=lambda row: tuple(sorted(row.items())))
    return {
        "r": r,
        "deficit_variables": "d_s=r-n_s>=0",
        "deficit_sum": deficit_sum,
        "negative_deficit_impossible": True,
        "negative_deficit_reason": (
            "one d_s=-h forces the other 30 deficits to sum at least 30h-1, "
            "so the total is at least 29h-1>=28>4-r"
        ),
        "infinity_endpoint_count_histograms": histograms,
        "transverse_cell_count": "N_st=d_s+d_t+1_{st=B_d}",
        "proved": bool(histograms),
    }


def opposite_type_final_contradiction_row(
    infinity_degree: int,
) -> dict[str, object]:
    """Eliminate one of the three rows using the phase-zero floor table."""
    if infinity_degree not in (28, 60, 92):
        raise ValueError("need one of the three surviving infinity degrees")
    floors = p31_even_floor_table()["phase_zero_floors"]
    compatible_b = [int(b) for b, floor in floors.items() if int(floor) <= 32]
    b2_congruence = (infinity_degree - 4) % 15 == 0
    b30_congruence = (infinity_degree - 5) % 15 == 0
    all_opposite_b_zero = not b2_congruence and not b30_congruence
    hard_baseline_b_sum = 15 * 2
    required_high_b = GLOBAL_B_SUM - hard_baseline_b_sum
    maximum_even_b = P - 1
    contradiction = all_opposite_b_zero and required_high_b > maximum_even_b
    return {
        "I": infinity_degree,
        "opposite_direction_count": TYPE_SIZE,
        "opposite_direction_mean": 32,
        "phase_zero_floor_compatible_b": compatible_b,
        "phase_zero_b2_congruence": "15 divides I-4",
        "phase_zero_b2_congruence_holds": b2_congruence,
        "phase_zero_b30_congruence": "15 divides I-5",
        "phase_zero_b30_congruence_holds": b30_congruence,
        "all_opposite_directions_have_b0": all_opposite_b_zero,
        "global_sum_b": GLOBAL_B_SUM,
        "hard_baseline_b_sum": hard_baseline_b_sum,
        "required_high_hard_direction_b": required_high_b,
        "maximum_possible_even_b": maximum_even_b,
        "contradiction": contradiction,
        "proved": compatible_b == [0, 2, 30] and contradiction,
    }


def p31_block_direction_upgrade_row(four_secants: int) -> dict[str, object]:
    """Upgrade 15.728's nonrich hard-direction floor from 4+y to 5+y."""
    block = endpoint_block_row(P, four_secants)
    y = int(block["four_secants_y"])
    rich_lines = int(block["trisecants_x"]) + y
    nonrich_floor = 15 - rich_lines
    return {
        "four_secants_y": y,
        "trisecants_x": block["trisecants_x"],
        "rich_line_count": rich_lines,
        "hard_b2_direction_count": 15,
        "nonrich_hard_b2_direction_floor": nonrich_floor,
        "nonrich_floor_formula": "5+y",
        "proved": bool(block["proved"] and nonrich_floor == 5 + y),
    }


def proposition_15733() -> dict[str, object]:
    """Package the symbolic p=31 endpoint exclusion."""
    rules = p31_baseline_coefficient_rules()
    upgrade = hard_mean30_b2_upgrade()
    short_collapse = short_type_collapse_rows()
    short_final = short_s8_parallel_rows()
    odd = odd_g_norm_obstruction()
    candidates = pre_phase_zero_candidate_rows()
    filtered = phase_zero_filter_rows()
    survivors = surviving_global_rows()
    profiles = [infinity_star_profile_options(r) for r in (1, 2, 3)]
    final_rows = [opposite_type_final_contradiction_row(I) for I in (28, 60, 92)]
    blocks = [p31_block_direction_upgrade_row(y) for y in range(6)]
    proved = bool(
        rules["proved"]
        and upgrade["proved"]
        and [row["s"] for row in short_collapse if not row["excluded"]] == [8]
        and all(bool(row["excluded"]) for row in short_final)
        and odd["proved"]
        and len(candidates) == 12
        and sum(not bool(row["excluded"]) for row in filtered) == 3
        and all(row["proved"] for row in survivors)
        and all(row["proved"] for row in profiles)
        and all(row["proved"] for row in final_rows)
        and all(row["proved"] for row in blocks)
    )
    return {
        "prop": "15.733",
        "title": "Simultaneous-baseline coefficient exclusion at p=31",
        "result_status": "proved p=31 endpoint exclusion",
        "hypotheses": (
            "the p=31,R=10 residual endpoint normal form and Paley-hard "
            "type of Proposition 15.728"
        ),
        "baseline_coefficient_rules": rules,
        "hard_mean30_upgrade": upgrade,
        "canonical_short_type_collapse": short_collapse,
        "canonical_short_final_rows": short_final,
        "odd_g_obstruction": odd,
        "pre_phase_zero_candidates": candidates,
        "phase_zero_filter": filtered,
        "optional_three_row_strengthening": survivors,
        "optional_infinity_star_profiles": profiles,
        "optional_three_row_final_checks": final_rows,
        "block_direction_upgrades": blocks,
        "finite_configuration_search_used": False,
        "p31_R10_endpoint_excluded": True,
        "p31_first_possible_positive_slack": 11,
        "first_unexcluded_endpoint_prime": 37,
        "endpoint_all_primes_closed": False,
        "p_plus_one_shell_closed": False,
        "non_walsh_residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "next_gate": (
            "historical at this stage: Proposition 15.734 subsequently closes "
            "every k=4p boundary for p>=13; the live residual front is p=11 "
            "sharp equality or even k>4p"
        ),
        "superseded_as_live_endpoint_gate_by": "Proposition 15.734",
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic coefficient certificate."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15733.json"
    payload = json.dumps(proposition_15733(), indent=2, sort_keys=True) + "\n"
    output.write_text(payload)
    return output


def main() -> None:
    result = proposition_15733()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.733 p=31 endpoint audit failed")
    path = write_evidence()
    print("Prop 15.733 p=31,R=10 endpoint: excluded")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
