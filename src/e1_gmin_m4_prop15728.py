#!/usr/bin/env python3
r"""Prop. 15.728 -- the Paley-hard normal form at the p=31 endpoint.

Assume the residual-(ii) affine hypotheses at ``p=31`` and ``|H|=4p+1``.
In an outside chart let the all-finite odd-degree boundary ``D`` have size
``p+1=32`` and endpoint pair slack ``R=10``.  Proposition 15.727 then gives
the disjoint trisecant/4-secant normal form.

For a direction ``d``, let ``b_d`` be the number of affine fibres meeting
``D`` oddly.  Since ``|D|`` is even, every ``b_d`` is even, and the exact
pair identity gives

    sum_d b_d = (p+1)+4R = 72.                         (1)

The directional parity sign from Proposition 15.632 simplifies here to

    (-1)^A = -eps_d*c_H,

because ``(|H|-3)/2=61`` is odd, infinity is not in ``D``, and ``b_d`` is
even.  Thus the sixteen directions with ``eps_d=c_H`` all have phase one;
call this the Paley-hard type.  Each of the two quadratic types has exact
scaled-mean budget 512.

At ``p=31`` the phase-one floors are 30 for ``b=2,30`` and 62 for every
other even ``b``.  Same-type means have one common residue modulo 32.  Write

    a_d = 2u+32k_d,       0<=u<16,       sum_d k_d=16-u. (2)

For ``1<=u<=14`` every direction needs ``k_d>=1``, already more than the
right side of (2).  At ``u=0``, only ``b=2,30`` could use ``k_d=1``.  In
either cell the pointwise parity baseline is a Boolean quadratic ``q_0`` of
scaled mean 30.  If ``A`` had scaled mean 32, then ``C=(A-q_0)/2`` would be
a nonzero nonnegative integral quadratic with ``4p E[C]=2``.  Proposition
15.688 instead gives ``4p E[C]>=p-3=28``.  Hence ``u=15``.  Exactly fifteen hard-type
directions have ``b in {2,30}`` and mean 30, and the remaining direction has
mean 62 (it may itself have ``b=2`` or ``30`` after one 32-unit elevation).

Equation (1) permits at most one hard-type ``b=30`` direction.  Therefore
at least fourteen directions of one Paley type have ``b=2``.  If ``y`` is
the number of 4-secants in Proposition 15.727, there are only ``10-y`` rich
lines, so at least ``4+y`` of these directions contain no rich line.  Each
such direction has occupancy profile

    14 empty fibres, 2 singleton fibres, 15 two-point fibres.       (3)

This is a proved necessary normal form, not an exclusion of the endpoint.
No arc classification and no finite configuration search is used.  The
``p+1`` shell, residual (ii), Type I, and the limit remain open.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15727 import endpoint_block_row


ROOT = Path(__file__).resolve().parents[1]
P = 31
PERIOD = P + 1
M = PERIOD // 2
R = 10
EDGE_COUNT = 4 * P + 1
TYPE_BUDGET = PERIOD * PERIOD // 2
EVEN_B = tuple(range(0, P + 1, 2))
SPECIAL_PHASE_ONE_B = (2, P - 1)


def p31_endpoint_odd_fibre_ledger() -> dict[str, object]:
    """Record parity and the exact global odd-fibre sum at ``R=10``."""
    total = PERIOD + 4 * R
    return {
        "p": P,
        "boundary_size": PERIOD,
        "endpoint_pair_slack_R": R,
        "direction_count": PERIOD,
        "possible_b_values": list(EVEN_B),
        "every_b_even": True,
        "pair_slack_identity": "R=(sum_d b_d-(p+1))/4",
        "sum_b": total,
        "proved": total == 72,
    }


def p31_type_phase_ledger() -> dict[str, object]:
    """Identify the phase-one quadratic type in the all-finite chart."""
    exponent = (EDGE_COUNT - 3) // 2
    return {
        "p": P,
        "edge_count": EDGE_COUNT,
        "boundary_contains_infinity": False,
        "all_b_even": True,
        "parity_exponent": exponent,
        "parity_exponent_is_odd": exponent % 2 == 1,
        "directional_parity_sign": "-eps_d*c_H",
        "phase_one_type": "eps_d=c_H",
        "phase_zero_type": "eps_d=-c_H",
        "directions_per_type": M,
        "budget_per_type": TYPE_BUDGET,
        "both_type_budgets_are_exact": True,
        "proved": exponent == 61 and M == 16 and TYPE_BUDGET == 512,
    }


def p31_even_floor_table() -> dict[str, object]:
    """Return both exact Paley-type floor tables for even ``b``."""
    floors = {
        phase: {b: full_symbolic_floor(P, b, phase) for b in EVEN_B}
        for phase in (0, 1)
    }
    expected_zero = {
        0: 0,
        2: 32,
        4: 56,
        **{b: 62 for b in range(6, 27, 2)},
        28: 56,
        30: 32,
    }
    expected_one = {
        b: (30 if b in SPECIAL_PHASE_ONE_B else 62) for b in EVEN_B
    }
    if floors[0] != expected_zero or floors[1] != expected_one:
        raise ArithmeticError("the p=31 even floor table changed")
    return {
        "p": P,
        "even_b_values": list(EVEN_B),
        "phase_zero_floors": floors[0],
        "phase_one_floors": floors[1],
        "phase_one_special_b": list(SPECIAL_PHASE_ONE_B),
        "phase_one_special_floor": 30,
        "phase_one_other_floor": 62,
        "proved": True,
    }


def p31_endpoint_floor_plus_two_obstruction(b: int) -> dict[str, object]:
    """Exclude scaled mean 32 in either phase-one endpoint cell.

    Pointwise parity makes the Boolean baseline ``q_0`` a lower bound for
    every nonnegative integral slack quadratic ``A``.  Hence
    ``C=(A-q_0)/2`` is another nonnegative integral quadratic.  A two-unit
    scaled excess would give ``4p E[C]=2``, contradicting Proposition
    15.688's sharp nonzero-lift floor ``p-3=28``.
    """
    if b not in SPECIAL_PHASE_ONE_B:
        raise ValueError("the endpoint obstruction is for b=2 or b=30")
    lift = sharp_integral_quadratic_lift_floor(P)
    baseline = (
        "q_0=(1-x_i-x_j)^2"
        if b == 2
        else "q_0=1-x_i after complementing the 30-set"
    )
    baseline_scaled_mean = full_symbolic_floor(P, b, 1)
    proposed_scaled_mean = baseline_scaled_mean + 2
    induced_lift_scaled_mean = proposed_scaled_mean - baseline_scaled_mean
    lift_floor = int(lift["sharp_scaled_floor"])
    forbidden = induced_lift_scaled_mean < lift_floor
    return {
        "p": P,
        "b": b,
        "phase": 1,
        "pointwise_boolean_baseline": baseline,
        "parity_forces_A_at_least_q0": True,
        "baseline_scaled_mean_2p_E_q0": baseline_scaled_mean,
        "proposed_scaled_mean_2p_E_A": proposed_scaled_mean,
        "lift": "C=(A-q_0)/2",
        "lift_is_nonnegative_integral_quadratic": True,
        "induced_scaled_lift_mass_4p_E_C": induced_lift_scaled_mean,
        "prop_15_688_scaled_lift_floor": lift_floor,
        "floor_plus_two_forbidden": forbidden,
        "proved": forbidden and lift_floor == P - 3,
    }


def p31_phase_one_residue_ledger() -> dict[str, object]:
    """Classify the common residue in the phase-one quadratic type."""
    floor_two = p31_endpoint_floor_plus_two_obstruction(2)
    floor_thirty = p31_endpoint_floor_plus_two_obstruction(30)
    floor_plus_two_forbidden = bool(
        floor_two["floor_plus_two_forbidden"]
        and floor_thirty["floor_plus_two_forbidden"]
    )
    if not floor_plus_two_forbidden:
        raise ArithmeticError("a phase-one endpoint floor-plus-two cell reopened")

    rows = []
    for u in range(M):
        quotient_sum = M - u
        if u == 0:
            minimum_quotient_per_direction = 2
            minimum_total = 2 * M
            feasible = False
            reason = (
                "k=1 is possible only at b=2,30 and is the forbidden "
                "floor-plus-two lift"
            )
        elif u < M - 1:
            minimum_quotient_per_direction = 1
            minimum_total = M
            feasible = False
            reason = "every direction needs k>=1 but sum k=16-u<16"
        else:
            minimum_quotient_per_direction = 0
            minimum_total = 0
            feasible = True
            reason = (
                "sum k=1: fifteen special directions have k=0 and the "
                "remaining direction has k=1"
            )
        rows.append(
            {
                "u": u,
                "common_residue_2u": 2 * u,
                "quotient_sum": quotient_sum,
                "minimum_quotient_per_direction": minimum_quotient_per_direction,
                "minimum_quotient_sum": minimum_total,
                "feasible": feasible,
                "reason": reason,
            }
        )

    feasible_u = [int(row["u"]) for row in rows if row["feasible"]]
    if feasible_u != [15]:
        raise ArithmeticError("the p=31 phase-one residue changed")
    return {
        "p": P,
        "period": PERIOD,
        "directions": M,
        "type_budget": TYPE_BUDGET,
        "same_type_mean_form": "a_d=2u+32*k_d",
        "quotient_identity": "sum_d k_d=16-u",
        "floor_plus_two_cells": {"b=2": floor_two, "b=30": floor_thirty},
        "residue_rows": rows,
        "unique_u": 15,
        "common_residue": 30,
        "mean_multiset": {30: 15, 62: 1},
        "baseline_direction_count": 15,
        "baseline_b_values": list(SPECIAL_PHASE_ONE_B),
        "high_direction_mean": 62,
        "high_direction_b": "any even b; b=2 or 30 uses one 32-unit elevation",
        "proved": True,
    }


def p31_phase_zero_residue_parity() -> dict[str, object]:
    """Use the opposite type and the even infinity degree to halve its rows.

    The exact directional mean is

    ``a_d=I+32*P_d-eps_d*T-93``.

    Adding the residues in the two types cancels ``T``.  The phase-one
    residue is 30, so if the phase-zero residue is ``2*u0``, then
    ``u0=I+4 (mod 16)``.  Infinity is outside the all-finite boundary and
    therefore has even degree ``I``.
    """
    candidates = list(range(0, M, 2))
    return {
        "exact_directional_mean": "a_d=I+32*P_d-eps_d*T-93",
        "phase_one_residue": 30,
        "phase_zero_residue": "2*u0",
        "type_residue_sum": "30+2*u0=2*I+6 (mod 32)",
        "halved_congruence": "u0=I+4 (mod 16)",
        "infinity_not_in_boundary": True,
        "infinity_degree_I_even": True,
        "possible_phase_zero_u0": candidates,
        "possible_phase_zero_residues": [2 * u for u in candidates],
        "proved": candidates == [0, 2, 4, 6, 8, 10, 12, 14],
    }


def p31_hard_type_b_split_cases() -> dict[str, object]:
    """List the exact global ``sum b=72`` splits after hard-type rigidity.

    These are necessary arithmetic cases, not claimed realizable profiles.
    ``high_b`` is the odd-fibre count in the unique mean-62 direction.
    """
    cases: list[dict[str, object]] = []

    # No b=30 direction: all fifteen mean-30 directions have b=2.
    for high_b in range(0, 29, 2):
        hard_sum = 30 + high_b
        cases.append(
            {
                "b30_location": "none",
                "mean_30_profile": {2: 15, 30: 0},
                "high_b": high_b,
                "hard_type_sum_b": hard_sum,
                "other_type_sum_b": 72 - hard_sum,
            }
        )

    # The sole b=30 direction is one of the fifteen mean-30 baselines.
    for high_b in range(0, 15, 2):
        hard_sum = 30 + 14 * 2 + high_b
        cases.append(
            {
                "b30_location": "mean_30_baseline",
                "mean_30_profile": {2: 14, 30: 1},
                "high_b": high_b,
                "hard_type_sum_b": hard_sum,
                "other_type_sum_b": 72 - hard_sum,
            }
        )

    # Or it is the unique mean-62 direction.
    cases.append(
        {
            "b30_location": "mean_62_high_direction",
            "mean_30_profile": {2: 15, 30: 0},
            "high_b": 30,
            "hard_type_sum_b": 60,
            "other_type_sum_b": 12,
        }
    )
    if any(int(row["other_type_sum_b"]) < 0 for row in cases):
        raise ArithmeticError("a negative other-type b sum entered the ledger")
    return {
        "global_sum_b": 72,
        "necessary_case_count": len(cases),
        "cases": cases,
        "at_most_one_hard_type_b30_direction": True,
        "reason": (
            "two b=30 directions plus the other thirteen required special "
            "directions at b=2 would already contribute 86>72"
        ),
        "minimum_hard_type_b2_directions": 14,
        "proved": len(cases) == 24,
    }


def p31_b2_direction_occupancy(
    trisecants_in_direction: int, four_secants_in_direction: int
) -> dict[str, object]:
    """Recover every fibre count in a direction with ``b_d=2``."""
    r3 = trisecants_in_direction
    r4 = four_secants_in_direction
    if r3 < 0 or r4 < 0 or r3 > 2 or r3 + 2 * r4 > 15:
        raise ValueError("invalid rich-line counts for a b=2 direction")
    counts = {
        0: 14 + r3 + r4,
        1: 2 - r3,
        2: 15 - r3 - 2 * r4,
        3: r3,
        4: r4,
    }
    point_count = sum(occupancy * count for occupancy, count in counts.items())
    fibre_count = sum(counts.values())
    odd_fibres = counts[1] + counts[3]
    slack = r3 + 2 * r4
    pair_count = sum(
        occupancy * (occupancy - 1) // 2 * count
        for occupancy, count in counts.items()
    )
    proved = bool(
        point_count == PERIOD
        and fibre_count == P
        and odd_fibres == 2
        and pair_count == 15 + 2 * slack
        and all(count >= 0 for count in counts.values())
    )
    return {
        "p": P,
        "b": 2,
        "trisecants_in_direction": r3,
        "four_secants_in_direction": r4,
        "direction_slack": slack,
        "occupancy_counts": counts,
        "point_count": point_count,
        "fibre_count": fibre_count,
        "pair_count": pair_count,
        "nonrich_profile": r3 == 0 and r4 == 0,
        "proved": proved,
    }


def p31_endpoint_block_paley_rows() -> list[dict[str, object]]:
    """Add the hard-type pairing-direction lower bound to each block row."""
    rows = []
    for y in range(R // 2 + 1):
        block = endpoint_block_row(P, y)
        rich_lines = int(block["trisecants_x"]) + int(block["four_secants_y"])
        nonrich_b2 = 14 - rich_lines
        rows.append(
            {
                "four_secants_y": y,
                "trisecants_x": int(block["trisecants_x"]),
                "rich_line_count": rich_lines,
                "maximum_rich_directions": rich_lines,
                "minimum_hard_type_b2_directions": 14,
                "minimum_nonrich_hard_type_b2_directions": nonrich_b2,
                "minimum_nonrich_formula": "4+y",
                "nonrich_b2_occupancy": p31_b2_direction_occupancy(0, 0)[
                    "occupancy_counts"
                ],
                "endpoint_block": block,
                "proved": nonrich_b2 == 4 + y and bool(block["proved"]),
            }
        )
    return rows


def proposition_15728() -> dict[str, object]:
    """Package the necessary p=31 endpoint normal form honestly."""
    odd = p31_endpoint_odd_fibre_ledger()
    phases = p31_type_phase_ledger()
    floors = p31_even_floor_table()
    residue = p31_phase_one_residue_ledger()
    phase_zero = p31_phase_zero_residue_parity()
    split = p31_hard_type_b_split_cases()
    blocks = p31_endpoint_block_paley_rows()
    proved = bool(
        odd["proved"]
        and phases["proved"]
        and floors["proved"]
        and residue["proved"]
        and phase_zero["proved"]
        and split["proved"]
        and all(row["proved"] for row in blocks)
    )
    return {
        "prop": "15.728",
        "title": "Paley-hard profile rigidity at the p=31 endpoint",
        "result_status": "proved necessary normal form",
        "hypotheses": {
            "residual_affine_separator": True,
            "p": P,
            "edge_count": EDGE_COUNT,
            "all_finite_boundary_size": PERIOD,
            "outside_pair_slack": R,
            "prop_15_727_disjoint_block_normal_form": True,
        },
        "odd_fibre_ledger": odd,
        "type_phase_ledger": phases,
        "floor_table": floors,
        "phase_one_residue_normal_form": residue,
        "phase_zero_residue_parity": phase_zero,
        "hard_type_b_split_cases": split,
        "endpoint_block_paley_rows": blocks,
        "conclusion": {
            "one_paley_type": "eps_d=c_H",
            "mean_multiset": {30: 15, 62: 1},
            "minimum_b2_directions": 14,
            "possible_phase_zero_half_residues": [0, 2, 4, 6, 8, 10, 12, 14],
            "minimum_nonrich_b2_directions": "4+y when y is the 4-secant count",
            "nonrich_b2_profile": {0: 14, 1: 2, 2: 15, 3: 0, 4: 0},
        },
        "finite_configuration_search_used": False,
        "arc_classification_used": False,
        "p31_endpoint_excluded": False,
        "p_plus_one_shell_closed": False,
        "non_walsh_residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "next_gate": (
            "exclude a 32-point endpoint set having at least 4+y same-Paley "
            "nonrich directions of profile 14/2/15"
        ),
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic arithmetic certificate."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15728.json"
    output.write_text(json.dumps(proposition_15728(), indent=2, sort_keys=True) + "\n")
    return output


def main() -> None:
    result = proposition_15728()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.728 audit failed")
    path = write_evidence()
    print("Prop 15.728 p=31 Paley-hard endpoint normal form: proved")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
