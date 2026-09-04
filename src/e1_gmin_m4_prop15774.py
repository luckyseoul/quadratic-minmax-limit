#!/usr/bin/env python3
"""Prop. 15.774: sharp small-mass two-type capacity and two residual layers.

This is a bounded-size infinite-family theorem, NOT a global E1 or limit
proof. Scalar profiles at the first uncovered layer are not graph witnesses.
See evidence/NOTE_2026-09-04_SMALL_MASS_TWO_TYPE_BRIDGE.md.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15750 import type_I_multilevel_bad_case_closed_all_primes
from e1_gmin_m4_prop15764 import official_unit_entry_ledger, parity_bridge_ledger
from e1_gmin_m4_prop15773 import joint_layer_exclusion
from e1_gmin_m4_small_mass_spectrum import affine_parity_small_mass_spectrum
from io_atomic import write_json_atomic

ROOT = Path(__file__).resolve().parents[1]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _check(p: int, r: int = 3) -> None:
    if type(p) is not int or p < 29 or not is_prime(p):
        raise ValueError("need prime p>=29")
    if type(r) is not int or r not in (3, 4, 5):
        raise ValueError("need shell floor r in {3,4,5}")


def quotient_floor(p: int, u: int) -> int:
    """Necessary lower bound only; zero means need not be realizable."""
    _check(p)
    m = (p + 1) // 2
    if type(u) is not int or not 0 <= u < m:
        raise ValueError("residue outside [0,m)")
    if u in (0, m - 2, m - 1):
        return 0
    return 2 if u <= m - 7 else 1


def capacity_exclusion(p: int, r: int) -> dict[str, object]:
    """Both signed shell floors r, with H congruent to r modulo two."""
    _check(p, r)
    m = (p + 1) // 2
    spectrum = affine_parity_small_mass_spectrum(p)
    _require(spectrum["proved"] and spectrum["strict_upper_mass"] == 2 * p - 10
             and spectrum["union_allowed_masses"] == [0, p - 3, p - 1, p + 1],
             "affine parity union spectrum dependency failed")
    t_max = p + (r + 1) // 2 if p >= 37 else p - 6
    h_max = r * p + 2 * t_max
    isolated = p * p + 1 - 2 * h_max
    _require(isolated > 0, "isolated chart not justified")
    floors = [quotient_floor(p, u) for u in range(m)]
    # These lists replay the symbolic residue-interval argument, not primes,
    # graphs, local cells, or an equality catalog.
    residues = [u for u, k in enumerate(floors) if u + m * k <= t_max]
    if p >= 37:
        s = (r + 1) // 2 - 1
        expected = list(range(s + 1)) + list(range(m - 6, m))
        _require(residues == expected and 2 * s < r and m - 12 > r
                 and m - 6 > r and s - 1 < r,
                 "low/high residue separation failed")
    else:
        _require(residues == [0, m - 2, m - 1], "small-prime capacity changed")
    allowed = set(residues)
    collisions = [[u, (r - u) % m] for u in residues if (r - u) % m in allowed]
    _require(not collisions, "two types have a surviving residue collision")
    return {"p": p, "r": r, "m": m, "maximum_t": t_max,
            "maximum_H": h_max, "H_parity": r % 2,
            "guaranteed_isolated_vertices": isolated,
            "direction_mass": "a_eps=(p+1)P_eps-eps*T-r*p=2u_eps+(p+1)k",
            "type_quotient_sum": "t-u_eps", "residue_sum_mod_m": r,
            "quotient_floors": floors, "allowed_residues_at_maximum_t": residues,
            "collisions": collisions, "all_smaller_t_by_monotonicity": True,
            "H_below_rp_excluded_by_frame_average": True,
            "both_types_use_union_spectrum": True,
            "even_r_phase_one_assumed": False, "proved": True}


def first_scalar_survivor(p: int, r: int) -> dict[str, object]:
    """Exact necessary quota data at the first uncovered eventual layer."""
    _check(p, r)
    if p < 37:
        raise ValueError("sharp eventual threshold needs p>=37")
    m, t = (p + 1) // 2, p + 1 + (r + 1) // 2
    edges = r * p + 2 * t
    rows = []
    for u in (r // 2, (r + 1) // 2):
        extra = t - u - 2 * m
        _require(extra in (0, 1), "scalar surplus changed")
        quotients = [2 + extra] + [2] * (m - 1)
        parallel = [r + k for k in quotients]
        signed_T = r - 2 * u
        masses = [2 * u + (p + 1) * k for k in quotients]
        _require(sum(quotients) == t - u and min(masses) >= 2 * p - 10
                 and all((p + 1) * P - signed_T - r * p == a
                         for P, a in zip(parallel, masses)), "scalar identities failed")
        rows.append({"u": u, "quotients": quotients, "parallel_counts": parallel,
                     "signed_T": signed_T, "masses": masses})
    _require(sum(sum(x["parallel_counts"]) for x in rows) == edges
             and rows[0]["signed_T"] == -rows[1]["signed_T"], "scalar type glue failed")
    return {"p": p, "r": r, "t": t, "H": edges, "types": rows,
            "satisfies_scalar_relaxation": True, "graph_realization_claimed": False,
            "full_local_row_realization_claimed": False, "proved": True}


def _carry_branch(p: int, s: int, old: dict[str, object]) -> dict[str, object]:
    """Carry a LOCAL equality and recompute every H-dependent quantity."""
    q, m, edges = (p - 1) // 2, (p + 1) // 2, 5 * p + 2 * s
    u, low, P = old["u"], old["low_quotient"], old["coefficient_offset"]
    total = p + s - u
    upper = (edges - total + m * low) // m
    candidates = [j for j in range(upper + 1) if (j - P) % q == 0]
    _require(candidates == [P] and upper <= 9 < q, "carried normalization failed")
    hard = m * (P - low) + total
    hT = (p + 1) * P - 3 * p - 2 * u - (p + 1) * low
    opposite = edges - hard
    next_Q = old["forced_next_Q"]
    surplus = opposite - m * next_Q
    count = m - surplus
    _require(old["proved"] and old["opposite_local_exclusion"]["proved"]
             and 2 * hard == edges + hT
             and (p + 1) * next_Q + hT - 3 * p == old["forced_next_scaled_mean"]
             and count == old["forced_next_row_count_at_least"] - s > 0
             and surplus >= 0, "carried opposite-row contradiction failed")
    return {"u": u, "low_quotient": low, "P": P, "quotient_sum": total,
            "common_low_parallel_candidates": candidates, "low_parallel_upper_bound": upper,
            "hard_edges": hard, "hT": hT, "opposite_edges": opposite,
            "forced_Q": next_Q, "forced_mass": old["forced_next_scaled_mean"],
            "surplus": surplus, "forced_count": count,
            "carried_local_exclusion_proved": True, "proved": True}


def _uncatalogued_branch(p: int, s: int, u: int, low: int,
                         spectrum: dict[str, object]) -> dict[str, object]:
    """A minimum quotient determines common P without a new catalog."""
    m, edges = (p + 1) // 2, 5 * p + 2 * s
    total = p + s - u
    extra = total - m * low
    _require(0 <= extra < m and low in (1, 2), "minimum quotient is not forced")
    # u=q,low=1 gives mass p+9; u<s,low=2 gives p+7-2u.
    small = 8 if low == 1 else 6 - 2 * u
    mass = p + 1 + small
    count = 5 - s if low == 1 else 4 - s - u
    _require(spectrum["proved"] and 0 < small < p - 3
             and p + 1 < mass < spectrum["strict_upper_mass"] and count > 0,
             "uncatalogued local mass exclusion failed")
    upper = (edges - extra) // m
    _require(upper == 9, "uncatalogued parallel upper bound changed")
    rows = []
    for P in range(upper + 1):
        hard, next_Q = m * P + extra, 9 - P
        hT = (p + 1) * P - 3 * p - 2 * u - (p + 1) * low
        surplus = edges - hard - m * next_Q
        _require(2 * hard == edges + hT and surplus == m - count
                 and (p + 1) * next_Q + hT - 3 * p == mass,
                 "uncatalogued row budget failed")
        rows.append({"P": P, "hard_edges": hard, "hT": hT,
                     "forced_Q": next_Q, "surplus": surplus})
    return {"u": u, "low_quotient": low, "excess": extra,
            "low_row_count_at_least": m - extra, "small_forbidden_mass": small,
            "forced_mass": mass, "forced_count": count, "parallel_cases": rows,
            "new_equality_catalog_used": False, "P9_uses_R_not_negative_Q": True,
            "proved": True}


def residual_two_layer_exclusion(p: int, s: int) -> dict[str, object]:
    _check(p)
    if type(s) is not int or s not in (1, 2):
        raise ValueError("only the two proved layers s=1,2")
    q, m = (p - 1) // 2, (p + 1) // 2
    prior, spectrum = joint_layer_exclusion(p), affine_parity_small_mass_spectrum(p)
    _require(prior["proved"] and spectrum["proved"]
             and spectrum["union_allowed_masses"] == [0, p - 3, p - 1, p + 1],
             "carried layer or spectrum dependency failed")
    carried = {name: _carry_branch(p, s, old)
               for name, old in prior["carried_branch_exclusions"].items()}
    residues = []
    for u in range(q):
        # When k=1 is absent, sum k>=2m implies u<=s-1.
        residues.append({"u": u, "quotient_sum": p + s - u,
                         "zero_quotient_excluded": True,
                         "if_quotient_one_then_existing_local_classification": True,
                         "if_no_quotient_one_possible": u < s,
                         "quotient_one_count_lower_bound": max(0, u + 1 - s)})
    fresh = [_uncatalogued_branch(p, s, u, 2, spectrum) for u in range(s)]
    last = _uncatalogued_branch(p, s, q, 1, spectrum)
    _require(all(x["proved"] for x in carried.values()) and all(x["proved"] for x in fresh)
             and last["proved"] and p * p + 1 - 2 * (5 * p + 2 * s) > 0,
             "two-layer branch dependency failed")
    return {"p": p, "s": s, "t": q + s, "k": 5 * p + 2 * s - 1,
            "H": 5 * p + 2 * s, "hard_sign_relative_to_transported_c_H": (-1) ** (q + s),
            "residue_cases_below_q": residues, "carried_local_branches": carried,
            "fresh_quotient_two_branches": fresh, "q_nozero_branch": last,
            "q_zero_case_retains_arbitrary_quotient_heights": True,
            "all_boundary_sizes_excluded": True, "new_prime_or_local_census_used": False,
            "proved": True}


def minimal_four_gap_consequences(p: int, layers: list[dict[str, object]] | None = None) -> dict[str, object]:
    _check(p)
    odd, even = official_unit_entry_ledger(p, True), official_unit_entry_ledger(p, False)
    _require(odd["official_entry_proved"] and even["official_entry_proved"]
             and parity_bridge_ledger(5 * p)["proved"]
             and parity_bridge_ledger(4 * p)["proved"], "official bridge dependency failed")
    even_cap, no_bridge = capacity_exclusion(p, 4), capacity_exclusion(p, 5)
    if layers is None:
        layers = [residual_two_layer_exclusion(p, s) for s in (1, 2)]
    _require(even_cap["proved"] and no_bridge["proved"]
             and type_I_multilevel_bad_case_closed_all_primes()
             and len(layers) == 2
             and all(row["proved"] and row["p"] == p and row["s"] == s
                     for s, row in zip((1, 2), layers)), "minimal consequence dependency failed")
    return {"p": p, "odd_minimal_H_lower_bound": 5 * p + 6,
            "even_minimal_H_lower_bound": even_cap["maximum_H"] + 2,
            "odd_without_level_three_H_lower_bound": no_bridge["maximum_H"] + 2,
            "even_level_two_branch_uses_proved_Type_I_15_750": True,
            "odd_bound_uses_the_whole_preexisting_residual_band_and_two_new_layers": True,
            "all_size_localization_proved": False, "eventual_E1_proved": False,
            "global_conclusions_remain_conditional_on_missing_all_size_bridge": True,
            "proved": True}


def proposition_15774() -> dict[str, object]:
    capacity = [capacity_exclusion(p, r) for p in (29, 31, 37, 43) for r in (3, 4, 5)]
    residual = [residual_two_layer_exclusion(p, s) for p in (29, 31, 37, 43) for s in (1, 2)]
    survivors = [first_scalar_survivor(37, r) for r in (3, 4, 5)]
    minimal = [minimal_four_gap_consequences(p, [row for row in residual if row["p"] == p])
               for p in (29, 31, 37, 43)]
    return {"prop": "15.774", "status": "PROVED_INFINITE_FAMILY",
            "proof_note": "evidence/NOTE_2026-09-04_SMALL_MASS_TWO_TYPE_BRIDGE.md",
            "local_note": "evidence/NOTE_2026-09-04_SHARP_SMALL_MASS_SPECTRUM.md",
            "residual_range": "p>=29,k in {5p+1,5p+3}",
            "new_generic_frontier": "p>=29,t>=q+3,k>=5p+5, q=(p-1)/2",
            "eventual_capacity_H_bounds": {"r3": "5p+4", "r4": "6p+4", "r5": "7p+6"},
            "eventual_capacity_minimum_prime": 37,
            "capacity_records": capacity, "residual_records": residual,
            "minimal_four_gap_bounds": minimal,
            "first_scalar_survivors": survivors,
            "records_are_identity_replays_not_a_prime_census": True,
            "proved": all(r["proved"] for r in capacity + residual + survivors + minimal),
            "residual_ii_closed_general": False, "minimal_four_gap_bridge_closed_general": False,
            "eventual_E1_proved": False, "e1_closed_general": False,
            "original_MO_limit_closed": False}


if __name__ == "__main__":
    result = proposition_15774()
    write_json_atomic(ROOT / "evidence/e1_gmin_m4_prop15774.json", result)
    print(json.dumps({"prop": "15.774", "proved": result["proved"]}, sort_keys=True))
