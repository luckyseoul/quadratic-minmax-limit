#!/usr/bin/env python3
"""Prop. 15.674 -- close the full infinity-plus-(p-2) shell for p>=17.

Proposition 15.673 assumed every directional odd-fibre count was an endpoint.
That hypothesis is unnecessary.  Put ``P=p+1=2m`` and write the exact means
in one quadratic direction type as

    a_d = 2u + P*k_d,          sum_d k_d = m-u.

Every directional floor is at least ``P-2``.  For ``1<=u<=m-2`` every
``k_d`` is at least one, contradicting ``sum k_d=m-u<m``.  At ``u=0`` the
two genuine Proposition 15.723 floor-plus-two cells at ``p=17`` are retained;
their quotient two still cannot fit the total quotient sum ``m``.  Thus all
``m`` directions have floor ``P`` and mean ``P``.  At ``u=m-1``, exactly
``m-1`` directions have floor and mean
``P-2`` and one arbitrary direction has mean ``2P-2``.  All intermediate
odd-fibre counts have floors in ``(P,2P-2]``, so they can occur only as that
single exceptional direction.

The endpoint baseline directions in a type are homogeneous when both kinds
have the same floor: otherwise equal means force equal parallel counts while
their coefficient congruences differ by one modulo ``q=(p-1)/2``.  Two
``b=1`` baseline types violate the pair-deficit budget, while two
``b=p-2`` baseline types leave at most two determined directions and hence
are collinear, already closed by Propositions 15.671--15.672.  Thus every
branch again has one baseline type of each endpoint kind and exactly the
same four arithmetic rows as Proposition 15.673.  Their three uniform
contradictions and the exact ``p=17`` norm obstruction therefore close the
entire shell, including every intermediate odd-fibre profile.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15723 import floor_excess_admissible
from e1_gmin_m4_prop15673 import (
    branch_arithmetic,
    branch_name,
    coefficient_ledger,
    p17_xnor_l1_ledger,
    parameters,
    symbolic_range_ledger,
)


ROOT = Path(__file__).resolve().parents[1]


def full_profile_floor_ledger(p: int, phase: int) -> dict[str, object]:
    """Partition every allowed odd fibre count by its exact floor band."""
    q, m, period = parameters(p)
    s = p - 2
    floors = {
        b: full_symbolic_floor(p, b, phase) for b in range(1, s + 1, 2)
    }
    floor_p = tuple(b for b, value in floors.items() if value == period)
    floor_p_minus_2 = tuple(
        b for b, value in floors.items() if value == period - 2
    )
    intermediate = tuple(b for b in floors if b not in (1, s))
    expected = {
        "p3_phase0": ((1, s), ()),
        "p1_phase1": ((s,), (1,)),
        "p3_phase1": ((), (1, s)),
        "p1_phase0": ((1,), (s,)),
    }[branch_name(p, phase)]
    if (floor_p, floor_p_minus_2) != expected:
        raise ArithmeticError("endpoint floor bands changed")
    if not all(
        period < floors[b] <= 2 * period - 2
        and 2 * period - 2 - floors[b] != 2
        for b in intermediate
    ):
        raise ArithmeticError("an intermediate floor left the residue band")
    return {
        "p": p,
        "q": q,
        "m": m,
        "P": period,
        "phase": phase,
        "branch": branch_name(p, phase),
        "finite_boundary_size": s,
        "floors": floors,
        "floor_P_counts": floor_p,
        "floor_P_minus_2_counts": floor_p_minus_2,
        "intermediate_counts": intermediate,
        "intermediate_floor_range": [
            min((floors[b] for b in intermediate), default=None),
            max((floors[b] for b in intermediate), default=None),
        ],
        "all_intermediate_floors_strictly_above_P": True,
        "all_intermediate_floors_at_most_2P_minus_2": True,
        "proved": True,
    }


def residue_classification_ledger(p: int, phase: int) -> dict[str, object]:
    """Prove that only the saturated and one-exception type forms survive."""
    floor = full_profile_floor_ledger(p, phase)
    m = int(floor["m"])
    period = int(floor["P"])
    rows = []
    for u in range(m):
        residue = 2 * u
        quotient_sum = m - u

        def minimum_k(b: int | None, direction_floor: int) -> int:
            k = max(0, (direction_floor - residue + period - 1) // period)
            excess = residue + period * k - direction_floor
            admissible = (
                excess >= 0 and excess != 2
                if b is None
                else floor_excess_admissible(p, b, phase, excess)
            )
            if not admissible:
                k += 1
            return k

        # The two 15.723 exceptions lie in the middle band, never at either
        # endpoint baseline floor.
        low_k = minimum_k(None, period - 2)
        high_k = minimum_k(None, period)
        intermediate_k = min(
            minimum_k(int(b), int(value))
            for b, value in floor["floors"].items()
            if int(b) not in (1, p - 2)
        )
        possible = bool(
            (u == 0 and floor["floor_P_counts"])
            or (u == m - 1 and floor["floor_P_minus_2_counts"])
        )
        rows.append(
            {
                "u": u,
                "residue": residue,
                "required_quotient_sum": quotient_sum,
                "minimum_k_by_floor_class": {
                    "P-2": low_k,
                    "P": high_k,
                    "intermediate": intermediate_k,
                },
                "possible_in_this_branch": possible,
            }
        )

    zero_minima = rows[0]["minimum_k_by_floor_class"]
    if not (
        zero_minima["P-2"] == 2
        and zero_minima["P"] == 1
        and zero_minima["intermediate"] >= 2
    ):
        raise ArithmeticError("u=0 quotient minima changed")
    if rows[-1]["minimum_k_by_floor_class"] != {
        "P-2": 0,
        "P": 1,
        "intermediate": 1,
    }:
        raise ArithmeticError("u=m-1 quotient minima changed")
    if not all(
        row["minimum_k_by_floor_class"]["P-2"] >= 1
        and row["minimum_k_by_floor_class"]["P"] >= 1
        and row["minimum_k_by_floor_class"]["intermediate"] >= 1
        and row["required_quotient_sum"] < m
        for row in rows[1:-1]
    ):
        raise ArithmeticError("an interior residue was not excluded")

    return {
        "p": p,
        "phase": phase,
        "branch": str(floor["branch"]),
        "exact_type_form": "a_d=2u+P*k_d, sum k_d=m-u",
        "residue_rows": rows,
        "only_possible_residues_before_branch_floor_availability": [
            0,
            period - 2,
        ],
        "surviving_residues_in_this_branch": [
            row["residue"] for row in rows if row["possible_in_this_branch"]
        ],
        "saturated_form": (
            "u=0: all m directions have floor P and mean P"
        ),
        "one_exception_form": (
            "u=m-1: m-1 directions have floor/mean P-2 and one arbitrary "
            "direction has mean 2P-2"
        ),
        "intermediate_direction_limit_per_type": 1,
        "floor_plus_two_cells_retained": sorted(
            [int(b), phase]
            for b, value in floor["floors"].items()
            if int(b) not in (1, p - 2)
            and floor_excess_admissible(
                p,
                int(b),
                phase,
                2,
            )
        ),
        "proved": True,
    }


def full_profile_type_pair_ledger(p: int, phase: int) -> dict[str, object]:
    """Use incidence geometry to force opposite endpoint baseline kinds."""
    floor = full_profile_floor_ledger(p, phase)
    residue = residue_classification_ledger(p, phase)
    m = int(floor["m"])
    s = p - 2
    branch = str(floor["branch"])
    forms = {
        "p3_phase0": {
            "b1": {"baseline_count": m, "exceptions": 0},
            "complement": {"baseline_count": m, "exceptions": 0},
        },
        "p1_phase1": {
            "b1": {"baseline_count": m - 1, "exceptions": 1},
            "complement": {"baseline_count": m, "exceptions": 0},
        },
        "p3_phase1": {
            "b1": {"baseline_count": m - 1, "exceptions": 1},
            "complement": {"baseline_count": m - 1, "exceptions": 1},
        },
        "p1_phase0": {
            "b1": {"baseline_count": m, "exceptions": 0},
            "complement": {"baseline_count": m - 1, "exceptions": 1},
        },
    }[branch]
    b1_pair_count = 2 * int(forms["b1"]["baseline_count"])
    b1_pair_deficit = b1_pair_count * (s - 1)
    pair_budget = s * (s - 1)
    complement_pair_determined_upper = 2 * int(forms["complement"]["exceptions"])
    offset = int(forms["b1"]["exceptions"]) + int(
        forms["complement"]["exceptions"]
    )
    expected_offset = {
        "p3_phase0": 0,
        "p1_phase1": 1,
        "p3_phase1": 2,
        "p1_phase0": 1,
    }[branch]
    proved = bool(
        residue["proved"]
        and b1_pair_deficit > pair_budget
        and complement_pair_determined_upper <= 2
        and offset == expected_offset
    )
    if not proved:
        raise ArithmeticError("full-profile type pairing failed")
    return {
        "p": p,
        "phase": phase,
        "branch": branch,
        "forms": forms,
        "equal_floor_baseline_homogeneity": (
            "equal means force equal P_d, but b=1 and b=p-2 baseline "
            "coefficient congruences differ by one modulo q"
        ),
        "two_b1_types": {
            "b1_baseline_directions": b1_pair_count,
            "required_deficit": b1_pair_deficit,
            "pair_budget": pair_budget,
            "contradiction": True,
        },
        "two_complement_types": {
            "maximum_determined_directions": complement_pair_determined_upper,
            "at_most_two_directions_implies_collinear": True,
            "collinear_case_closed_by": "15.671--15.672",
            "contradiction": True,
        },
        "forced_pair": ["b=1 baseline type", "b=p-2 baseline type"],
        "finite_edge_offset": offset,
        "intermediate_directions_total_upper_bound": offset,
        "coefficient_ledger": coefficient_ledger(p, phase),
        "proved": True,
    }


def full_profile_branch_exclusion(p: int, phase: int) -> dict[str, object]:
    """Apply Proposition 15.673's arithmetic to the full shell."""
    pair = full_profile_type_pair_ledger(p, phase)
    arithmetic = branch_arithmetic(p, phase)
    candidates = list(arithmetic["candidates"])
    branch = str(arithmetic["branch"])
    if branch == "p1_phase0" and p == 17:
        l1 = p17_xnor_l1_ledger()
        excluded = bool(
            len(candidates) == 1
            and candidates[0]["x"] == 0
            and candidates[0]["y"] == 7
            and l1["contradiction"]
        )
        method = "exact complementary-baseline inter-fibre l1 minimum"
    else:
        l1 = None
        excluded = not candidates or all(
            bool(row["boundary_contradiction"]) for row in candidates
        )
        method = (
            "no congruence candidate"
            if not candidates
            else "I<=|S|+2E boundary support inequality"
        )
    return {
        "p": p,
        "phase": phase,
        "c_H": 1 if phase == 0 else -1,
        "excluded": excluded,
        "method": method,
        "full_profile_type_pair_ledger": pair,
        "arithmetic": arithmetic,
        "p17_l1": l1,
    }


def theorem_record() -> dict[str, object]:
    samples = {
        str(p): {
            str(phase): full_profile_branch_exclusion(p, phase)
            for phase in (0, 1)
        }
        for p in (17, 19, 23, 29, 31, 37, 41, 101)
    }
    symbolic = symbolic_range_ledger()
    proved = bool(
        symbolic["covers_every_odd_p_at_least_17"]
        and all(
            row["excluded"]
            for by_phase in samples.values()
            for row in by_phase.values()
        )
    )
    return {
        "prop": "15.674",
        "title": "Complete infinity-plus-(p-2) boundary-shell exclusion",
        "proved": proved,
        "theorem": {
            "boundary": "infinity plus p-2 finite points",
            "all_directional_odd_fibre_profiles": True,
            "all_odd_primes_p_at_least_17": "EXCLUDED_FOR_BOTH_PRODUCT_SIGNS",
            "dependency": "15.671--15.673",
            "larger_infinity_present_boundaries": "OPEN",
            "all_finite_large_boundaries": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "generic_residue_argument": (
            "only residues 0 and p-1 survive; every intermediate direction "
            "is the unique mean-2p exception of its type"
        ),
        "symbolic_arithmetic_ledger": symbolic,
        "samples": samples,
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.674 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15674.json"
    destination.write_text(json.dumps(record, indent=2) + "\n")
    print("Prop 15.674 full infinity-plus-(p-2) shell: proved")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
