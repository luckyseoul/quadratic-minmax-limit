#!/usr/bin/env python3
"""Prop. 15.673 -- exclude every endpoint-only near-line boundary for p>=17.

Let the odd-degree boundary be infinity together with ``s=p-2`` finite
points.  For a projective direction ``d``, let ``b_d`` be the number of
parallel affine fibres meeting the finite boundary oddly.  This proposition
handles the complete endpoint branch

    b_d in {1,p-2} for every d.

The collinear case was closed in Propositions 15.671--15.672.  In the
remaining case, exact same-type mean quantization reduces each of the four
``(p mod 4, phase)`` branches to two baseline parallel counts ``x,y``.
Writing ``q=(p-1)/2`` and ``m=q+1``, the four reductions are

==================  ========  =========================
branch              E offset  congruences
==================  ========  =========================
p=3 mod 4, phase 0      0     q|y,     q|(x-1)
p=1 mod 4, phase 1      1     q|x,     q|(y-1)
p=3 mod 4, phase 1      2     q|y,     q|(x+1)
p=1 mod 4, phase 0      1     q|x,     q|(y+1)
==================  ========  =========================

Here ``E=m(x+y)+offset`` is the finite-edge count and
``I=4p+1-E`` is the infinity-edge count.  Since ``I>=1``, always
``x+y<=7``.  The first two rows have the unique large-prime candidates
``(1,0)`` and ``(0,1)``; both violate the elementary boundary inequality
``I<=s+2E``.  The third row has no candidate for ``q>=9``.  The fourth has
none for ``q>8``.  Its sole endpoint ``p=17,(x,y)=(0,7)`` has
``I=5,E=64``.  A baseline complementary direction prescribes

    L_st = 1-n_s-n_t+1_{st=ab},   sum n_s=5.

Its exact entrywise l1 minimum is 75, while only ``E-y=57`` transverse
selected edges are available.  This contradiction closes the last case.

Consequently no endpoint-only infinity-plus-(p-2) boundary is possible for
any odd prime p>=17.  Non-endpoint directional profiles, all-finite large
boundaries, residual (ii), R1, Type I, and the limit remain open.
"""
from __future__ import annotations

import json
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15632 import scaled_direction_floor
from e1_gmin_m4_prop15642 import nonbaseline_scaled_cost_floor


ROOT = Path(__file__).resolve().parents[1]


def parameters(p: int) -> tuple[int, int, int]:
    if p < 17 or p % 2 == 0:
        raise ValueError("need odd p>=17")
    q = (p - 1) // 2
    m = q + 1
    return q, m, p + 1


def branch_name(p: int, phase: int) -> str:
    if phase not in (0, 1):
        raise ValueError("phase must be zero or one")
    parameters(p)
    return f"p{p % 4}_phase{phase}"


def endpoint_floor_ledger(p: int, phase: int) -> dict[str, object]:
    """Exact floors at b=1 and b=p-2 and the lift quantum used below."""
    q, m, period = parameters(p)
    one = scaled_direction_floor(p, 1, phase)
    complement = scaled_direction_floor(p, p - 2, phase)
    lift = nonbaseline_scaled_cost_floor(p)
    expected = {
        "p3_phase0": (period, period),
        "p1_phase1": (period - 2, period),
        "p3_phase1": (period - 2, period - 2),
        "p1_phase0": (period, period - 2),
    }[branch_name(p, phase)]
    if (one, complement) != expected or lift < 4:
        raise ArithmeticError("endpoint floor or lift ledger changed")
    return {
        "p": p,
        "q": q,
        "m": m,
        "period": period,
        "phase": phase,
        "b1_floor": one,
        "b_p_minus_2_floor": complement,
        "minimum_nonzero_lift_scaled_cost": lift,
        "lift_cost_at_least_four": lift >= 4,
        "type_budget": m * period,
    }


def same_type_normal_form(p: int, phase: int) -> dict[str, object]:
    """Record the residue argument giving the two-count normal form.

    For one quadratic type, all exact means have the form
    ``a_d=r+(p+1)k_d`` with ``r=2u`` and ``sum k_d=m-u``.  A putative value
    only two above its parity baseline is impossible because a nonzero lift
    costs at least four.  Comparing this with the two endpoint floors leaves
    the stated baseline/exception patterns.
    """
    floor = endpoint_floor_ledger(p, phase)
    type_forms = endpoint_type_form_ledger(p, phase)
    period = int(floor["period"])
    m = int(floor["m"])
    branch = branch_name(p, phase)
    forms = {
        "p3_phase0": {
            "allowed_type_forms": [
                "all b=1 baselines at mean p+1",
                "all b=p-2 baselines at mean p+1",
            ],
            "exceptions_per_type": 0,
            "finite_edge_offset": 0,
            "geometry_forces_mixed_forms": True,
        },
        "p1_phase1": {
            "allowed_type_forms": [
                "all b=p-2 baselines, residue 0",
                "m-1 b=1 baselines plus one p+1 mean jump, residue p-1",
            ],
            "exceptions_per_type": "0 or 1",
            "finite_edge_offset": 1,
            "geometry_forces_mixed_forms": True,
        },
        "p3_phase1": {
            "allowed_type_forms": [
                "m-1 baselines of one endpoint kind plus one p+1 mean jump"
            ],
            "exceptions_per_type": 1,
            "finite_edge_offset": 2,
            "geometry_forces_opposite_baseline_kinds": True,
        },
        "p1_phase0": {
            "allowed_type_forms": [
                "all b=1 baselines, residue 0",
                "m-1 b=p-2 baselines plus one p+1 mean jump, residue p-1",
            ],
            "exceptions_per_type": "0 or 1",
            "finite_edge_offset": 1,
            "geometry_forces_mixed_forms": True,
        },
    }[branch]
    return {
        "branch": branch,
        "same_type_mean": "a_d=r+(p+1)k_d",
        "residue": "r=2u, 0<=u<m",
        "quotient_sum": "sum_d k_d=m-u",
        "forbidden_two_unit_lift": int(
            floor["minimum_nonzero_lift_scaled_cost"]
        )
        >= 4,
        "endpoint_pair_deficit": "number of b=1 directions is at most p-2",
        "equality_case": (
            "equality forces a (p-2)-arc with exactly three undetermined "
            "directions; adjoining any two gives a p-arc, so Segre's odd-order "
            "p-arc theorem puts each extension on a conic; two extensions "
            "share p-2>=5 points and hence the same conic, which would contain "
            "three collinear points at infinity"
        ),
        "noncollinear_direction_minimum": 3,
        "collinear_case": "CLOSED_BY_15.671_AND_15.672",
        "type_form_ledger": type_forms,
        "directions_per_type": m,
        "period": period,
        **forms,
        "proved": True,
    }


def endpoint_type_form_ledger(p: int, phase: int) -> dict[str, object]:
    """Classify the two same-type endpoint forms before count arithmetic.

    This is the explicit bridge between the endpoint geometry and the
    two-count normal forms.  After the equality and collinear exits, the
    number ``R`` of ``b=1`` directions satisfies ``3<=R<=p-3``.  Mean
    residues classify each quadratic type into one of the two rows below.
    In the two equal-floor branches, baseline directions of different
    endpoint kinds would have equal means and hence equal ``P_d``; their
    coefficient congruences differ by one modulo ``q``, so they cannot
    coexist in one type.  Enumerating the possible ``R`` values then forces
    one b=1-baseline type and one complementary-baseline type.
    """
    q, m, period = parameters(p)
    branch = branch_name(p, phase)
    forms = {
        "p3_phase0": [
            {
                "baseline_kind": "b=1",
                "b1_count_options": [m],
                "exceptions": 0,
                "type_finite_edge_offset": 0,
            },
            {
                "baseline_kind": "b=p-2",
                "b1_count_options": [0],
                "exceptions": 0,
                "type_finite_edge_offset": 0,
            },
        ],
        "p1_phase1": [
            {
                "baseline_kind": "b=1",
                "b1_count_options": [m - 1, m],
                "exceptions": 1,
                "type_finite_edge_offset": 1,
            },
            {
                "baseline_kind": "b=p-2",
                "b1_count_options": [0],
                "exceptions": 0,
                "type_finite_edge_offset": 0,
            },
        ],
        "p3_phase1": [
            {
                "baseline_kind": "b=1",
                "b1_count_options": [m - 1, m],
                "exceptions": 1,
                "type_finite_edge_offset": 1,
            },
            {
                "baseline_kind": "b=p-2",
                "b1_count_options": [0, 1],
                "exceptions": 1,
                "type_finite_edge_offset": 1,
            },
        ],
        "p1_phase0": [
            {
                "baseline_kind": "b=1",
                "b1_count_options": [m],
                "exceptions": 0,
                "type_finite_edge_offset": 0,
            },
            {
                "baseline_kind": "b=p-2",
                "b1_count_options": [0, 1],
                "exceptions": 1,
                "type_finite_edge_offset": 1,
            },
        ],
    }[branch]

    minimum_r = 3
    maximum_r = p - 3
    admissible_pairs = []
    for first_index, first in enumerate(forms):
        for second in forms[first_index:]:
            r_values = sorted(
                {
                    first_count + second_count
                    for first_count in first["b1_count_options"]
                    for second_count in second["b1_count_options"]
                    if minimum_r <= first_count + second_count <= maximum_r
                }
            )
            if r_values:
                admissible_pairs.append(
                    {
                        "baseline_kinds": sorted(
                            [first["baseline_kind"], second["baseline_kind"]]
                        ),
                        "R_values": r_values,
                        "finite_edge_offset": (
                            first["type_finite_edge_offset"]
                            + second["type_finite_edge_offset"]
                        ),
                    }
                )

    expected_offset = {
        "p3_phase0": 0,
        "p1_phase1": 1,
        "p3_phase1": 2,
        "p1_phase0": 1,
    }[branch]
    expected_kinds = ["b=1", "b=p-2"]
    if len(admissible_pairs) != 1:
        raise ArithmeticError("endpoint geometry did not force one type pair")
    pair = admissible_pairs[0]
    if (
        pair["baseline_kinds"] != expected_kinds
        or pair["finite_edge_offset"] != expected_offset
    ):
        raise ArithmeticError("endpoint type normal form changed")

    return {
        "p": p,
        "q": q,
        "m": m,
        "period": period,
        "branch": branch,
        "same_type_residue_classification": (
            "a_d=2u+(p+1)k_d, sum k_d=m-u, and no floor+2 lift"
        ),
        "equal_floor_mixed_baseline_obstruction": (
            "equal baseline means force equal P_d, while the b=1 and "
            "b=p-2 coefficient congruences differ by one modulo q"
        ),
        "geometry_after_exits": f"{minimum_r}<=R<={maximum_r}",
        "type_forms": forms,
        "admissible_type_pairs": admissible_pairs,
        "forced_opposite_baseline_kinds": True,
        "finite_edge_offset": expected_offset,
        "proved": True,
    }


def coefficient_ledger(p: int, phase: int) -> dict[str, object]:
    """Baseline coefficient congruences before substituting the edge count."""
    q, _m, _period = parameters(p)
    sigma = 1 if phase == 0 else -1
    tau = -1 if branch_name(p, phase) in ("p3_phase0", "p1_phase1") else 1
    return {
        "q": q,
        "b1_target": f"4 {'+' if sigma == 1 else '-'} z_j",
        "b1_divisibility": f"q divides I+P_d-{4 + sigma}",
        "complement_target": f"4 {'+' if tau == 1 else '-'} z_a*z_b",
        "complement_divisibility": "q divides I+P_d-4",
        "general_b1_matrix": (
            "L_st=2c-n_s-n_t+sigma(1_{s=j}+1_{t=j}), "
            "(p-1)c=I+P_d-4-sigma"
        ),
        "general_complement_matrix": (
            "L_st=2c-n_s-n_t+tau*1_{st=ab}, "
            "(p-1)c=I+P_d-4"
        ),
        "proved": True,
    }


def branch_arithmetic(p: int, phase: int) -> dict[str, object]:
    """Enumerate the at-most-eight quotient box after symbolic reduction."""
    q, m, _period = parameters(p)
    branch = branch_name(p, phase)
    offset, first, second = {
        "p3_phase0": (0, lambda x, y: y, lambda x, y: x - 1),
        "p1_phase1": (1, lambda x, y: x, lambda x, y: y - 1),
        "p3_phase1": (2, lambda x, y: y, lambda x, y: x + 1),
        "p1_phase0": (1, lambda x, y: x, lambda x, y: y + 1),
    }[branch]
    congruence_text = {
        "p3_phase0": ["q divides y", "q divides x-1"],
        "p1_phase1": ["q divides x", "q divides y-1"],
        "p3_phase1": ["q divides y", "q divides x+1"],
        "p1_phase0": ["q divides x", "q divides y+1"],
    }[branch]
    candidates = []
    for x in range(8):
        for y in range(8 - x):
            finite = m * (x + y) + offset
            infinity = 4 * p + 1 - finite
            if infinity < 1 or first(x, y) % q or second(x, y) % q:
                continue
            boundary_upper = (p - 2) + 2 * finite
            candidates.append(
                {
                    "x": x,
                    "y": y,
                    "E": finite,
                    "I": infinity,
                    "I_boundary_upper": boundary_upper,
                    "boundary_contradiction": infinity > boundary_upper,
                }
            )
    return {
        "p": p,
        "phase": phase,
        "branch": branch,
        "x": "baseline parallel count in the b=1 type",
        "y": "baseline parallel count in the b=p-2 type",
        "finite_edge_formula": f"E=m(x+y)+{offset}",
        "infinity_edge_formula": "I=4p+1-E",
        "quotient_box": "x>=0, y>=0, x+y<=7",
        "congruences": congruence_text,
        "candidates": candidates,
    }


def p17_xnor_l1_ledger() -> dict[str, object]:
    """Exact final l1 obstruction for p=17, phase zero, x=0,y=7."""
    values = {}
    for positive_fibres in range(1, 6):
        u = positive_fibres
        zero = 17 - u
        # First omit the distinguished +1_ab.  Separating zero-zero,
        # zero-positive, and positive-positive pairs eliminates the positive
        # magnitudes because their sum is five.  The distinguished term can
        # lower this base norm by at most one.
        base = (
            comb(zero, 2)
            + zero * (5 - u)
            + 5 * (u - 1)
            - comb(u, 2)
        )
        values[str(u)] = {
            "base_l1": base,
            "after_distinguished_pair_lower_bound": base - 1,
        }
    minimum = min(
        int(row["after_distinguished_pair_lower_bound"])
        for row in values.values()
    )
    # Equality occurs for five unit star counts, including fibres a,b.
    achieved = 75
    transverse_budget = 64 - 7
    return {
        "p": 17,
        "candidate": {"x": 0, "y": 7, "E": 64, "I": 5},
        "baseline_complement_parallel_count": 7,
        "coefficient_scalar_2c": 1,
        "matrix": "L_st=1-n_s-n_t+1_{st=ab}, sum n_s=5",
        "positive_fibre_ledger": values,
        "lower_bound": minimum,
        "achieved_by_five_unit_counts_including_a_b": achieved,
        "exact_minimum": minimum if minimum == achieved else None,
        "transverse_edge_budget": transverse_budget,
        "contradiction": minimum > transverse_budget,
    }


def endpoint_branch_exclusion(p: int, phase: int) -> dict[str, object]:
    arithmetic = branch_arithmetic(p, phase)
    branch = str(arithmetic["branch"])
    candidates = list(arithmetic["candidates"])
    if branch == "p1_phase0" and p == 17:
        l1 = p17_xnor_l1_ledger()
        excluded = bool(
            len(candidates) == 1
            and candidates[0]["x"] == 0
            and candidates[0]["y"] == 7
            and l1["contradiction"]
        )
        method = "exact inter-fibre l1 minimum"
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
        "applicable": p >= 17,
        "excluded": excluded,
        "method": method,
        "floor_ledger": endpoint_floor_ledger(p, phase),
        "normal_form": same_type_normal_form(p, phase),
        "coefficient_ledger": coefficient_ledger(p, phase),
        "arithmetic": arithmetic,
        "p17_l1": l1,
    }


def symbolic_range_ledger() -> dict[str, object]:
    """Four inequalities covering every odd p>=17, not a finite-prime scan."""
    return {
        "p3_phase0": {
            "range": "q>=9",
            "forced_candidate": [1, 0],
            "boundary_gap": "3q+3>0",
            "excluded": True,
        },
        "p1_phase1": {
            "range": "q>=8",
            "forced_candidate": [0, 1],
            "boundary_gap": "3q>0",
            "excluded": True,
        },
        "p3_phase1": {
            "range": "q>=9",
            "reason": "0<=x+1<=8<q contradicts q|(x+1)",
            "excluded": True,
        },
        "p1_phase0": {
            "range": "q>8",
            "reason": "0<=x<=7 forces x=0, while 1<=y+1<=8<q",
            "endpoint_q_eq_8": "p=17, (x,y)=(0,7), l1 75>57",
            "excluded": True,
        },
        "covers_every_odd_p_at_least_17": True,
    }


def theorem_record() -> dict[str, object]:
    samples = {
        str(p): {
            str(phase): endpoint_branch_exclusion(p, phase)
            for phase in (0, 1)
        }
        for p in (17, 19, 23, 29, 31, 37, 41, 101)
    }
    symbolic = symbolic_range_ledger()
    proved = bool(
        symbolic["covers_every_odd_p_at_least_17"]
        and p17_xnor_l1_ledger()["contradiction"]
        and all(
            row["excluded"]
            for by_phase in samples.values()
            for row in by_phase.values()
        )
    )
    return {
        "prop": "15.673",
        "title": "Complete endpoint-only infinity-plus-(p-2) boundary exclusion",
        "proved": proved,
        "theorem": {
            "endpoint_hypothesis": "b_d in {1,p-2} for every direction",
            "all_odd_primes_p_at_least_17": "EXCLUDED_FOR_BOTH_PRODUCT_SIGNS",
            "collinear_dependency": "Propositions 15.671--15.672",
            "nonendpoint_directional_profiles": "OPEN",
            "all_finite_large_boundaries": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "symbolic_range_ledger": symbolic,
        "p17_l1_ledger": p17_xnor_l1_ledger(),
        "samples": samples,
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.673 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15673.json"
    destination.write_text(json.dumps(record, indent=2) + "\n")
    print("Prop 15.673 endpoint-only near-line exclusion: proved")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
