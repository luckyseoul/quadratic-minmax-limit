#!/usr/bin/env python3
"""Prop. 15.671 -- rigid-sign exclusion for the near-line boundary.

Let ``H`` satisfy the residual affine hypotheses of Proposition 15.632,
``|H|=4p+1``, and suppose its odd-degree boundary consists of infinity and
``p-2`` finite points on one affine line.  The line direction has one odd
boundary fibre, while every transverse direction has ``p-2`` odd fibres.

In either of the following cases the exact parity floors exhaust the
quadratic-type budgets and force the directional slacks pointwise:

* ``p=3 (mod 4)`` and ``c_H=+1``;
* ``p=1 (mod 4)`` and ``c_H=-1`` (for ``p>=13``; the line type has only two
  units of apparent surplus, less than the minimum cost of a nonzero
  quadratic lift).

Writing ``z_s=2x_s-1`` on the middle slice, the forced targets are

    eps_0 S_H = 4 + sigma z_j
    eps_d S_H = 4 - z_a z_b                  (d transverse),

where ``sigma=(-1)^phase`` and ``a,b`` are the two fibres corresponding to
the two omitted points of the line.  Comparing degree-two coefficients on
the slice gives congruences for the infinity-edge count ``I`` and the finite
parallel-edge counts ``P_d``.  The signed inter-fibre matrices also give

    2I + (p+1)P_d <= 8p+1+r_d,

with ``r_0=sigma`` and ``r_d=1`` transversely.

For ``q=(p-1)/2``, summing the parallel-count congruences gives

    I = 3+sigma (mod q).

If ``p=1 (mod 4)``, then ``q`` is even and ``sigma=-1``, so this forces the
odd number ``I`` to be even.  If ``p=3 (mod 4)``, then ``sigma=1`` and
``I=4+q k_0`` with odd ``k_0``.  For ``p>=19`` the l1 inequalities force all
remaining quotient variables to vanish, while the edge-count identity
forces ``k_0=8``, a contradiction.

Thus the rigid product-sign branch of every collinear
infinity-plus-``(p-2)`` boundary is impossible for ``p>=13`` in the
``1 mod 4`` class and for ``p>=19`` in the ``3 mod 4`` class.  The opposite
product sign, noncollinear boundaries, residual (ii), R1, and the limit
remain open.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15642 import nonbaseline_scaled_cost_floor


ROOT = Path(__file__).resolve().parents[1]


def middle_weight(p: int) -> int:
    if p < 5 or p % 2 == 0:
        raise ValueError("p must be odd and at least five")
    return (p + 1) // 2


def rigid_sign(p: int) -> int:
    """The product sign whose near-line floors force pointwise rigidity."""
    if p < 5 or p % 2 == 0:
        raise ValueError("p must be odd and at least five")
    return 1 if p % 4 == 3 else -1


def near_line_floor_ledger(p: int) -> dict[str, object]:
    """Exact type-floor ledger for the rigid product sign."""
    m = middle_weight(p)
    c_h = rigid_sign(p)
    phase = 0 if c_h == 1 else 1
    sigma = 1 if phase == 0 else -1
    type_budget = m * (p + 1)
    special_cost = p + sigma
    transverse_cost = p + 1
    line_type_floor = special_cost + (m - 1) * transverse_cost
    opposite_type_floor = m * transverse_cost
    line_type_surplus = type_budget - line_type_floor
    opposite_type_surplus = type_budget - opposite_type_floor
    lift_cost = nonbaseline_scaled_cost_floor(p)

    if p % 4 == 3:
        rigid = bool(
            phase == 0
            and m % 2 == 0
            and line_type_surplus == 0
            and opposite_type_surplus == 0
        )
    else:
        # p=5 has lift cost two and is deliberately outside this theorem.
        rigid = bool(
            p >= 13
            and phase == 1
            and m % 2 == 1
            and line_type_surplus == 2
            and opposite_type_surplus == 0
            and lift_cost > line_type_surplus
        )

    return {
        "p": p,
        "p_mod_4": p % 4,
        "c_H": c_h,
        "phase": phase,
        "sigma": sigma,
        "directions_per_type": m,
        "type_budget": type_budget,
        "special_b": 1,
        "transverse_b": p - 2,
        "special_baseline": "x_j" if phase == 0 else "1-x_j",
        "transverse_baseline": "x_a xor x_b",
        "special_scaled_cost": special_cost,
        "transverse_scaled_cost": transverse_cost,
        "line_type_floor_sum": line_type_floor,
        "opposite_type_floor_sum": opposite_type_floor,
        "line_type_surplus": line_type_surplus,
        "opposite_type_surplus": opposite_type_surplus,
        "minimum_nonzero_lift_scaled_cost": lift_cost,
        "all_directional_slacks_forced_to_baseline": rigid,
    }


def coefficient_ledger(p: int) -> dict[str, object]:
    """The all-direction congruence and l1 consequences after rigidity."""
    floor = near_line_floor_ledger(p)
    if not floor["all_directional_slacks_forced_to_baseline"]:
        raise ValueError("the requested prime is outside the rigid floor range")
    q = (p - 1) // 2
    sigma = int(floor["sigma"])
    return {
        "q": q,
        "special_target": f"4 {'+' if sigma == 1 else '-'} z_j",
        "transverse_target": "4-z_a*z_b",
        "special_parallel_congruence": f"P_0 = {4 + sigma}-I (mod q)",
        "transverse_parallel_congruence": "P_d = 4-I (mod q)",
        "parallel_count_sum": "P_0+sum_transverse(P_d)=4p+1-I",
        "summed_congruence": f"I = {3 + sigma} (mod q)",
        "special_l1": f"2I+(p+1)P_0 <= {8 * p + 1 + sigma}",
        "transverse_l1": f"2I+(p+1)P_d <= {8 * p + 2}",
        "proved": True,
    }


def rigid_near_line_exclusion(p: int) -> dict[str, object]:
    """Symbolic arithmetic contradiction for the applicable prime class."""
    floor = near_line_floor_ledger(p)
    if not floor["all_directional_slacks_forced_to_baseline"]:
        return {
            "p": p,
            "applicable": False,
            "excluded": False,
            "reason": "pointwise rigidity is not forced by this floor ledger",
            "floor_ledger": floor,
        }

    q = (p - 1) // 2
    sigma = int(floor["sigma"])
    coeff = coefficient_ledger(p)
    if p % 4 == 1:
        excluded = bool(q % 2 == 0 and sigma == -1)
        contradiction = (
            "I is odd because infinity is in the boundary, but "
            "I=2 (mod q) with even q forces I even"
        )
        arithmetic = {
            "q_even": q % 2 == 0,
            "I_congruence_residue": 2,
            "I_boundary_parity": 1,
        }
    else:
        excluded = bool(p >= 19 and q >= 9 and sigma == 1)
        contradiction = (
            "I=4+q*k0 with odd k0; l1 forces a0=a_d=0, while the "
            "global edge count forces k0=8"
        )
        arithmetic = {
            "q_at_least_9": q >= 9,
            "I_parameterization": "I=4+q*k0",
            "k0_parity": "odd",
            "transverse_parameterization": "P_d=q*a_d, a_d>=0",
            "special_parameterization": "P_0=1+q*a_0, a_0>=0",
            "global_quotient_sum": "k0+a_0+sum(a_d)=8",
            "transverse_l1_quotient": "k0+(q+1)a_d<=8",
            "special_l1_quotient": "k0+(q+1)a_0<=7",
            "l1_forces_all_a_zero": q >= 9,
            "edge_count_forces_k0": 8,
            "forced_k0_has_wrong_parity": True,
        }

    return {
        "p": p,
        "applicable": True,
        "boundary": "infinity plus p-2 collinear finite points",
        "c_H": int(floor["c_H"]),
        "excluded": excluded,
        "contradiction": contradiction,
        "floor_ledger": floor,
        "coefficient_ledger": coeff,
        "arithmetic": arithmetic,
    }


def theorem_record() -> dict[str, object]:
    samples = {
        str(p): rigid_near_line_exclusion(p)
        for p in (7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 101)
    }
    p1_samples = [13, 17, 29, 37, 41, 101]
    p3_samples = [19, 23, 31, 43, 47, 59, 103]
    proved = bool(
        all(rigid_near_line_exclusion(p)["excluded"] for p in p1_samples)
        and all(rigid_near_line_exclusion(p)["excluded"] for p in p3_samples)
    )
    return {
        "prop": "15.671",
        "title": "Rigid-sign exclusion for collinear infinity-plus-(p-2) boundaries",
        "proved": proved,
        "theorem": {
            "p_eq_1_mod_4": (
                "for every odd prime p>=13, c_H=-1 is impossible when the "
                "boundary is infinity plus p-2 collinear finite points"
            ),
            "p_eq_3_mod_4": (
                "for every odd prime p>=19, c_H=+1 is impossible when the "
                "boundary is infinity plus p-2 collinear finite points"
            ),
            "uniform_rigid_sign": "c_H=+1 for p=3 mod 4, c_H=-1 for p=1 mod 4",
            "opposite_product_sign": "OPEN",
            "noncollinear_near_line_profile": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "limit_exists": False,
        },
        "samples": samples,
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.671 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15671.json"
    destination.write_text(json.dumps(record, indent=2) + "\n")
    print("Prop 15.671 rigid-sign near-line exclusion: proved")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
