#!/usr/bin/env python3
"""Prop. 15.675 -- close the first all-finite survivor for p=3,5 mod 8.

Proposition 15.669 closes even all-finite boundary sizes through
``3(p-1)/4`` using directional floors and the pair-deficit budget.  Let
``s`` be the first even integer above that threshold.  This proposition adds
the exact same-type mean congruence

    a_d = 2u + (p+1) k_d,       sum k_d = m-u,

where ``m=(p+1)/2``.  In phase one it forces ``m-1`` directions with
``b=2`` and one middle direction.  In phase zero the unique minimum has
residue four and uses quotient weights zero, one, and two at ``b=0,2,s``.
The resulting exact relaxed gap is

    p mod 8:     1          3          5          7
    gap:      -(p-1)/4   (p+1)/2   (p-1)/2   -(p-7)/4.

Thus the first floor-plus-pair survivor is actually impossible for every
prime ``p>=19`` congruent to 3 or 5 modulo 8.  The other two residue classes,
the next even size, the infinity-present remainder, residual (ii), R1,
global QVAR, Type I, and the limit remain open.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15669 import full_symbolic_floor


ROOT = Path(__file__).resolve().parents[1]


def first_even_survivor(p: int) -> int:
    """First even integer strictly above 3(p-1)/4."""
    if p < 19 or p % 2 == 0:
        raise ValueError("need odd p>=19")
    residue_formula = {
        1: (3 * p + 5) // 4,
        3: (3 * p - 1) // 4,
        5: (3 * p + 1) // 4,
        7: (3 * p + 3) // 4,
    }[p % 8]
    if residue_formula % 2 or not 4 * residue_formula > 3 * (p - 1):
        raise ArithmeticError("first-survivor formula changed")
    if 4 * (residue_formula - 2) > 3 * (p - 1):
        raise ArithmeticError("the preceding even size is still outside")
    return residue_formula


def phase_one_quantized_minimum(p: int) -> dict[str, object]:
    """Exact phase-one deficit after same-type residue quantization."""
    s = first_even_survivor(p)
    m = (p + 1) // 2
    period = p + 1
    deficit = (m - 1) * (s - 2)
    return {
        "p": p,
        "s": s,
        "m": m,
        "P": period,
        "only_residue": period - 2,
        "u": m - 1,
        "quotient_sum": 1,
        "profile": {"b=2": m - 1, "b=s": 1},
        "minimum_deficit": deficit,
        "proof": (
            "for u<=m-2 every k_d>=1 but sum k_d=m-u<m; at u=0 "
            "the floor-P+2 prohibition raises every minimum further"
        ),
        "proved": True,
    }


def phase_zero_quantized_minimum(p: int) -> dict[str, object]:
    """Exact phase-zero deficit after same-type residue quantization."""
    s = first_even_survivor(p)
    m = (p + 1) // 2
    period = p + 1
    # For u>=2, the quotient weights at the deficit-optimal b=0,2,s
    # are 0,1,2.  The deficit increases strictly with u, so u=2 wins.
    u = 2
    b2_count = m & 1
    b0_count = (m + 2 - b2_count) // 2
    middle_count = (m - 2 - b2_count) // 2
    deficit = b0_count * s + b2_count * (s - 2)

    # At u=0, quotient weights 0,1,2,3 occur at b=0,2,4,s.  A weight-three
    # middle direction saves s deficit units, while each residual quotient
    # unit saves only two.  This gives the displayed closed form.
    u0_middle_count, u0_remainder = divmod(m, 3)
    u0_deficit = (m - u0_middle_count) * s - 2 * u0_remainder
    if not u0_deficit > deficit:
        raise ArithmeticError("the u=0 phase-zero branch reached the minimum")

    return {
        "p": p,
        "s": s,
        "m": m,
        "P": period,
        "winning_residue": 4,
        "u": u,
        "quotient_sum": m - u,
        "profile": {
            "b=0": b0_count,
            "b=2": b2_count,
            "b=s": middle_count,
        },
        "minimum_deficit": deficit,
        "u_zero_minimum": u0_deficit,
        "u_zero_strict_gap": u0_deficit - deficit,
        "interior_residue_monotonicity": (
            "D(u+1)-D(u) alternates between s-2 and 2 for u>=2"
        ),
        "proved": True,
    }


def quantized_type_minimum_dp(p: int, phase: int) -> dict[str, object]:
    """Independent relaxed DP over every common residue, b, and lift."""
    if phase not in (0, 1):
        raise ValueError("phase must be zero or one")
    s = first_even_survivor(p)
    period = p + 1
    m = period // 2
    best: tuple[int, int, tuple[int, ...]] | None = None
    for u in range(m):
        residue = 2 * u
        target = m - u
        options = []
        for b in range(0, s + 1, 2):
            floor = full_symbolic_floor(p, b, phase)
            for k in range(target + 1):
                excess = residue + period * k - floor
                if excess >= 0 and excess != 2:
                    options.append((k, s - b, b))
        states: dict[int, tuple[int, tuple[int, ...]]] = {0: (0, ())}
        for _ in range(m):
            next_states: dict[int, tuple[int, tuple[int, ...]]] = {}
            for quotient, (deficit, profile) in states.items():
                for k, added_deficit, b in options:
                    new_quotient = quotient + k
                    if new_quotient > target:
                        continue
                    candidate = (deficit + added_deficit, profile + (b,))
                    old = next_states.get(new_quotient)
                    if old is None or candidate[0] < old[0]:
                        next_states[new_quotient] = candidate
            states = next_states
        if target in states:
            deficit, profile = states[target]
            candidate = (deficit, u, profile)
            if best is None or candidate[0] < best[0]:
                best = candidate
    if best is None:
        raise ArithmeticError("quantized type DP found no state")
    return {
        "p": p,
        "s": s,
        "phase": phase,
        "minimum_deficit": best[0],
        "winning_u": best[1],
        "profile": {
            str(key): value for key, value in sorted(Counter(best[2]).items())
        },
    }


def first_survivor_gap(p: int) -> dict[str, object]:
    """Exact quantized-floor lower bound against the pair budget."""
    s = first_even_survivor(p)
    phase_zero = phase_zero_quantized_minimum(p)
    phase_one = phase_one_quantized_minimum(p)
    required = int(phase_zero["minimum_deficit"]) + int(
        phase_one["minimum_deficit"]
    )
    pair_budget = s * (s - 1)
    gap = required - pair_budget
    expected = {
        1: -(p - 1) // 4,
        3: (p + 1) // 2,
        5: (p - 1) // 2,
        7: -(p - 7) // 4,
    }[p % 8]
    if gap != expected:
        raise ArithmeticError("mod-eight gap formula changed")
    return {
        "p": p,
        "p_mod_8": p % 8,
        "s": s,
        "phase_zero": phase_zero,
        "phase_one": phase_one,
        "required_total_deficit": required,
        "pair_deficit_budget": pair_budget,
        "gap": gap,
        "closed": gap > 0,
    }


def symbolic_mod8_ledger() -> dict[str, object]:
    return {
        "p=1 mod 8": {"gap": "-(p-1)/4", "closed": False},
        "p=3 mod 8": {"gap": "(p+1)/2", "closed": True},
        "p=5 mod 8": {"gap": "(p-1)/2", "closed": True},
        "p=7 mod 8": {"gap": "-(p-7)/4", "closed": False},
        "proof": (
            "substitute the four formulas for the first even integer above "
            "3(p-1)/4 into the two exact quantized deficits"
        ),
    }


def theorem_record() -> dict[str, object]:
    samples = {
        str(p): first_survivor_gap(p)
        for p in (19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 101)
    }
    dp_samples = {
        str(p): {
            str(phase): quantized_type_minimum_dp(p, phase)
            for phase in (0, 1)
        }
        for p in (19, 23, 29, 31, 37, 41, 43)
    }
    dp_agrees = all(
        int(dp_samples[str(p)][str(phase)]["minimum_deficit"])
        == int(
            (
                phase_zero_quantized_minimum(p)
                if phase == 0
                else phase_one_quantized_minimum(p)
            )["minimum_deficit"]
        )
        for p in (19, 23, 29, 31, 37, 41, 43)
        for phase in (0, 1)
    )
    proved = bool(
        dp_agrees
        and all(
            row["closed"] == (int(row["p_mod_8"]) in (3, 5))
            for row in samples.values()
        )
    )
    return {
        "prop": "15.675",
        "title": "Quantized close of the first all-finite survivor",
        "proved": proved,
        "theorem": {
            "all_odd_primes_p_at_least_19_p_mod_8_in_3_5": (
                "the first even s>3(p-1)/4 is excluded"
            ),
            "p_mod_8_in_1_7": "OPEN_AT_THIS_RELAXATION",
            "next_even_boundary_size": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "symbolic_mod8_ledger": symbolic_mod8_ledger(),
        "independent_dp_agrees": dp_agrees,
        "samples": samples,
        "dp_samples": dp_samples,
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.675 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15675.json"
    destination.write_text(json.dumps(record, indent=2) + "\n")
    print("Prop 15.675 first all-finite survivor: closed for p=3,5 mod 8")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
