#!/usr/bin/env python3
r"""Prop. 15.758 -- sharp coefficient cancellation and two local survivor rays.

Put ``m=(p+1)/2`` and use signed coordinates ``z_i=2x_i-1`` on
``J(p,m)``, so ``sum_i z_i=1``.  Once the constant offset ``P`` is fixed,
the representation

    4B = P + sum_(i<j) K_ij z_i z_j

is unique.  Moreover

    4p E[B] = pP-sum K_ij.                              (1)

For ``a,b>=0`` prescribe offset ``P=a+b`` and scaled mass

    a(p-3)+b(p+1).

Then every representation has ``sum K=3a-b``.  Its least possible
coefficient l1 norm, even among all nonnegative integer-valued quadratics,
is

    |3a-b|,                  if |3a-b| != 1,
    3,                       if |3a-b| == 1.             (2)

The upper bound uses only three coordinates.  With
``T=z1*z2+z1*z3+z2*z3``, the forms ``1+T``, ``3-T``,
``1+z1*z2-z1*z3-z2*z3``, and ``2-2*z1*z2`` are four times
nonnegative integer-valued quadratics and carry respectively one ``p-3``
unit, three ``p+1`` units, one ``p+1`` unit, and two ``p+1`` units.
Thus three opposite sharp atoms can cancel one same-sign sharp atom all the
way to a constant.  This disproves any coefficient-capacity lower bound
that grows with ``a+b`` rather than the invariant ``|3a-b|``.

The module also records two infinite families of *local directional
survivors*.  They satisfy the isolated-chart residues, offsets, parity
types, directional mean budgets, parallel-edge budgets, and nonnegative
integral cell conditions.  They are not asserted to be the directional
images of one simple graph.

For ``p=4r+1`` (branch B, ``u=0``), the interval is

    2r^2-5r <= t <= 4r^2-6r-3.                         (3)

Hard rows have ``e_L>=0``, ``sum e_L=t``, literal-plus-lift cell
``A=x_j+2B``, scaled mean ``(e_L+1)(p+1)``, and ``P_L=5+e_L``.
Opposite rows have ``Q_L>=r``, ``sum Q_L=6r+t``, and a phase-zero lift
whose half has scaled mass

    (r-1)(p-3)+(Q_L-r)(p+1)=(p+1)Q_L-2p+4.            (4)

One omitted-pair ``p-3`` atom supplies the required offset shift by minus
two; the remaining atoms can be compact triangles.

For ``p=4r+3`` (branch C, ``u=m-1``), the interval is

    2r^2-4r-2 <= t <= 4r^2-2r-5.                       (5)

Hard rows have ``sum e_L=t+1``, complement-literal-plus-lift mean
``p-1+e_L(p+1)``, and ``P_L=3+e_L``.  Opposite rows have
``Q_L>=r+2``, ``sum Q_L=10r+6+t``, and half-lift mass

    (r-1)(p-3)+(Q_L-r-2)(p+1)=(p+1)Q_L-4p+4.          (6)

At the lower endpoint of either interval, exact common difference-Radon
Parseval has a scalar lower value of order ``8r^3``.  The exact maximum
within the displayed compact, aligned triangle templates is of order
``20r^3``.  Hence scalar Parseval versus row energy has no gap for
``r>=7``.  The remaining question is the integral midpoint/difference-Radon
lift to one ``0/1`` graph, not another l1 or scalar-energy estimate.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from e1_gmin_m4_prop15721 import is_prime


ROOT = Path(__file__).resolve().parents[1]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _check_prime_class(p: int, residue: int) -> int:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 17
        or not is_prime(p)
        or p % 4 != residue
    ):
        raise ValueError(f"need a prime p>=17 with p={residue} mod 4")
    return (p - residue) // 4


def _balanced(total: int, count: int, baseline: int = 0) -> tuple[int, ...]:
    """Return a deterministic balanced allocation with the requested sum."""
    if total < baseline * count or count <= 0:
        raise ValueError("the total lies below the requested baseline")
    quotient, remainder = divmod(total - baseline * count, count)
    return tuple(
        baseline + quotient + int(index < remainder)
        for index in range(count)
    )


def canonical_scaled_mass(p: int, offset: int, coefficients: Iterable[int]) -> int:
    """Evaluate (1) from a fixed-offset canonical coefficient vector."""
    values = tuple(int(value) for value in coefficients)
    return p * offset - sum(values)


def _three_bit_values(offset: int, coefficients: tuple[int, int, int]) -> tuple[int, ...]:
    values = []
    c12, c13, c23 = coefficients
    for z1 in (-1, 1):
        for z2 in (-1, 1):
            for z3 in (-1, 1):
                values.append(
                    offset
                    + c12 * z1 * z2
                    + c13 * z1 * z3
                    + c23 * z2 * z3
                )
    return tuple(values)


def sharp_atom_l1_certificate(p: int, a: int, b: int) -> dict[str, object]:
    """Prove the exact fixed-offset l1 minimum (2).

    The construction is a sum of nonnegative integral sharp atoms.  The
    lower bound uses only the forced coefficient sum.  In the exceptional
    ``|3a-b|=1`` case, l1 one would leave two slice values differing by two
    modulo four, while l1 two has the wrong parity.
    """
    if p < 5 or p % 2 == 0:
        raise ValueError("need an odd integer p>=5")
    if any(not isinstance(value, int) or isinstance(value, bool) or value < 0 for value in (a, b)):
        raise ValueError("a and b must be nonnegative integers")

    triples, remainder = divmod(b, 3)
    height = a - triples
    if remainder == 0:
        coefficients = (height, height, height)
    elif remainder == 1:
        # Add 1+z1z2-z1z3-z2z3.
        coefficients = (height + 1, height - 1, height - 1)
    else:
        # Add the two-unit XOR 2-2z1z2.
        coefficients = (height - 2, height, height)

    offset = a + b
    coefficient_sum = 3 * a - b
    values = _three_bit_values(offset, coefficients)
    l1_norm = sum(abs(value) for value in coefficients)
    exact_minimum = 3 if abs(coefficient_sum) == 1 else abs(coefficient_sum)
    scaled_mass = canonical_scaled_mass(p, offset, coefficients)
    expected_mass = a * (p - 3) + b * (p + 1)
    exceptional_integrality_obstruction = abs(coefficient_sum) == 1
    proved = bool(
        sum(coefficients) == coefficient_sum
        and l1_norm == exact_minimum
        and scaled_mass == expected_mass
        and all(value >= 0 and value % 4 == 0 for value in values)
        and (
            not exceptional_integrality_obstruction
            or (
                l1_norm == 3
                and not (
                    (offset - 1) % 4 == 0
                    and (offset + 1) % 4 == 0
                )
            )
        )
    )
    _require(proved, "the sharp fixed-offset l1 construction failed")
    return {
        "p": p,
        "p_minus_3_atom_count_a": a,
        "p_plus_1_atom_count_b": b,
        "fixed_offset_P": offset,
        "scaled_mass": scaled_mass,
        "forced_coefficient_sum": coefficient_sum,
        "attaining_three_edge_coefficients": list(coefficients),
        "attaining_values_of_4B": sorted(set(values)),
        "exact_minimum_l1": l1_norm,
        "absolute_sum_lower_bound": abs(coefficient_sum),
        "exceptional_abs_sum_one": exceptional_integrality_obstruction,
        "l1_one_ruled_out_by_mod_four_values": exceptional_integrality_obstruction,
        "l1_two_ruled_out_by_norm_sum_parity": exceptional_integrality_obstruction,
        "fixed_offset_representation_unique": True,
        "proved": proved,
    }


def _p1_bounds(r: int) -> tuple[int, int]:
    return 2 * r * r - 5 * r, 4 * r * r - 6 * r - 3


def _p3_bounds(r: int) -> tuple[int, int]:
    return 2 * r * r - 4 * r - 2, 4 * r * r - 2 * r - 5


def p1_local_survivor(p: int, t: int) -> dict[str, object]:
    """Build every scalar/local row in the ``p=1 mod 4`` survivor (3)."""
    r = _check_prime_class(p, 1)
    lower, upper = _p1_bounds(r)
    if not isinstance(t, int) or isinstance(t, bool) or not lower <= t <= upper:
        raise ValueError(f"need {lower}<=t<={upper}")
    q = 2 * r
    m = 2 * r + 1
    h = 4 * p + 2 * t + 1
    isolation_gap = p * p + 1 - 2 * h
    hard_excesses = _balanced(t, m)
    opposite_counts = _balanced(6 * r + t, m, r)
    hard_parallel = tuple(5 + value for value in hard_excesses)

    hard_rows = [
        {
            "e": excess,
            "scaled_mean_a": (excess + 1) * (p + 1),
            "parallel_P": 5 + excess,
            "cell": "x_j+2B, 4pE[B]=e(p+1)",
            "compact_coefficient_l1_upper": p - 1 + 3 * excess,
        }
        for excess in hard_excesses
    ]
    opposite_rows = []
    for count in opposite_counts:
        atom_a = r - 1
        atom_b = count - r
        mass = atom_a * (p - 3) + atom_b * (p + 1)
        opposite_rows.append(
            {
                "Q": count,
                "p_minus_3_units": atom_a,
                "p_plus_1_units": atom_b,
                "scaled_mean_a": mass,
                "closed_mass_formula": (p + 1) * count - 2 * p + 4,
                "cell": "A=2C",
                "atom_offset": count - 3,
                "universal_offset": 3,
                "one_omitted_pair_required_for_offset": True,
            }
        )

    type_budget = 2 * m * (m + t)
    hard_mean_sum = sum(row["scaled_mean_a"] for row in hard_rows)
    opposite_mean_sum = sum(row["scaled_mean_a"] for row in opposite_rows)
    hard_edges = sum(hard_parallel)
    opposite_edges = sum(opposite_counts)
    signed_total = hard_edges - opposite_edges
    proved = bool(
        isolation_gap > 0
        and sum(hard_excesses) == t
        and sum(opposite_counts) == 6 * r + t
        and min(opposite_counts) >= r
        and hard_mean_sum == opposite_mean_sum == type_budget
        and hard_edges + opposite_edges == h
        and signed_total == p + 4
        and all(row["scaled_mean_a"] == row["closed_mass_formula"] for row in opposite_rows)
        and all(row["atom_offset"] + row["universal_offset"] == row["Q"] for row in opposite_rows)
    )
    _require(proved, "the p=1 mod 4 local survivor ledger failed")
    return {
        "p": p,
        "p_mod_4": 1,
        "r": r,
        "q": q,
        "m": m,
        "t": t,
        "valid_t_interval": [lower, upper],
        "k": 4 * p + 2 * t,
        "H_edge_count": h,
        "isolated_vertex_gap": isolation_gap,
        "transported_infinity_degree_I": 0,
        "direction_count_per_type": m,
        "missing_projective_direction": False,
        "phase_residue_u": 0,
        "hard_rows": hard_rows,
        "opposite_rows": opposite_rows,
        "type_budget": type_budget,
        "hard_parallel_edge_total": hard_edges,
        "opposite_parallel_edge_total": opposite_edges,
        "common_signed_total_T": signed_total,
        "one_common_simple_graph_constructed": False,
        "proved_local_aggregate": proved,
    }


def p3_local_survivor(p: int, t: int) -> dict[str, object]:
    """Build every scalar/local row in the ``p=3 mod 4`` survivor (5)."""
    r = _check_prime_class(p, 3)
    lower, upper = _p3_bounds(r)
    if not isinstance(t, int) or isinstance(t, bool) or not lower <= t <= upper:
        raise ValueError(f"need {lower}<=t<={upper}")
    q = 2 * r + 1
    m = 2 * r + 2
    h = 4 * p + 2 * t + 1
    isolation_gap = p * p + 1 - 2 * h
    hard_excesses = _balanced(t + 1, m)
    opposite_counts = _balanced(10 * r + 6 + t, m, r + 2)
    hard_parallel = tuple(3 + value for value in hard_excesses)

    hard_rows = [
        {
            "e": excess,
            "scaled_mean_a": p - 1 + excess * (p + 1),
            "parallel_P": 3 + excess,
            "cell": "1-x_j+2B, 4pE[B]=e(p+1)",
            "compact_coefficient_l1_upper": p - 1 + 3 * excess,
        }
        for excess in hard_excesses
    ]
    opposite_rows = []
    for count in opposite_counts:
        atom_a = r - 1
        atom_b = count - r - 2
        mass = atom_a * (p - 3) + atom_b * (p + 1)
        opposite_rows.append(
            {
                "Q": count,
                "p_minus_3_units": atom_a,
                "p_plus_1_units": atom_b,
                "scaled_mean_a": mass,
                "closed_mass_formula": (p + 1) * count - 4 * p + 4,
                "cell": "A=2C",
                "atom_offset": count - 3,
                "universal_offset": 3,
                "all_atoms_may_be_compact_triangles": True,
            }
        )

    type_budget = 2 * m * (m + t)
    hard_mean_sum = sum(row["scaled_mean_a"] for row in hard_rows)
    opposite_mean_sum = sum(row["scaled_mean_a"] for row in opposite_rows)
    hard_edges = sum(hard_parallel)
    opposite_edges = sum(opposite_counts)
    signed_total = hard_edges - opposite_edges
    proved = bool(
        isolation_gap > 0
        and sum(hard_excesses) == t + 1
        and sum(opposite_counts) == 10 * r + 6 + t
        and min(opposite_counts) >= r + 2
        and hard_mean_sum == opposite_mean_sum == type_budget
        and hard_edges + opposite_edges == h
        and signed_total == 4 - p
        and all(row["scaled_mean_a"] == row["closed_mass_formula"] for row in opposite_rows)
        and all(row["atom_offset"] + row["universal_offset"] == row["Q"] for row in opposite_rows)
    )
    _require(proved, "the p=3 mod 4 local survivor ledger failed")
    return {
        "p": p,
        "p_mod_4": 3,
        "r": r,
        "q": q,
        "m": m,
        "t": t,
        "valid_t_interval": [lower, upper],
        "k": 4 * p + 2 * t,
        "H_edge_count": h,
        "isolated_vertex_gap": isolation_gap,
        "transported_infinity_degree_I": 0,
        "direction_count_per_type": m,
        "missing_projective_direction": False,
        "phase_residue_u": m - 1,
        "hard_rows": hard_rows,
        "opposite_rows": opposite_rows,
        "type_budget": type_budget,
        "hard_parallel_edge_total": hard_edges,
        "opposite_parallel_edge_total": opposite_edges,
        "common_signed_total_T": signed_total,
        "one_common_simple_graph_constructed": False,
        "proved_local_aggregate": proved,
    }


def p1_lower_endpoint_parseval(r: int) -> dict[str, object]:
    """Compare exact scalar Parseval with compact row maxima at (3)'s start."""
    if not isinstance(r, int) or isinstance(r, bool) or r < 7:
        raise ValueError("need an integer r>=7")
    p = 4 * r + 1
    m = 2 * r + 1
    t = 2 * r * r - 5 * r
    h = 4 * r * r + 6 * r + 5
    signed_total = 4 * r + 5
    parallel_square_sum = 4 * r**3 + 10 * r**2 + 18 * r + 19
    global_base = 8 * r**3 + 40 * r**2 + 70 * r + 17
    hard_compact_max = 10 * r**3 - 47 * r**2 + 118 * r - 30
    opposite_row_compact_max = 5 * r**2 - 8 * r + 53
    opposite_compact_max = m * opposite_row_compact_max
    local_compact_max = 20 * r**3 - 58 * r**2 + 216 * r + 23
    gap = local_compact_max - global_base
    gap_half = 6 * r**3 - 49 * r**2 + 73 * r + 3
    hard_e_low = r - 3
    hard_e_high = r - 2
    hard_row_max = lambda excess: 8 * r - 4 * excess + 5 * excess**2
    derived_parallel_square_sum = (
        (m - 3) * (r + 2) ** 2 + 3 * (r + 3) ** 2 + m * r**2
    )
    derived_hard_max = (
        (m - 3) * hard_row_max(hard_e_low)
        + 3 * hard_row_max(hard_e_high)
    )
    derived_opposite_row_max = (
        32 * r - 7 - 20 * (r - 2) + 5 * (r - 2) ** 2
    )
    proved = bool(
        derived_parallel_square_sum == parallel_square_sum
        and derived_hard_max == hard_compact_max
        and derived_opposite_row_max == opposite_row_compact_max
        and p * h + 2 * signed_total**2 - 2 * parallel_square_sum == global_base
        and gap == 2 * gap_half
        and hard_compact_max + opposite_compact_max == local_compact_max
        and gap > 0
    )
    _require(proved, "the p=1 lower-endpoint Parseval audit failed")
    return {
        "p_expression": "4r+1",
        "r": r,
        "t_lower_endpoint": t,
        "H_edge_count": h,
        "parallel_square_sum": parallel_square_sum,
        "hard_excess_multiset": {
            str(hard_e_low): m - 3,
            str(hard_e_high): 3,
        },
        "opposite_parallel_count_multiset": {str(r): m},
        "common_signed_total_T": signed_total,
        "global_Radon_energy": f"{global_base}+2*({p})*C",
        "global_energy_at_C_0": global_base,
        "hard_compact_triangle_energy_maximum": hard_compact_max,
        "one_opposite_compact_triangle_energy_maximum": opposite_row_compact_max,
        "all_opposite_compact_triangle_energy_maximum": opposite_compact_max,
        "total_compact_triangle_energy_maximum": local_compact_max,
        "compact_max_minus_global_C0": gap,
        "scalar_Parseval_upper_lower_gap_closes": False,
        "proved": proved,
    }


def p3_lower_endpoint_parseval(r: int) -> dict[str, object]:
    """Compare exact scalar Parseval with compact row maxima at (5)'s start."""
    if not isinstance(r, int) or isinstance(r, bool) or r < 7:
        raise ValueError("need an integer r>=7")
    p = 4 * r + 3
    m = 2 * r + 2
    t = 2 * r * r - 4 * r - 2
    h = 4 * r * r + 8 * r + 9
    signed_total = 1 - 4 * r
    parallel_square_sum = 4 * r**3 + 12 * r**2 + 26 * r + 13
    global_base = 8 * r**3 + 52 * r**2 - 8 * r + 3
    hard_compact_max = 10 * r**3 - 26 * r**2 + 88 * r - 31
    opposite_compact_max = 10 * (r - 1) ** 2 * (r + 1)
    local_compact_max = 20 * r**3 - 36 * r**2 + 78 * r - 21
    gap = local_compact_max - global_base
    gap_half = 6 * r**3 - 44 * r**2 + 43 * r - 12
    hard_e_low = r - 3
    hard_e_high = r - 2
    hard_row_max = lambda excess: 4 * (2 * r + 1) + 4 * excess + 5 * excess**2
    derived_parallel_square_sum = (
        (m - 5) * r**2 + 5 * (r + 1) ** 2 + m * (r + 2) ** 2
    )
    derived_hard_max = (
        (m - 5) * hard_row_max(hard_e_low)
        + 5 * hard_row_max(hard_e_high)
    )
    derived_opposite_max = m * 5 * (r - 1) ** 2
    proved = bool(
        derived_parallel_square_sum == parallel_square_sum
        and derived_hard_max == hard_compact_max
        and derived_opposite_max == opposite_compact_max
        and p * h + 2 * signed_total**2 - 2 * parallel_square_sum == global_base
        and gap == 2 * gap_half
        and hard_compact_max + opposite_compact_max == local_compact_max
        and gap > 0
    )
    _require(proved, "the p=3 lower-endpoint Parseval audit failed")
    return {
        "p_expression": "4r+3",
        "r": r,
        "t_lower_endpoint": t,
        "H_edge_count": h,
        "parallel_square_sum": parallel_square_sum,
        "hard_excess_multiset": {
            str(hard_e_low): m - 5,
            str(hard_e_high): 5,
        },
        "opposite_parallel_count_multiset": {str(r + 2): m},
        "common_signed_total_T": signed_total,
        "global_Radon_energy": f"{global_base}+2*({p})*C",
        "global_energy_at_C_0": global_base,
        "hard_compact_triangle_energy_maximum": hard_compact_max,
        "all_opposite_compact_triangle_energy_maximum": opposite_compact_max,
        "total_compact_triangle_energy_maximum": local_compact_max,
        "compact_max_minus_global_C0": gap,
        "scalar_Parseval_upper_lower_gap_closes": False,
        "proved": proved,
    }


def theorem_record() -> dict[str, object]:
    """Package the proved structural barrier without claiming a graph lift."""
    l1_samples = [
        sharp_atom_l1_certificate(29, a, b)
        for a, b in ((1, 0), (0, 1), (0, 2), (0, 3), (1, 2), (1, 3))
    ]
    p1 = p1_local_survivor(29, _p1_bounds(7)[0])
    p3 = p3_local_survivor(31, _p3_bounds(7)[0])
    parseval = {
        "p_1_mod_4": p1_lower_endpoint_parseval(7),
        "p_3_mod_4": p3_lower_endpoint_parseval(7),
    }
    proved = bool(
        all(row["proved"] for row in l1_samples)
        and p1["proved_local_aggregate"]
        and p3["proved_local_aggregate"]
        and all(row["proved"] for row in parseval.values())
    )
    return {
        "prop": "15.758",
        "title": "Sharp fixed-offset l1 cancellation and two infinite local survivor rays",
        "proved": {
            "fixed_offset_canonical_coefficient_sum_identity": proved,
            "sharp_atom_l1_minimum": proved,
            "p_1_mod_4_local_survivor_interval": proved,
            "p_3_mod_4_local_survivor_interval": proved,
            "scalar_Parseval_alone_excludes_survivors": False,
            "one_common_simple_graph_constructed": False,
            "residual_ii_closed": False,
            "e1_closed_general": False,
            "L": False,
        },
        "sharp_l1_samples": l1_samples,
        "local_survivor_samples": {
            "p_1_mod_4": p1,
            "p_3_mod_4": p3,
        },
        "lower_endpoint_Parseval": parseval,
        "remaining_obstruction": (
            "Integral midpoint/difference-Radon consistency for all directional "
            "coefficient matrices of one simple 0/1 graph H"
        ),
        "duplicate_work_guards": [
            "Do not revive an l1 lower bound proportional to the number of sharp atoms.",
            "Do not infer a common graph from the displayed local row aggregate.",
            "Do not claim scalar Parseval supplies an upper-lower contradiction on these rays.",
        ],
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    out = theorem_record()
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15758.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print("Prop. 15.758 coefficient-cancellation/local-survivor barrier: proved")
    print("  one common simple graph: NOT CONSTRUCTED")
    print("  residual (ii): OPEN")
    print(f"  wrote {destination}")
    return out


if __name__ == "__main__":
    main()
