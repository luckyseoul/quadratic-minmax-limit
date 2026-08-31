#!/usr/bin/env python3
r"""Prop. 15.730 -- simultaneous endpoint repairs and their line census.

Assume the endpoint normal form of Proposition 15.727.  Thus ``D`` is an
affine ``(p+1)``-set, ``p=3R+c`` with ``c in {1,2}``, and its rich lines are
pairwise disjoint as subsets of ``D``.  There are ``x`` trisecants and ``y``
4-secants, with ``x+2y=R``.  Every rich block contributes two points to an
arc repair and one or two points to its complement.

The disjoint-block completion is stronger than the single construction in
Proposition 15.729.  Every choice of two points on each rich block, together
with every point outside the rich blocks, is a maximum arc in ``D``.  These
are all the maximum arcs in ``D``, so their exact number is

    3**x * 6**y.

For every such repair ``A`` its complement ``T=D\A`` is an ``R``-arc.  Each
point of ``T`` has exactly one ``A``-secant, namely its rich block.  The
selected secants form a matching on ``A`` and their fibres in ``T`` have
sizes one and two on the trisecant and 4-secant blocks, respectively.

Writing ``k=|A|=p+1-R`` and

    n_ij = #{lines ell : |ell intersect A|=i, |ell intersect T|=j},

the common completion forces the exact bivariate projective-line census

    n_00 = p(p-1)/2-R-y,       n_01 = 2R+2y,
    n_02 = binom(R,2)-y,
    n_10 = k+2R,               n_11 = R(k-2),       n_12 = 0,
    n_20 = binom(k,2)-R+y,     n_21 = x,            n_22 = y.

In particular, this is not merely a collection of global line moments.  It
also gives exact line-type signatures at every point.  For each complement
point ``z``, ``A union {z}`` is one of ``R`` simultaneous affine
unique-trisecant 3-arcs.  Deleting one endpoint of a selected secant exposes,
within ``D``, two individually valid co-tangent extensions on a trisecant
block and three on a 4-secant block.  Points outside ``D`` may give further
extensions and are not classified here.

This is a proved necessary normal form.  It does not exclude an endpoint or
assert that the resulting incidence data are realizable without the full
Proposition 15.727 hypotheses.
"""
from __future__ import annotations

import json
from math import comb
from pathlib import Path
from typing import Mapping

from e1_gmin_m4_prop15727 import endpoint_block_row, endpoint_residue_data


ROOT = Path(__file__).resolve().parents[1]
LINE_TYPE_KEYS = tuple(f"a{a}_t{t}" for a in range(3) for t in range(3))


def _line_key(a_points: int, t_points: int) -> str:
    return f"a{a_points}_t{t_points}"


def bivariate_line_census(p: int, four_secants: int) -> dict[str, int]:
    """Return the forced projective line counts for one endpoint block row."""
    block = endpoint_block_row(p, four_secants)
    R = int(block["R"])
    x = int(block["trisecants_x"])
    y = int(block["four_secants_y"])
    k = p + 1 - R
    return {
        "a0_t0": p * (p - 1) // 2 - R - y,
        "a0_t1": 2 * R + 2 * y,
        "a0_t2": comb(R, 2) - y,
        "a1_t0": k + 2 * R,
        "a1_t1": R * (k - 2),
        "a1_t2": 0,
        "a2_t0": comb(k, 2) - R + y,
        "a2_t1": x,
        "a2_t2": y,
    }


def verify_bivariate_line_census(
    p: int,
    four_secants: int,
    census: Mapping[str, int] | None = None,
) -> dict[str, object]:
    """Check every point, pair, cross-pair, and occupancy moment.

    ``census`` is injectable so focused tests can establish that the theorem
    predicate fails when an individual claimed count is corrupted.
    """
    block = endpoint_block_row(p, four_secants)
    R = int(block["R"])
    x = int(block["trisecants_x"])
    y = int(block["four_secants_y"])
    k = p + 1 - R
    counts = dict(bivariate_line_census(p, four_secants) if census is None else census)

    exact_keys = set(counts) == set(LINE_TYPE_KEYS)
    integral_nonnegative = bool(
        exact_keys
        and all(
            isinstance(counts[key], int)
            and not isinstance(counts[key], bool)
            and counts[key] >= 0
            for key in LINE_TYPE_KEYS
        )
    )
    if not integral_nonnegative:
        return {
            "exact_keys": exact_keys,
            "integral_nonnegative": False,
            "proved": False,
        }

    total_lines = sum(counts.values())
    a_incidence = sum(
        a * counts[_line_key(a, t)] for a in range(3) for t in range(3)
    )
    t_incidence = sum(
        t * counts[_line_key(a, t)] for a in range(3) for t in range(3)
    )
    a_pairs = sum(
        comb(a, 2) * counts[_line_key(a, t)]
        for a in range(3)
        for t in range(3)
    )
    t_pairs = sum(
        comb(t, 2) * counts[_line_key(a, t)]
        for a in range(3)
        for t in range(3)
    )
    cross_pairs = sum(
        a * t * counts[_line_key(a, t)]
        for a in range(3)
        for t in range(3)
    )
    a_tangents = sum(counts[_line_key(1, t)] for t in range(3))
    t_tangents = sum(counts[_line_key(a, 1)] for a in range(3))
    occupancy_counts = {
        occupancy: sum(
            counts[_line_key(a, t)]
            for a in range(3)
            for t in range(3)
            if a + t == occupancy
        )
        for occupancy in range(5)
    }
    expected_occupancy = {
        int(occupancy): int(count)
        for occupancy, count in dict(
            block["projective_line_occupancy_counts"]
        ).items()
    }
    checks = {
        "total_projective_lines": total_lines == p * p + p + 1,
        "A_point_line_incidences": a_incidence == k * (p + 1),
        "T_point_line_incidences": t_incidence == R * (p + 1),
        "A_pair_moment": a_pairs == comb(k, 2),
        "T_pair_moment": t_pairs == comb(R, 2),
        "A_T_cross_pair_moment": cross_pairs == k * R,
        "A_tangent_count": a_tangents == k * (R + 1),
        "T_tangent_count": t_tangents == R * (k + 1),
        "forbidden_a1_t2_absent": counts["a1_t2"] == 0,
        "trisecants_are_a2_t1": counts["a2_t1"] == x,
        "four_secants_are_a2_t2": counts["a2_t2"] == y,
        "D_occupancy_census": occupancy_counts == expected_occupancy,
    }
    return {
        "exact_keys": exact_keys,
        "integral_nonnegative": integral_nonnegative,
        "actual": {
            "total_projective_lines": total_lines,
            "A_point_line_incidences": a_incidence,
            "T_point_line_incidences": t_incidence,
            "A_pair_moment": a_pairs,
            "T_pair_moment": t_pairs,
            "A_T_cross_pair_moment": cross_pairs,
            "A_tangent_count": a_tangents,
            "T_tangent_count": t_tangents,
            "D_occupancy_census": occupancy_counts,
        },
        "expected": {
            "total_projective_lines": p * p + p + 1,
            "A_point_line_incidences": k * (p + 1),
            "T_point_line_incidences": R * (p + 1),
            "A_pair_moment": comb(k, 2),
            "T_pair_moment": comb(R, 2),
            "A_T_cross_pair_moment": k * R,
            "A_tangent_count": k * (R + 1),
            "T_tangent_count": R * (k + 1),
            "D_occupancy_census": expected_occupancy,
        },
        "checks": checks,
        "proved": all(checks.values()),
    }


def point_signature_rows(p: int, four_secants: int) -> list[dict[str, object]]:
    """Return every point type and its exact incident line-type signature."""
    block = endpoint_block_row(p, four_secants)
    R = int(block["R"])
    c = int(block["c"])
    x = int(block["trisecants_x"])
    y = int(block["four_secants_y"])
    k = p + 1 - R
    return [
        {
            "set": "A",
            "opposite_colour_points_on_rich_block": 0,
            "description": "point outside every rich block",
            "multiplicity": c + 1 + 2 * y,
            "incident_line_types": {
                "a2_t0": k - 1,
                "a1_t1": R,
                "a1_t0": 1,
            },
        },
        {
            "set": "A",
            "opposite_colour_points_on_rich_block": 1,
            "description": "A-point on a trisecant block",
            "multiplicity": 2 * x,
            "incident_line_types": {
                "a2_t1": 1,
                "a2_t0": k - 2,
                "a1_t1": R - 1,
                "a1_t0": 2,
            },
        },
        {
            "set": "A",
            "opposite_colour_points_on_rich_block": 2,
            "description": "A-point on a 4-secant block",
            "multiplicity": 2 * y,
            "incident_line_types": {
                "a2_t2": 1,
                "a2_t0": k - 2,
                "a1_t1": R - 2,
                "a1_t0": 3,
            },
        },
        {
            "set": "T",
            "opposite_colour_points_on_rich_block": 2,
            "description": "T-point on a trisecant block",
            "multiplicity": x,
            "incident_line_types": {
                "a2_t1": 1,
                "a1_t1": k - 2,
                "a0_t2": R - 1,
                "a0_t1": 2,
            },
        },
        {
            "set": "T",
            "opposite_colour_points_on_rich_block": 2,
            "description": "T-point on a 4-secant block",
            "multiplicity": 2 * y,
            "incident_line_types": {
                "a2_t2": 1,
                "a1_t1": k - 2,
                "a0_t2": R - 2,
                "a0_t1": 3,
            },
        },
    ]


def verify_point_signatures(p: int, four_secants: int) -> dict[str, object]:
    """Reaggregate the point signatures to the bivariate line census."""
    block = endpoint_block_row(p, four_secants)
    R = int(block["R"])
    k = p + 1 - R
    census = bivariate_line_census(p, four_secants)
    rows = point_signature_rows(p, four_secants)
    A_rows = [row for row in rows if row["set"] == "A"]
    T_rows = [row for row in rows if row["set"] == "T"]

    A_multiplicity = sum(int(row["multiplicity"]) for row in A_rows)
    T_multiplicity = sum(int(row["multiplicity"]) for row in T_rows)
    signatures_have_p_plus_one_lines = all(
        sum(dict(row["incident_line_types"]).values()) == p + 1
        for row in rows
        if int(row["multiplicity"]) > 0
    )

    A_reaggregation = {
        key: sum(
            int(row["multiplicity"])
            * int(dict(row["incident_line_types"]).get(key, 0))
            for row in A_rows
        )
        for key in LINE_TYPE_KEYS
    }
    T_reaggregation = {
        key: sum(
            int(row["multiplicity"])
            * int(dict(row["incident_line_types"]).get(key, 0))
            for row in T_rows
        )
        for key in LINE_TYPE_KEYS
    }
    expected_A = {key: int(key[1]) * census[key] for key in LINE_TYPE_KEYS}
    expected_T = {key: int(key[4]) * census[key] for key in LINE_TYPE_KEYS}
    checks = {
        "A_point_types_partition_A": A_multiplicity == k,
        "T_point_types_partition_T": T_multiplicity == R,
        "each_signature_has_p_plus_one_lines": signatures_have_p_plus_one_lines,
        "A_signatures_reaggregate": A_reaggregation == expected_A,
        "T_signatures_reaggregate": T_reaggregation == expected_T,
    }
    return {
        "rows": rows,
        "A_point_count": A_multiplicity,
        "T_point_count": T_multiplicity,
        "checks": checks,
        "proved": all(checks.values()),
    }


def common_completion_row(p: int, four_secants: int) -> dict[str, object]:
    """Package the simultaneous repair family forced by one endpoint row."""
    block = endpoint_block_row(p, four_secants)
    R = int(block["R"])
    c = int(block["c"])
    x = int(block["trisecants_x"])
    y = int(block["four_secants_y"])
    k = p + 1 - R
    rich_blocks = x + y
    points_on_rich_blocks = 3 * x + 4 * y
    singleton_points = p + 1 - points_on_rich_blocks
    maximum_repairs = 3**x * 6**y
    complement_points = x + 2 * y
    repair_points = singleton_points + 2 * rich_blocks
    line_census = bivariate_line_census(p, four_secants)
    census_audit = verify_bivariate_line_census(p, four_secants, line_census)
    signature_audit = verify_point_signatures(p, four_secants)

    two_extension_bases = 2 * x
    three_extension_bases = 2 * y
    one_extension_bases = singleton_points
    all_one_point_deletion_bases = (
        one_extension_bases + two_extension_bases + three_extension_bases
    )
    proved = bool(
        block["proved"]
        and block["rich_lines_pairwise_D_disjoint"]
        and x + 2 * y == R
        and singleton_points == c + 1 + 2 * y
        and repair_points == k
        and complement_points == R
        and maximum_repairs == int(block["maximum_arc_choice_count"])
        and all_one_point_deletion_bases == k
        and census_audit["proved"]
        and signature_audit["proved"]
    )
    return {
        "p": p,
        "R": R,
        "c": c,
        "repair_arc_size_k": k,
        "trisecants_x": x,
        "four_secants_y": y,
        "rich_block_count": rich_blocks,
        "points_on_rich_blocks": points_on_rich_blocks,
        "singleton_points": singleton_points,
        "repair_family": {
            "construction": (
                "retain all singleton points and choose exactly two points "
                "on every rich block"
            ),
            "maximum_D_subarc_count": maximum_repairs,
            "maximum_D_subarc_count_formula": "3^x 6^y",
            "these_are_all_maximum_D_subarcs": True,
            "reason": (
                "an arc takes at most two points from each rich block; "
                "block disjointness makes the displayed upper bound exact"
            ),
        },
        "complement_family": {
            "complement_size": complement_points,
            "complement_is_an_arc": True,
            "every_complement_point_has_A_secant_index": 1,
            "unique_secant_blocks_form_matching_on_A": True,
            "unique_secant_fibre_counts": {1: x, 2: y},
            "reason": (
                "a second A-secant through a complement point would be a "
                "second D-rich block through that point"
            ),
        },
        "simultaneous_unique_trisecants": {
            "count": R,
            "sets": "A union {z} for every z in T",
            "size": k + 1,
            "size_formula": "p+2-R",
            "all_affine": True,
            "each_has_exactly_one_trisecant": True,
        },
        "cotangent_deletion_bases": {
            "scope": "for each fixed maximum repair A; extension counts are within D",
            "base_size": k - 1,
            "base_size_formula": "p-R",
            "one_within_D_extension_base_count": one_extension_bases,
            "two_within_D_cotangent_extension_base_count": two_extension_bases,
            "three_within_D_cotangent_extension_base_count": three_extension_bases,
            "all_one_point_deletion_bases": all_one_point_deletion_bases,
            "displayed_D_extensions_individually_preserve_the_arc": True,
            "displayed_D_extensions_on_each_tangent_are_pairwise_incompatible": True,
            "three_displayed_D_extensions_occur_at_4_secant_A_endpoints": True,
            "additional_extension_points_outside_D_excluded": False,
        },
        "line_type_census": line_census,
        "line_type_census_audit": census_audit,
        "point_signature_audit": signature_audit,
        "finite_configuration_search_used": False,
        "endpoint_excluded": False,
        "top_level_gates_changed": False,
        "result_status": "proved simultaneous necessary normal form",
        "proved": proved,
    }


def direction_refinement_row(
    p: int,
    four_secants: int,
    *,
    A_secants: int,
    T_secants: int,
    rich_trisecants: int,
    rich_four_secants: int,
    ordinary_A_T_lines: int,
) -> dict[str, object]:
    """Audit a necessary bivariate census for one affine direction.

    This checks only the integer consequences of a proposed direction row;
    it does not claim that every accepted row is geometrically realizable.
    """
    block = endpoint_block_row(p, four_secants)
    R = int(block["R"])
    x = int(block["trisecants_x"])
    y = int(block["four_secants_y"])
    k = p + 1 - R
    sigma = A_secants
    tau = T_secants
    r3 = rich_trisecants
    r4 = rich_four_secants
    m = ordinary_A_T_lines
    parameters = (sigma, tau, r3, r4, m)
    integral_nonnegative_parameters = all(
        isinstance(value, int)
        and not isinstance(value, bool)
        and value >= 0
        for value in parameters
    )
    if not integral_nonnegative_parameters:
        raise ValueError("direction parameters must be nonnegative integers")

    counts = {
        "a0_t0": sigma + tau + r3 + r4 + m - 1,
        "a0_t1": R - 2 * tau - r3 - m,
        "a0_t2": tau - r4,
        "a1_t0": k - 2 * sigma - m,
        "a1_t1": m,
        "a1_t2": 0,
        "a2_t0": sigma - r3 - r4,
        "a2_t1": r3,
        "a2_t2": r4,
    }
    counts_nonnegative = all(value >= 0 for value in counts.values())
    b = p + 1 - 2 * (sigma + tau + m)
    slack = r3 + 2 * r4
    D_profile = {
        0: (p - 1 - b) // 2 + r3 + r4,
        1: b - r3,
        2: (p + 1 - b) // 2 - r3 - 2 * r4,
        3: r3,
        4: r4,
    }
    aggregated_profile = {
        occupancy: sum(
            count
            for key, count in counts.items()
            if int(key[1]) + int(key[4]) == occupancy
        )
        for occupancy in range(5)
    }
    checks = {
        "rich_counts_within_global_row": r3 <= x and r4 <= y,
        "line_types_nonnegative": counts_nonnegative,
        "p_affine_lines": sum(counts.values()) == p,
        "p_plus_one_D_points": sum(
            occupancy * count for occupancy, count in aggregated_profile.items()
        )
        == p + 1,
        "A_secant_count": sum(
            counts[_line_key(2, t)] for t in range(3)
        )
        == sigma,
        "T_secant_count": sum(
            counts[_line_key(a, 2)] for a in range(3)
        )
        == tau,
        "odd_fibre_count_even": b >= 0 and b % 2 == 0,
        "profile_formula": D_profile == aggregated_profile,
    }
    return {
        "p": p,
        "R": R,
        "four_secants_y": y,
        "line_scope": (
            "the p affine lines of one direction; the common empty line at "
            "infinity belongs only to the global projective a0_t0 count"
        ),
        "parameters": {
            "A_secants_sigma": sigma,
            "T_secants_tau": tau,
            "rich_trisecants_r3": r3,
            "rich_four_secants_r4": r4,
            "ordinary_A_T_lines_m": m,
        },
        "affine_line_type_census": counts,
        "odd_D_fibres_b": b,
        "direction_slack": slack,
        "D_fibre_profile": D_profile,
        "checks": checks,
        "geometric_realization_claimed": False,
        "proved": all(checks.values()),
    }


def proposition_15730() -> dict[str, object]:
    """Package sample rows while keeping the all-prime statement symbolic."""
    sample_primes = (31, 41)
    rows = {
        str(p): [
            common_completion_row(p, y)
            for y in range(int(endpoint_residue_data(p)["R"]) // 2 + 1)
        ]
        for p in sample_primes
    }
    proved = all(row["proved"] for prime_rows in rows.values() for row in prime_rows)
    return {
        "prop": "15.730",
        "title": "Simultaneous endpoint repair family and bivariate line census",
        "result_status": "proved simultaneous necessary normal form",
        "statement": (
            "every 15.727 disjoint-block completion has exactly 3^x 6^y "
            "maximum arc repairs, each with an R-arc index-one complement, "
            "the displayed exact 3-by-3 line census, and simultaneous "
            "two/three co-tangent extension bases"
        ),
        "source_correction": {
            "retracted_claim": (
                "Bartoli--Storme Corollary 2.7 is a universal odd-order "
                "unique-trisecant 3-arc size ceiling"
            ),
            "correct_reading": (
                "d<=2(q+2)/3+2 is the upper endpoint under the other "
                "hypotheses, including d>3+2sqrt(q) and existence of the "
                "configuration, for which the associated hyperplane "
                "arrangement is second-smallest; it is not a 3-arc "
                "nonexistence bound"
            ),
            "proposition_15_729_core_affected": False,
        },
        "sample_all_block_rows": rows,
        "finite_configuration_search_used": False,
        "endpoint_excluded": False,
        "p_plus_one_shell_closed": False,
        "non_walsh_residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "next_gate": (
            "exclude the common completion under residual direction, phase, "
            "and lift constraints; a nontrivial repair-cycle identity from "
            "Proposition 15.731 is a proposed attack, not a proved condition"
        ),
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic structural certificate."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15730.json"
    payload = json.dumps(proposition_15730(), indent=2, sort_keys=True) + "\n"
    output.write_text(payload)
    return output


def main() -> None:
    result = proposition_15730()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.730 audit failed")
    path = write_evidence()
    print("Prop 15.730 complementary-arc repair ensemble: proved")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
