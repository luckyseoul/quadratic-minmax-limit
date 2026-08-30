#!/usr/bin/env python3
"""Prop. 15.684 -- retracted low-mass reduction at p=23.

At the second all-finite endpoint ``p=23,s=20``, corrected quotient
arithmetic leaves phase-one residue ``u_1=11`` and phase-zero residues
``u_0=0,2,3,4,5,6,8,9``.  Proposition 15.681's universal scaled mass floor
removes ``u_0=2,3,4,5``; the valid low-value argument below removes
``u_0=6,8``.  The newly restored ``u_0=9`` row has scaled mass ``c=18`` and
is not covered by either argument.

If ``B`` is the resulting nonzero nonnegative integral quadratic on
``J(23,12)``, put ``c=4p E[B]`` and ``H=max B``.  The stabilizer identity
gives ``c>=4H``.  A paired 11-cube through a maximum has mean

    (H+c/4)/24.

The degree-two cube distance, its equality case, and a two-bit
Reed--Muller decomposition exclude ``c=12`` and the cases ``H<=3`` at
``c=16``.  If ``c=16,H=4``, equality in the stabilizer identity makes ``B``
vanish on the shell ``|X cap Y|=6``.  The kernel of restriction to that
shell has dimension 23 and equals ``(|X cap Y|-6)V_1``.  Affine
parallelograms at one and two replacements then make an integer value
congruent to ``-40/15`` modulo ``4/5``, an impossibility.

Conditional on residue zero, exact completion-bounded enumeration gives
1,247 phase-labelled profiles.  Segre's tangent envelope excludes all 363
arc profiles.  Coolsaet--Sticker's complete-arc classification of
``PG(2,23)`` and an off-conic secant count exclude every profile of pair
slack 4 and 8, all but one profile of slack 12, and all but one profile of
slack 16.  Exactly 203 profiles remain, all explicitly accounted for.

The residue-zero sublemmas remain valid, but the claim that only residue
zero remains is retracted.  The ``u_0=9`` row, the two low-slack exceptions,
and every residue-zero profile of slack at least 20 remain open.
"""
from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from functools import lru_cache

from e1_gmin_m4_prop15642 import stabilizer_mass_certificate
from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15681 import (
    exact_type_rows,
    paired_cube_integral_quadratic_floor,
    pair_slack_divisibility,
)
from e1_gmin_m4_prop15683 import tangent_envelope_input
from e1_gmin_m4_prop15723 import (
    backward_floor_plus_two_cell,
    floor_excess_admissible,
)


P = 23
M = 12
S = 20
PERIOD = 24
PAIR_DEFICIT_BUDGET = S * (S - 1)


def reed_muller_distance(degree: int) -> Fraction:
    """Elementary cube support density for a nonzero degree-d function."""
    if degree < 0:
        raise ValueError("degree must be nonnegative")
    return Fraction(1, 2**degree)


def low_value_cube_certificate() -> dict[str, object]:
    """The three cube facts needed at scaled masses 12 and 16."""
    degree_two = reed_muller_distance(2)
    degree_four = reed_muller_distance(4)
    value_three_floor = degree_two + 2 * degree_four
    if degree_two != Fraction(1, 4) or value_three_floor != Fraction(3, 8):
        raise ArithmeticError("cube distance arithmetic changed")
    return {
        "integer_mobius_coefficients": (
            "an integer-valued Boolean-cube function has integer multilinear "
            "Mobius coefficients"
        ),
        "nonzero_degree_two_mean_floor": degree_two,
        "degree_two_equality_case": (
            "mean 1/4 forces support density 1/4 and every nonzero value 1; "
            "therefore a value 2 is impossible"
        ),
        "value_three_decomposition": (
            "for f in {0,1,2,3}, f=(f mod 2)+2*(binom(f,2) mod 2)"
        ),
        "first_bit_degree": 2,
        "second_bit_degree": 4,
        "first_bit_support_floor": degree_two,
        "second_bit_support_floor": degree_four,
        "value_three_mean_floor": value_three_floor,
        "proof_method": "Reed--Muller distance by induction on cube dimension",
        "proved": True,
    }


def p23_small_mass_exclusion() -> dict[str, object]:
    """Exclude nonnegative integral quadratics of scaled mass 12 or 16."""
    stabilizer = stabilizer_mass_certificate(P)
    if (
        stabilizer["nodes"] != (5, 6, 12)
        or stabilizer["weights"]
        != (Fraction(0), Fraction(22, 23), Fraction(1, 23))
    ):
        raise ArithmeticError("p=23 stabilizer identity changed")

    cube = low_value_cube_certificate()
    case_rows = []
    for c in (12, 16):
        for height in range(1, c // 4 + 1):
            cube_mean = Fraction(height + c // 4, 24)
            if height == 1:
                lower = Fraction(1, 4)
                reason = "nonzero degree-two cube support"
            elif height == 2:
                lower = Fraction(1, 4)
                reason = (
                    "support floor, with equality impossible because the origin "
                    "has value two"
                )
            elif height == 3:
                lower = Fraction(3, 8)
                reason = "two-bit Reed--Muller decomposition at a value three"
            else:
                lower = None
                reason = "handled by the shell-kernel factor argument"
            excluded_by_cube = lower is not None and (
                cube_mean < lower or (height == 2 and cube_mean == lower)
            )
            case_rows.append(
                {
                    "scaled_mass": c,
                    "maximum_height": height,
                    "paired_cube_mean": cube_mean,
                    "cube_lower_bound": lower,
                    "reason": reason,
                    "excluded_by_cube": excluded_by_cube,
                }
            )

    # Restriction V_2(J(23,12)) -> V_2(J(12,6) x J(11,6)).
    v2_dimension = math.comb(23, 2)
    first_harmonics = (1, 11, 54)
    second_harmonics = (1, 10, 44)
    shell_image_dimension = sum(
        first_harmonics[i] * second_harmonics[j]
        for i in range(3)
        for j in range(3)
        if i + j <= 2
    )
    kernel_dimension = v2_dimension - shell_image_dimension
    v1_dimension = math.comb(23, 1)
    if (
        v2_dimension != 253
        or shell_image_dimension != 230
        or kernel_dimension != 23
        or v1_dimension != 23
    ):
        raise ArithmeticError("p=23 shell restriction dimensions changed")

    # If B=(t-6)L, B(X)=4 gives L(X)=2/3.  For one replacements
    # B(X-i+j)=b_ij and L=(b_ij)/5.  At two replacements integrality asks
    # 15 | 4(3(b_1+b_2)-10), impossible already modulo three.
    modular_obstruction = all(
        (4 * (3 * total - 10)) % 15 != 0 for total in range(0, 9)
    )
    if not modular_obstruction:
        raise ArithmeticError("two-replacement congruence acquired a solution")

    shell_case = {
        "scaled_mass": 16,
        "maximum_height": 4,
        "mean": Fraction(4, 23),
        "stabilizer_identity": "E[B]=(22/23)q(6)+B(X)/23",
        "forced_shell_average": 0,
        "forced_vanishing_shell": "B(Y)=0 whenever |X cap Y|=6",
        "domain_v2_dimension": v2_dimension,
        "shell_product": "J(12,6) x J(11,6)",
        "shell_filtered_harmonic_dimensions": {
            "J(12,6)": first_harmonics,
            "J(11,6)": second_harmonics,
        },
        "restriction_image_dimension": shell_image_dimension,
        "restriction_kernel_dimension": kernel_dimension,
        "factor_subspace_dimension": v1_dimension,
        "factorization": "B(Y)=(|X cap Y|-6)L(Y), with L affine",
        "factor_map_injective": (
            "an affine L vanishing on the t=5 and t=7 shells has equal "
            "coefficients in each block and then two roots as a linear "
            "function of t, hence L=0"
        ),
        "maximum_value_identity": "L(X)=2/3",
        "one_replacement_identity": "L(X-i+j)=B(X-i+j)/5",
        "two_replacement_identity": (
            "B(X-i1-i2+j1+j2)=4*(3*(b_i1j1+b_i2j2)-10)/15"
        ),
        "integrality_obstruction": (
            "15 cannot divide 4*(3n-10), since 3n-10 is nonzero modulo 3"
        ),
        "proved": modular_obstruction,
    }

    excluded_cube_cases = {
        (int(row["scaled_mass"]), int(row["maximum_height"]))
        for row in case_rows
        if row["excluded_by_cube"]
    }
    expected_cube_cases = {
        (12, 1),
        (12, 2),
        (12, 3),
        (16, 1),
        (16, 2),
        (16, 3),
    }
    if excluded_cube_cases != expected_cube_cases or not shell_case["proved"]:
        raise ArithmeticError("small-mass case split failed")
    return {
        "p": P,
        "middle_weight": M,
        "scaled_mass_definition": "c=4p*E[B]",
        "stabilizer_maximum_inequality": "c>=4H",
        "paired_cube_dimension": 11,
        "paired_cube_mean_at_maximum": "(H+c/4)/24",
        "cube_certificate": cube,
        "case_rows": case_rows,
        "height_four_shell_case": shell_case,
        "excluded_scaled_masses": [12, 16],
        "proved": True,
    }


def p23_u9_open_profile() -> dict[str, object]:
    """Verify an exact slack-zero profile in the restored ``u0=9`` row."""
    specifications = {
        "0": (0, 9, ((0, 0, 9), (20, 1, 3))),
        "1": (1, 11, ((2, 0, 11), (18, 1, 1))),
    }
    phases: dict[str, dict[str, object]] = {}
    all_admissible = True
    for key, (phase, u, entries) in specifications.items():
        quotient_sum = 0
        deficit = 0
        direction_count = 0
        rows = []
        for b, quotient, multiplicity in entries:
            floor_value = full_symbolic_floor(P, b, phase)
            excess = 2 * u + PERIOD * quotient - floor_value
            admissible = floor_excess_admissible(P, b, phase, excess)
            all_admissible = all_admissible and admissible
            rows.append(
                {
                    "b": b,
                    "quotient": quotient,
                    "multiplicity": multiplicity,
                    "floor": floor_value,
                    "excess": excess,
                    "admissible": admissible,
                }
            )
            quotient_sum += multiplicity * quotient
            deficit += multiplicity * (S - b)
            direction_count += multiplicity
        phases[key] = {
            "phase": phase,
            "u": u,
            "entries": rows,
            "direction_count": direction_count,
            "quotient_sum": quotient_sum,
            "deficit": deficit,
        }

    zero = phases["0"]
    one = phases["1"]
    total_deficit = int(zero["deficit"]) + int(one["deficit"])
    if (
        zero["direction_count"] != M
        or zero["quotient_sum"] != M - 9
        or zero["deficit"] != 180
        or one["direction_count"] != M
        or one["quotient_sum"] != M - 11
        or one["deficit"] != 200
        or not all_admissible
        or total_deficit != PAIR_DEFICIT_BUDGET
    ):
        raise ArithmeticError("p=23 u0=9 open profile changed")
    return {
        "u0": 9,
        "u1": 11,
        "scaled_mass_c": 18,
        "phases": phases,
        "total_deficit": total_deficit,
        "pair_slack": PAIR_DEFICIT_BUDGET - total_deficit,
        "restored_floor_plus_two_cell": backward_floor_plus_two_cell(P, 20, 0),
        "excluded": False,
        "proved_feasible_in_floor_relaxation": True,
    }


def p23_endpoint_residue_ledger() -> dict[str, object]:
    """Return the corrected exact type minima and positive-residue status."""
    phase_zero = exact_type_rows(P, 0)
    phase_one = exact_type_rows(P, 1)
    pair_rows = []
    for zero in phase_zero:
        for one in phase_one:
            required = int(zero["minimum_deficit"]) + int(
                one["minimum_deficit"]
            )
            if required <= PAIR_DEFICIT_BUDGET:
                pair_rows.append(
                    {
                        "u0": int(zero["u"]),
                        "u1": int(one["u"]),
                        "phase_zero_minimum_deficit": int(
                            zero["minimum_deficit"]
                        ),
                        "phase_one_minimum_deficit": int(
                            one["minimum_deficit"]
                        ),
                        "relaxed_pair_slack": PAIR_DEFICIT_BUDGET - required,
                        "phase_zero_profile": zero["profile"],
                        "phase_one_profile": one["profile"],
                    }
                )
    survivors = [int(row["u0"]) for row in pair_rows]
    if survivors != [0, 2, 3, 4, 5, 6, 8, 9]:
        raise ArithmeticError("p=23 pair-surviving residue list changed")
    if len(phase_one) != 1 or int(phase_one[0]["u"]) != 11:
        raise ArithmeticError("p=23 phase-one residue changed")

    universal = Fraction(
        paired_cube_integral_quadratic_floor(P)["universal_scaled_mass_floor"]
    )
    if universal != 12:
        raise ArithmeticError("p=23 universal scaled mass floor changed")
    special = p23_small_mass_exclusion()
    positive = []
    for row in pair_rows:
        u = int(row["u0"])
        if not u:
            continue
        scaled_mass = 2 * u
        if scaled_mass < universal:
            reason = "Proposition 15.681 universal integral quadratic floor"
        elif scaled_mass in special["excluded_scaled_masses"]:
            reason = "Proposition 15.684 low-value cube/shell exclusion"
        else:
            reason = "not excluded"
        positive.append(
            {
                "u0": u,
                "quotient_sum": M - u,
                "forces_quotient_zero": M - u < M,
                "least_positive_b_floor": min(
                    full_symbolic_floor(P, b, 0) for b in range(2, S + 1, 2)
                ),
                "therefore_b_zero": scaled_mass < 24,
                "factorization": "A=2B with B nonzero nonnegative integral quadratic",
                "scaled_mass_c": scaled_mass,
                "exclusion_reason": reason,
                "excluded": reason != "not excluded",
            }
        )
    positive_residues_all_excluded = all(
        bool(row["excluded"]) for row in positive
    )
    open_positive_residues = [
        int(row["u0"]) for row in positive if not bool(row["excluded"])
    ]
    if open_positive_residues != [9]:
        raise ArithmeticError("corrected p=23 open residue list changed")
    return {
        "p": P,
        "s": S,
        "pair_deficit_budget": PAIR_DEFICIT_BUDGET,
        "phase_zero_rows": phase_zero,
        "phase_one_rows": phase_one,
        "pair_survivors": pair_rows,
        "positive_residue_rows": positive,
        "positive_residues_all_excluded": positive_residues_all_excluded,
        "open_positive_residues": open_positive_residues,
        "open_positive_residue_witness": p23_u9_open_profile(),
        "residue_zero_remains": True,
        "small_mass_exclusion": special,
        "proved": True,
    }


@lru_cache(maxsize=None)
def _profile_rows(
    phase: int, u: int, deficit_cap: int
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Enumerate exact p=23 profiles with completion lower-bound pruning."""
    if phase not in (0, 1) or not 0 <= u < M:
        raise ValueError("phase must be 0/1 and 0<=u<12")
    target = M - u
    options: list[tuple[int, int, int]] = []
    for b in range(0, S + 1, 2):
        floor_value = full_symbolic_floor(P, b, phase)
        for quotient in range(target + 1):
            excess = 2 * u + PERIOD * quotient - floor_value
            if floor_excess_admissible(P, b, phase, excess):
                options.append((quotient, S - b, b))

    infinity = deficit_cap + S * M + 1
    completion = [[infinity] * (target + 1) for _ in range(M + 1)]
    completion[0][0] = 0
    for count in range(1, M + 1):
        for quotient_sum in range(target + 1):
            completion[count][quotient_sum] = min(
                (
                    added + completion[count - 1][quotient_sum - quotient]
                    for quotient, added, _b in options
                    if quotient <= quotient_sum
                ),
                default=infinity,
            )

    states: set[tuple[int, int, tuple[int, ...]]] = {(0, 0, ())}
    for count in range(M):
        next_states: set[tuple[int, int, tuple[int, ...]]] = set()
        for used, deficit, profile in states:
            for quotient, added, b in options:
                new_used = used + quotient
                new_deficit = deficit + added
                remaining_count = M - count - 1
                remaining_sum = target - new_used
                if (
                    new_used <= target
                    and new_deficit <= deficit_cap
                    and new_deficit
                    + completion[remaining_count][remaining_sum]
                    <= deficit_cap
                ):
                    next_states.add(
                        (new_used, new_deficit, tuple(sorted(profile + (b,))))
                    )
        states = next_states
    return tuple(
        sorted(
            (deficit, profile)
            for used, deficit, profile in states
            if used == target
        )
    )


@lru_cache(maxsize=1)
def _residue_zero_candidates() -> tuple[tuple[object, ...], ...]:
    phase_zero = _profile_rows(0, 0, PAIR_DEFICIT_BUDGET - 198)
    phase_one = _profile_rows(1, 11, PAIR_DEFICIT_BUDGET - 120)
    candidates = []
    for deficit_zero, profile_zero in phase_zero:
        for deficit_one, profile_one in phase_one:
            slack = PAIR_DEFICIT_BUDGET - deficit_zero - deficit_one
            if slack < 0 or slack % 4:
                continue
            secants = Counter((S - b) // 2 for b in profile_zero + profile_one)
            candidates.append(
                (
                    deficit_zero,
                    profile_zero,
                    deficit_one,
                    profile_one,
                    slack,
                    tuple(sorted(secants.items())),
                )
            )
    return tuple(sorted(candidates))


def _histogram(values: tuple[int, ...]) -> dict[str, int]:
    return {str(key): value for key, value in sorted(Counter(values).items())}


def p23_residue_zero_profile_census() -> dict[str, object]:
    """Aggregate the exact 1,247-profile residue-zero census."""
    phase_zero = _profile_rows(0, 0, PAIR_DEFICIT_BUDGET - 198)
    phase_one = _profile_rows(1, 11, PAIR_DEFICIT_BUDGET - 120)
    candidates = _residue_zero_candidates()
    slack_histogram = Counter(int(row[4]) for row in candidates)
    expected_slack = {
        0: 363,
        4: 264,
        8: 189,
        12: 136,
        16: 94,
        20: 68,
        24: 49,
        28: 35,
        32: 21,
        36: 13,
        40: 7,
        44: 4,
        48: 1,
        52: 1,
        56: 1,
        60: 1,
    }
    shapes_by_slack: dict[int, set[tuple[tuple[int, int], ...]]] = {}
    t0_by_slack: dict[int, Counter[int]] = {}
    for row in candidates:
        slack = int(row[4])
        secants = tuple(row[5])
        shapes_by_slack.setdefault(slack, set()).add(secants)
        t0_by_slack.setdefault(slack, Counter())[dict(secants).get(0, 0)] += 1
    expected_shape_counts = {
        0: 124,
        4: 95,
        8: 72,
        12: 54,
        16: 39,
        20: 30,
        24: 23,
        28: 17,
        32: 11,
        36: 8,
        40: 5,
        44: 3,
        48: 1,
        52: 1,
        56: 1,
        60: 1,
    }
    if (
        len(phase_zero) != 426
        or len(phase_one) != 11
        or len(candidates) != 1247
        or dict(slack_histogram) != expected_slack
        or {key: len(value) for key, value in shapes_by_slack.items()}
        != expected_shape_counts
    ):
        raise ArithmeticError("p=23 residue-zero profile census changed")

    canonical = [
        [row[0], list(row[1]), row[2], list(row[3]), row[4], list(row[5])]
        for row in candidates
    ]
    fingerprint = hashlib.sha256(
        json.dumps(canonical, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "p": P,
        "s": S,
        "pair_deficit_budget": PAIR_DEFICIT_BUDGET,
        "phase_zero_row_count": len(phase_zero),
        "phase_one_row_count": len(phase_one),
        "profile_count": len(candidates),
        "distinct_global_shape_count": sum(
            len(value) for value in shapes_by_slack.values()
        ),
        "pair_slack_histogram": dict(sorted(slack_histogram.items())),
        "shape_count_by_slack": {
            key: len(value) for key, value in sorted(shapes_by_slack.items())
        },
        "undetermined_direction_histogram_by_slack": {
            key: dict(sorted(value.items()))
            for key, value in sorted(t0_by_slack.items())
        },
        "canonical_profile_sha256": fingerprint,
        "pair_slack_divisibility": pair_slack_divisibility(),
        "proved": True,
    }


def p23_complete_arc_classification() -> dict[str, object]:
    """Record the exhaustive complete-arc classification used below."""
    return {
        "external_dependency": True,
        "source": (
            "K. Coolsaet and H. Sticker, A full classification of the "
            "complete k-arcs of PG(2,23) and PG(2,25), Journal of "
            "Combinatorial Designs 17 (2009), 459-477"
        ),
        "doi": "10.1002/jcd.20211",
        "independent_thesis_table": (
            "H. Sticker, Classification of Arcs in Small Desarguesian "
            "Projective Planes, Ghent PhD thesis (2012), Section 5.1"
        ),
        "complete_arc_sizes_pg2_23": [10, 12, 13, 14, 15, 16, 17, 24],
        "complete_arc_counts": {
            10: 1,
            12: 112449,
            13: 4341514,
            14: 1828196,
            15: 58361,
            16: 564,
            17: 5,
            24: 1,
        },
        "no_complete_arc_sizes": [18, 19, 20, 21, 22, 23],
        "unique_24_arc": "a nondegenerate conic",
        "finite_extension_consequence": (
            "every arc of size at least 18 extends to the unique-size class "
            "24 and is contained in a nondegenerate conic"
        ),
        "proved_conditional_on_external_classification": True,
    }


def p23_arc_envelope_exclusion() -> dict[str, object]:
    """Exclude all 363 slack-zero profiles."""
    signatures: Counter[tuple[int, int]] = Counter()
    for row in _residue_zero_candidates():
        if int(row[4]) != 0:
            continue
        secants = dict(row[5])
        high_count = sum(count for t, count in secants.items() if t <= 4)
        high_edges = sum(t * count for t, count in secants.items() if t <= 4)
        signatures[(high_count, high_edges)] += 1
    expected = {
        (3, 0): 43,
        (3, 1): 23,
        (3, 2): 11,
        (3, 3): 4,
        (4, 1): 3,
        (4, 2): 8,
        (4, 3): 24,
        (4, 4): 51,
        (4, 5): 67,
        (4, 6): 60,
        (4, 7): 20,
        (4, 8): 4,
        (5, 11): 34,
        (5, 12): 11,
    }
    if dict(signatures) != expected:
        raise ArithmeticError("p=23 arc envelope signatures changed")
    envelope = tangent_envelope_input()
    classification = p23_complete_arc_classification()
    return {
        "profile_count": sum(signatures.values()),
        "arc_size": 20,
        "tau": 5,
        "envelope_degree": 10,
        "high_direction_definition": "a direction with at most four secants",
        "minimum_high_direction_count": min(key[0] for key in signatures),
        "signature_histogram": [
            {"high_directions": key[0], "high_edges": key[1], "profiles": value}
            for key, value in sorted(signatures.items())
        ],
        "profiles_with_high_edge": sum(
            value for key, value in signatures.items() if key[1] > 0
        ),
        "high_edge_component_argument": (
            "d high direction squares leave degree r=10-2d.  At an endpoint "
            "of e>=1 high secants, the remaining 5-d+e tangents give "
            "r+2e>r zeros with multiplicity on its point-pencil, forcing that "
            "line.  After one point-pencil is removed, every other arc point "
            "has at least r surviving multiplicity on degree r-1, forcing "
            "more line components than the residual degree."
        ),
        "zero_high_edge_profiles": signatures[(3, 0)],
        "zero_high_edge_structure": (
            "exactly three high directions, all undetermined"
        ),
        "zero_high_edge_conic_argument": (
            "adjoin pairs from three undetermined infinity points.  Each "
            "22-arc is conic-contained by the complete-arc classification; "
            "the two conics share the original 20-arc and coincide, forcing "
            "one conic through three collinear infinity points."
        ),
        "tangent_envelope": envelope,
        "complete_arc_classification": classification,
        "excluded": True,
        "proved": True,
    }


def line_pair_slack(occupancy: int) -> int:
    """Contribution of one line to twice the pair-floor gap."""
    if occupancy < 0:
        raise ValueError("occupancy must be nonnegative")
    return 2 * (math.comb(occupancy, 2) - occupancy // 2)


def conic_core_repair_lemma() -> dict[str, object]:
    """Delete bad-line points, reach a conic, then count off-conic secants."""
    line_rows = [
        {
            "occupancy": n,
            "pair_slack": line_pair_slack(n),
            "deletions_to_occupancy_two": max(0, n - 2),
        }
        for n in range(3, 9)
    ]
    if not all(
        row["deletions_to_occupancy_two"] <= row["pair_slack"] // 4
        for row in line_rows
    ):
        raise ArithmeticError("bad-line deletion charge failed")
    off_conic = {
        h: {
            "off_conic_points": h,
            "minimum_secant_incidence_count": h * (7 - h),
            "pair_slack_floor": 4 * h * (7 - h),
        }
        for h in range(1, 5)
    }
    if min(row["pair_slack_floor"] for row in off_conic.values()) != 24:
        raise ArithmeticError("off-conic slack floor changed")
    return {
        "bad_line_charge_rows": line_rows,
        "arc_repair_bound": (
            "pair slack 4r permits deleting at most r points to obtain an arc; "
            "process every line of occupancy at least three and charge n-2 "
            "deletions to its line slack"
        ),
        "classification_threshold": 18,
        "conic_size": 24,
        "off_conic_full_secant_minimum": 11,
        "off_conic_count_rows": off_conic,
        "counting_argument": (
            "if h<=4 of the 20 original points are off a conic, each sees at "
            "least 11-(24-(20-h))=7-h secants whose two conic points remain. "
            "A line counted for r off-conic points has occupancy at least 2+r "
            "and contributes at least 4r pair slack."
        ),
        "consequence": (
            "once an arc repaired from a profile of slack below 24 is shown "
            "conic-contained, the original 20-point set is impossible unless "
            "it is itself an arc"
        ),
        "proved": True,
    }


def p23_low_slack_conic_exclusion() -> dict[str, object]:
    """Exclude slack 4,8 and almost all slack 12,16 profiles."""
    candidates = _residue_zero_candidates()
    rules = {
        4: {"delete_at_most": 1, "required_t0": 0, "adjoin": 0},
        8: {"delete_at_most": 2, "required_t0": 0, "adjoin": 0},
        12: {"delete_at_most": 3, "required_t0": 1, "adjoin": 1},
        16: {"delete_at_most": 4, "required_t0": 2, "adjoin": 2},
    }
    rows = []
    for slack, rule in rules.items():
        profiles = [row for row in candidates if int(row[4]) == slack]
        qualifying = [
            row
            for row in profiles
            if dict(row[5]).get(0, 0) >= int(rule["required_t0"])
        ]
        repaired_size_floor = S - int(rule["delete_at_most"])
        classified_arc_size_floor = repaired_size_floor + int(rule["adjoin"])
        if classified_arc_size_floor < 18:
            raise ArithmeticError("conic classification threshold not reached")
        rows.append(
            {
                "pair_slack": slack,
                "profile_count": len(profiles),
                "delete_at_most": rule["delete_at_most"],
                "repaired_arc_size_floor": repaired_size_floor,
                "required_undetermined_directions": rule["required_t0"],
                "adjoined_undetermined_infinity_points": rule["adjoin"],
                "classified_arc_size_floor": classified_arc_size_floor,
                "excluded_profile_count": len(qualifying),
                "remaining_profile_count": len(profiles) - len(qualifying),
            }
        )
    expected = {
        4: (264, 264, 0),
        8: (189, 189, 0),
        12: (136, 135, 1),
        16: (94, 93, 1),
    }
    observed = {
        int(row["pair_slack"]): (
            int(row["profile_count"]),
            int(row["excluded_profile_count"]),
            int(row["remaining_profile_count"]),
        )
        for row in rows
    }
    if observed != expected:
        raise ArithmeticError("p=23 conic-repair profile counts changed")
    return {
        "rows": rows,
        "repair_lemma": conic_core_repair_lemma(),
        "complete_arc_classification": p23_complete_arc_classification(),
        "excluded_profile_count": sum(
            int(row["excluded_profile_count"]) for row in rows
        ),
        "remaining_low_slack_profile_count": 2,
        "proved": True,
    }


def _exceptional_low_slack_profiles() -> list[dict[str, object]]:
    out = []
    for row in _residue_zero_candidates():
        deficit_zero, profile_zero, deficit_one, profile_one, slack, secants = row
        t0 = dict(secants).get(0, 0)
        if (slack == 12 and t0 == 0) or (slack == 16 and t0 < 2):
            out.append(
                {
                    "phase_deficits": {"0": deficit_zero, "1": deficit_one},
                    "phase_profiles_b": {
                        "0": _histogram(profile_zero),
                        "1": _histogram(profile_one),
                    },
                    "pair_slack": slack,
                    "global_floor_secant_distribution": {
                        str(key): value for key, value in secants
                    },
                    "undetermined_directions": t0,
                }
            )
    if len(out) != 2:
        raise ArithmeticError("exceptional low-slack profile count changed")
    return out


def p23_reduction_theorem() -> dict[str, object]:
    """Record the retracted reduction and its still-valid sublemmas."""
    residue = p23_endpoint_residue_ledger()
    census = p23_residue_zero_profile_census()
    arcs = p23_arc_envelope_exclusion()
    conic = p23_low_slack_conic_exclusion()
    excluded = int(arcs["profile_count"]) + int(conic["excluded_profile_count"])
    remaining_histogram = dict(census["pair_slack_histogram"])
    remaining_histogram[0] = 0
    remaining_histogram[4] = 0
    remaining_histogram[8] = 0
    remaining_histogram[12] = 1
    remaining_histogram[16] = 1
    remaining_histogram = {
        key: value for key, value in sorted(remaining_histogram.items()) if value
    }
    remaining = sum(remaining_histogram.values())
    if excluded != 1044 or remaining != 203 or excluded + remaining != 1247:
        raise ArithmeticError("p=23 reduction accounting failed")
    return {
        "prop": "15.684",
        "title": "Retracted low-mass and conic-core reduction at p=23",
        "record_status": "OPEN_RETRACTED_REDUCTION",
        "proved": False,
        "former_only_residue_zero_claim_retracted": True,
        "retraction_reason": (
            "the corrected floor-plus-two ledger restores u0=9 with scaled "
            "mass c=18, outside the c=12 and c=16 exclusions"
        ),
        "positive_residues_excluded": [2, 3, 4, 5, 6, 8],
        "open_positive_residues": [9],
        "only_residue_zero_remains": False,
        "residue_zero_reduction_proved_conditionally": True,
        "residue_zero_profile_count_before": 1247,
        "residue_zero_profiles_excluded": excluded,
        "residue_zero_profile_count_after": remaining,
        "remaining_pair_slack_histogram": remaining_histogram,
        "exceptional_low_slack_profiles": _exceptional_low_slack_profiles(),
        "p23_second_all_finite_endpoint_closed": False,
        "remaining_same_boundary_primes": [17, 19, 23],
        "top_level_gates_changed": False,
        "open_after_this_proposition": [
            "the restored p=23 phase-zero residue u0=9 (scaled mass c=18)",
            "the two exceptional p=23 profiles of slack 12 and 16",
            "the 201 p=23 profiles of pair slack at least 20",
            "the p=17 and p=19 second all-finite endpoints",
            "later all-finite sizes",
            "the infinity-present strict-deficit remainder",
            "residual (ii)",
            "R1",
            "global QVAR",
            "Type I",
            "the quadratic min-max limit",
        ],
        "residue_ledger": residue,
        "profile_census": census,
        "arc_exclusion": arcs,
        "low_slack_conic_exclusion": conic,
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return {
            "numerator": value.numerator,
            "denominator": value.denominator,
            "value": float(value),
        }
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def build_evidence() -> dict[str, object]:
    return _jsonable(p23_reduction_theorem())


def main() -> dict[str, object]:
    theorem = p23_reduction_theorem()
    print(
        "Prop. 15.684: OPEN (former only-residue-zero reduction retracted); "
        "u0=9 survives"
    )
    print("canonical evidence intentionally not regenerated")
    return theorem


if __name__ == "__main__":
    main()
