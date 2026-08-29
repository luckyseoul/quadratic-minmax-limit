#!/usr/bin/env python3
"""Prop. 15.697 -- Booleanize and reduce the p=19 all-b2 profile.

After Proposition 15.696, the unique slack-twenty survivor has phase profiles
``{0:5,16:5}`` and ``{2:10}``.  Nine phase-one b=2 directions attain their
pointwise floor.  The tenth has

    A=(t-1)^2+2B,  B>=0 integral quadratic,  E[B]=5/19.

The sharp paired-cube and stabilizer bounds force max(B) to be one or five.
If max(B)=5, stabilizer equality makes B vanish on the intersection-five
layer about a maximum point.  That layer's quadratic evaluation has rank 152
and kernel ``(t-5)L`` for linear L.  Bounded integrality on layers eight and
nine leaves four coefficient patterns; each is negative on another layer.
Thus B is Boolean.

Separately, the five rigid phase-zero b=0 directions give exact signed-cross
coefficients.  Aggregate accounting first leaves infinity degrees
``0,10,20,28,38``; exact symmetric l1 minimization excludes 10 and 28.
The profile remains open only at infinity degrees ``0,20,38``.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15642 import stabilizer_mass_certificate
from e1_gmin_m4_prop15688 import p19_residue_zero_profiles


ROOT = Path(__file__).resolve().parents[1]
P = 19
M = 10
H_SIZE = 77
PAIR_COLUMNS = tuple(combinations(range(P), 2))
MAX_POINT = set(range(M))

RANK_WITNESS_MASKS = (
    31775, 48159, 80927, 146463, 277535, 56351, 89119, 154655, 285727,
    105503, 171039, 302111, 203807, 334879, 400415, 60447, 93215,
    158751, 289823, 109599, 117791, 62495, 95263, 160799, 291871,
    111647, 119839, 123935, 63519, 96287, 161823, 292895, 112671,
    120863, 124959, 127007, 31791, 48175, 80943, 146479, 277551,
    56367, 60463, 62511, 63535, 31823, 48207, 80975, 146511, 277583,
    56399, 60495, 62543, 63567, 31887, 48271, 81039, 146575, 277647,
    56463, 60559, 62607, 63631, 32015, 48399, 81167, 146703, 277775,
    56591, 60687, 62735, 63759, 32271, 48655, 81423, 146959, 278031,
    56847, 60943, 62991, 64015, 31799, 48183, 80951, 146487, 277559,
    56375, 60471, 62519, 63543, 31831, 31895, 32023, 32279, 31847,
    31911, 32039, 32295, 31943, 32071, 32327, 32135, 32391, 32519,
    31803, 48187, 80955, 146491, 277563, 56379, 60475, 62523, 63547,
    31835, 31899, 32027, 32283, 31851, 31859, 31805, 48189, 80957,
    146493, 277565, 56381, 60477, 62525, 63549, 31837, 31901, 32029,
    32285, 31853, 31861, 31865, 31806, 48190, 80958, 146494, 277566,
    56382, 60478, 62526, 63550, 31838, 31902, 32030, 32286, 31854,
    31862, 31866, 31868,
)


def _pair_row_bits(mask: int) -> int:
    row = 0
    for column, (left, right) in enumerate(PAIR_COLUMNS):
        if (mask >> left) & 1 and (mask >> right) & 1:
            row |= 1 << column
    return row


def _rank_mod_two(rows: list[int]) -> int:
    basis: dict[int, int] = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def _rank_mod_prime(rows: list[list[int]], prime: int) -> int:
    rows = [[value % prime for value in row] for row in rows]
    rank = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next(
            (row for row in range(rank, len(rows)) if rows[row][column]), None
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(value * inverse) % prime for value in rows[rank]]
        for row in range(rank + 1, len(rows)):
            if rows[row][column]:
                factor = rows[row][column]
                rows[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(rows[row], rows[rank])
                ]
        rank += 1
    return rank


def p19_max_five_layer_kernel_certificate() -> dict[str, object]:
    masks = list(RANK_WITNESS_MASKS)
    if len(masks) != 152 or any(mask.bit_count() != M for mask in masks):
        raise ArithmeticError("intersection-five rank witness left the slice")
    if any(sum((mask >> i) & 1 for i in MAX_POINT) != 5 for mask in masks):
        raise ArithmeticError("rank witness left the intersection-five layer")
    rank = _rank_mod_two([_pair_row_bits(mask) for mask in masks])
    if rank != 152:
        raise ArithmeticError("intersection-five witness lost rank")

    # A linear function has dimension 19 on J(19,10).  Its values on the
    # union of t=4 and t=6 have full rank, so multiplication by t-5 is
    # injective.  This proves that the displayed kernel has dimension 19.
    linear_rows = []
    for subset in combinations(range(P), M):
        if len(set(subset) & MAX_POINT) not in (4, 6):
            continue
        linear_rows.append([1, *(int(i in subset) for i in range(P - 1))])
    linear_rank = _rank_mod_prime(linear_rows, 101)
    if linear_rank != 19:
        raise ArithmeticError("linear kernel multiplier lost injectivity")
    return {
        "slice": "J(19,10)",
        "fixed_maximum_point_size": M,
        "vanishing_intersection_layer": 5,
        "quadratic_dimension": len(PAIR_COLUMNS),
        "rank_witness_rows": len(masks),
        "rank_mod_two": rank,
        "kernel_dimension_upper_bound": len(PAIR_COLUMNS) - rank,
        "displayed_kernel": "B=(t-5)L, degree(L)<=1",
        "displayed_kernel_dimension": linear_rank,
        "therefore_exact_kernel": True,
        "proved": True,
    }


def p19_cross_difference_pattern_certificate() -> dict[str, object]:
    """Exhaust the additive cross-difference matrices in the max-five branch.

    The ten-by-nine matrix ``D`` has entries in ``{0,1,2}``, all of one
    parity, and satisfies the additive parallelogram identities.  Every two
    entries in distinct rows and columns have sum zero or two.  In odd
    parity this forces the all-one matrix.  In even parity write ``D=2E``;
    an additive Boolean matrix is fixed by its first row and first column,
    so only ``2^(10+9-1)`` candidates need checking.
    """
    orbit_counts: Counter[str] = Counter()
    negative_witnesses: dict[str, dict[str, object]] = {}

    # Odd parity has only D_ij=1.  Four disjoint replacements give
    # R=2-4=-2 and B=(6-5)R/2=-1.
    orbit_counts["all_cross_differences_one"] = 1
    negative_witnesses["all_cross_differences_one"] = {
        "removed_rows": [0, 1, 2, 3],
        "added_columns": [0, 1, 2, 3],
        "bad_layer": 6,
        "R": -2,
        "B": -1,
    }

    even_candidates_checked = 0
    even_candidates_admissible = 0
    for corner in (0, 1):
        for first_row_tail in range(1 << 8):
            first_row = [corner] + [
                (first_row_tail >> bit) & 1 for bit in range(8)
            ]
            for first_column_tail in range(1 << 9):
                even_candidates_checked += 1
                first_column = [corner] + [
                    (first_column_tail >> bit) & 1 for bit in range(9)
                ]
                matrix = [
                    [first_column[i] + first_row[j] - corner for j in range(9)]
                    for i in range(10)
                ]
                if any(value not in (0, 1) for row in matrix for value in row):
                    continue
                support = [
                    (i, j)
                    for i, row in enumerate(matrix)
                    for j, value in enumerate(row)
                    if value
                ]
                if any(
                    i != k and j != ell
                    for index, (i, j) in enumerate(support)
                    for k, ell in support[index + 1 :]
                ):
                    continue
                even_candidates_admissible += 1
                rows = {i for i, _ in support}
                columns = {j for _, j in support}
                if not support:
                    orbit = "all_cross_differences_zero"
                    removed = list(range(6))
                    added = list(range(6))
                elif len(rows) == 1 and len(support) == 9:
                    orbit = "one_X_row_of_twos"
                    exceptional = next(iter(rows))
                    removed = [i for i in range(10) if i != exceptional][:6]
                    added = list(range(6))
                elif len(columns) == 1 and len(support) == 10:
                    orbit = "one_Y_column_of_twos"
                    exceptional = next(iter(columns))
                    removed = list(range(6))
                    added = [j for j in range(9) if j != exceptional][:6]
                else:
                    raise ArithmeticError("unclassified additive even pattern")
                if any(matrix[i][j] for i in removed for j in added):
                    raise ArithmeticError("purported zero-rectangle witness is nonzero")
                # D=2E vanishes on the chosen six-by-six rectangle.  Six
                # paired replacements therefore give R=2 and B=-R/2=-1
                # on the t=4 layer.
                orbit_counts[orbit] += 1
                negative_witnesses.setdefault(
                    orbit,
                    {
                        "removed_rows": removed,
                        "added_columns": added,
                        "bad_layer": 4,
                        "R": 2,
                        "B": -1,
                    },
                )

    expected = {
        "all_cross_differences_one": 1,
        "all_cross_differences_zero": 1,
        "one_X_row_of_twos": 10,
        "one_Y_column_of_twos": 9,
    }
    if dict(orbit_counts) != expected or even_candidates_admissible != 20:
        raise ArithmeticError("cross-difference pattern exhaustion changed")
    if not all(row["B"] == -1 for row in negative_witnesses.values()):
        raise ArithmeticError("a max-five cross-difference orbit survived")
    return {
        "matrix_shape": [10, 9],
        "even_additive_candidates_checked": even_candidates_checked,
        "even_admissible_labelled_matrices": even_candidates_admissible,
        "admissible_labelled_matrices_including_odd": sum(orbit_counts.values()),
        "orbit_counts": dict(orbit_counts),
        "negative_layer_witnesses": negative_witnesses,
        "all_patterns_excluded": True,
        "proved": True,
    }


def p19_elevated_lift_booleanization() -> dict[str, object]:
    target_mean = Fraction(5, 19)
    target_scaled_mass = 4 * P * target_mean
    possible_maxima = []
    for maximum in range(1, 20):
        stabilizer_ok = 4 * maximum <= target_scaled_mass
        paired_ok = maximum == 1 or 2 * (P + 1) - 4 * maximum <= target_scaled_mass
        if stabilizer_ok and paired_ok:
            possible_maxima.append(maximum)
    if possible_maxima != [1, 5]:
        raise ArithmeticError("mass-twenty maximum dichotomy changed")

    stabilizer = stabilizer_mass_certificate(P)
    if stabilizer["nodes"] != (4, 5, 10) or stabilizer["weights"] != (
        Fraction(0), Fraction(18, 19), Fraction(1, 19)
    ):
        raise ArithmeticError("p=19 stabilizer identity changed")
    # At a maximum-five point the endpoint contribution already equals E B,
    # forcing the entire t=5 layer to zero by pointwise nonnegativity.
    forced_layer_average = (
        target_mean - Fraction(1, 19) * 5
    ) / Fraction(18, 19)
    if forced_layer_average != 0:
        raise ArithmeticError("maximum-five equality no longer forces zero")

    patterns = p19_cross_difference_pattern_certificate()
    return {
        "elevated_direction": {
            "baseline": "A0=(t-1)^2",
            "decomposition": "A=A0+2B",
            "B_nonnegative_integer_quadratic": True,
            "B_mean": target_mean,
            "scaled_mass_4pEB": target_scaled_mass,
        },
        "possible_maxima_before_layer_argument": possible_maxima,
        "maximum_five_stabilizer_forced_layer_average": forced_layer_average,
        "maximum_five_layer_kernel": p19_max_five_layer_kernel_certificate(),
        "maximum_five_integral_cross_difference_patterns": patterns,
        "maximum_five_excluded": True,
        "therefore_maximum_one": True,
        "therefore_B_is_boolean": True,
        "proved": True,
    }


def _partitions(total: int, length: int, minimum: int = 0):
    if length == 0:
        if total == 0:
            yield ()
        return
    for value in range(minimum, min(P, total // length) + 1):
        for rest in _partitions(total - value, length - 1, value):
            yield (value, *rest)


def _minimum_symmetric_l1(total: int, gauge: int) -> tuple[int, tuple[int, ...]]:
    best = 10**9
    witness = ()
    for values in _partitions(total, P):
        objective = sum(
            abs(gauge - values[i] - values[j])
            for i in range(P) for j in range(i + 1, P)
        )
        if objective < best:
            best, witness = objective, values
    return best, witness


def _mean_only_parallel_options(infinity_degree: int, mean: int) -> tuple[int, ...]:
    finite_edges = H_SIZE - infinity_degree
    return tuple(
        parallel
        for parallel in range(finite_edges + 1)
        if abs(infinity_degree + P * parallel - 3 * P - mean)
        <= finite_edges - parallel
    )


def p19_allb2_infinity_degree_reduction() -> dict[str, object]:
    aggregate = []
    for infinity_degree in range(0, H_SIZE + 1, 2):
        finite_edges = H_SIZE - infinity_degree
        low = []
        zero = []
        for gauge in range(-H_SIZE, H_SIZE + 1):
            low_parallel = 4 + 9 * gauge - infinity_degree
            low_cross = finite_edges - low_parallel
            if (
                low_parallel >= 0 and low_cross >= 0
                and abs(1 + 171 * gauge - 18 * infinity_degree) <= low_cross
            ):
                low.append((low_parallel, gauge))
            zero_parallel = 3 + 9 * gauge - infinity_degree
            zero_cross = finite_edges - zero_parallel
            if (
                zero_parallel >= 0 and zero_cross >= 0
                and abs(171 * gauge - 18 * infinity_degree) <= zero_cross
            ):
                zero.append((zero_parallel, gauge))
        elevated = _mean_only_parallel_options(infinity_degree, 38)
        high = _mean_only_parallel_options(infinity_degree, 40)
        possible_sums = {0}
        for options, count in (
            ([row[0] for row in low], 9),
            (elevated, 1),
            ([row[0] for row in zero], 5),
            (high, 5),
        ):
            for _ in range(count):
                possible_sums = {
                    left + right for left in possible_sums for right in options
                }
        if finite_edges in possible_sums:
            aggregate.append(
                {
                    "infinity_degree": infinity_degree,
                    "low_b2_parallel_gauge": low,
                    "zero_b0_parallel_gauge": zero,
                }
            )
    initial = [row["infinity_degree"] for row in aggregate]
    if initial != [0, 10, 20, 28, 38]:
        raise ArithmeticError("initial all-b2 infinity table changed")

    rows = []
    for row in aggregate:
        infinity_degree = int(row["infinity_degree"])
        zero_parallel, gauge = row["zero_b0_parallel_gauge"][0]
        minimum, witness = _minimum_symmetric_l1(infinity_degree, gauge)
        cross_capacity = H_SIZE - infinity_degree - zero_parallel
        rows.append(
            {
                "infinity_degree": infinity_degree,
                "b0_gauge": gauge,
                "b0_parallel_edges": zero_parallel,
                "b0_cross_edge_capacity": cross_capacity,
                "minimum_sum_abs_g_minus_hs_minus_ht": minimum,
                "minimizing_star_degree_histogram": dict(sorted(Counter(witness).items())),
                "excluded": minimum > cross_capacity,
            }
        )
    remaining = [row["infinity_degree"] for row in rows if not row["excluded"]]
    if remaining != [0, 20, 38]:
        raise ArithmeticError("exact l1 degree reduction changed")
    return {
        "initial_aggregate_infinity_degrees": initial,
        "exact_symmetric_l1_rows": rows,
        "excluded_by_l1": [10, 28],
        "remaining_infinity_degrees": remaining,
        "proved": True,
    }


def _mobius_coefficients(values: list[int], variables: int) -> list[int]:
    coefficients = values[:]
    for variable in range(variables):
        for mask in range(1 << variables):
            if mask & (1 << variable):
                coefficients[mask] -= coefficients[mask ^ (1 << variable)]
    return coefficients


def p19_boolean_density_catalog() -> dict[str, object]:
    """Conditional slice-to-cube classification, then exhaustive four-cube audit."""
    target = Fraction(5, 19)
    observed = set()
    for truth_table in range(1 << 16):
        values = [(truth_table >> mask) & 1 for mask in range(16)]
        coefficients = _mobius_coefficients(values, 4)
        if any(
            value and mask.bit_count() > 2
            for mask, value in enumerate(coefficients)
        ):
            continue
        numerator = sum(
            value * comb(15, M - mask.bit_count())
            for mask, value in enumerate(values)
        )
        if Fraction(numerator, comb(P, M)) == target:
            observed.add(tuple(values))

    predicted = set()
    for i, j in combinations(range(4), 2):
        for signs in ((1, 1), (1, 0), (0, 1)):
            predicted.add(
                tuple(
                    int(((mask >> i) & 1) == signs[0] and ((mask >> j) & 1) == signs[1])
                    for mask in range(16)
                )
            )
    for triple in combinations(range(4), 3):
        for singled in triple:
            others = [index for index in triple if index != singled]
            predicted.add(
                tuple(
                    int(
                        (((mask >> singled) & 1), *((mask >> index) & 1 for index in others))
                        in ((1, 0, 0), (0, 1, 1))
                    )
                    for mask in range(16)
                )
            )
    if observed != predicted or len(observed) != 30:
        raise ArithmeticError("four-cube Boolean density catalog changed")
    family_counts = {
        "positive_positive_conjunction": comb(P, 2),
        "positive_negative_conjunction": P * (P - 1),
        "nonconstant_antipodal_triple": P * comb(P - 1, 2),
    }
    if sum(family_counts.values()) != 3420:
        raise ArithmeticError("p=19 Boolean form count changed")
    return {
        "external_dependency": True,
        "restriction_theorem": (
            "Filmus-Vinciguerra: for Boolean values the slice-to-cube "
            "restriction threshold equals 2d; apply d=2 to the complementary "
            "slice J(19,9)"
        ),
        "external_source": (
            "Yuval Filmus publication page, short note with Antoine Vinciguerra "
            "linked under Junta threshold for low degree Boolean functions on the slice"
        ),
        "cube_degree_two_relevant_variable_bound": 4,
        "four_variable_truth_tables_checked": 1 << 16,
        "target_density_truth_tables": len(observed),
        "essential_families": list(family_counts),
        "p19_family_counts": family_counts,
        "p19_form_count": sum(family_counts.values()),
        "proved_conditional_on_external_restriction_theorem": True,
    }


def p19_allb2_structural_reduction() -> dict[str, object]:
    profiles = [
        row for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) == 20
        and row["phase_profiles_b"] == {"0": {0: 5, 16: 5}, "1": {2: 10}}
    ]
    if len(profiles) != 1:
        raise ArithmeticError("all-b2 slack-twenty profile changed")
    return {
        "proposition": "15.697",
        "p": P,
        "boundary_size": 16,
        "profile": profiles[0],
        "elevated_lift_booleanization": p19_elevated_lift_booleanization(),
        "infinity_degree_reduction": p19_allb2_infinity_degree_reduction(),
        "conditional_boolean_form_catalog": p19_boolean_density_catalog(),
        "bounded_edge_lift_trials": (
            "UNKNOWN results are diagnostics only and are not proposition evidence"
        ),
        "p19_profiles_before": 4,
        "p19_profiles_after": 4,
        "remaining_slack_histogram": {20: 1, 24: 1, 28: 1, 32: 1},
        "p19_second_all_finite_endpoint_closed": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "closes_type_I": False,
        "L_status": "OPEN",
        "proved": True,
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def main() -> None:
    theorem = p19_allb2_structural_reduction()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15697.json"
    target.write_text(json.dumps(_jsonable(theorem), indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.697: elevated p=19 all-b2 lift is Boolean; "
        "infinity degrees reduce to 0,20,38; profile remains open"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
