#!/usr/bin/env python3
r"""Prop. 15.748 -- exact interpolation in the p13 ``P=5`` branch.

Proposition 15.747 makes every minimum ``Q=3`` opposite cell in the
all-equal-triple branch an exact literal.  If ``z`` is the number of those
minimum cells, their directions are common roots of the homogeneous forms
``M_2,M_4,M_6``.  The opposite excesses are nonnegative integers with sum
five, so ``z>=2``.

Five literal roots would make ``M_4`` identically zero, contradicting the
nonzero fourth moment of every hard baseline-pair/all-equal-triple cell.
The exact interpolation runner checks the remaining ``z=4`` and ``z=3``
forms against the complete 69-element hard moment alphabet and finds none,
for either hard sign.  At ``z=2`` it finds exactly 336 moment-level
survivors for each sign.  Therefore only excess partition ``1^5`` remains;
this is a proved exact reduction, not a close of the branch or residue.
"""
from __future__ import annotations

import hashlib
import json
import sys
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.p13_p5_literal_interpolation import (  # noqa: E402
    FOURTH_ALPHABET,
    P,
    POINTS,
    audit,
    hard_moment_alphabet,
    quartic_code_contains,
    quartic_code_membership_setup,
    root_product_values,
)

from e1_gmin_m4_prop15747 import proposition_15747  # noqa: E402
from e1_gmin_m4_prop15746 import t4_u4_catalog_consequence  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402


HARD_DIRECTION_COUNT = 7
OPPOSITE_DIRECTION_COUNT = 7
OPPOSITE_EXCESS_SUM = 5
EXPECTED_PAYLOAD_SHA256 = (
    "894c087d4acae7ff0722ba236b1fac494984b9b331431e6117b2edbde0afbbec"
)
EXPECTED_ALPHABET_SHA256 = (
    "5088c5d586f8b651b25c5b2c15df6b940a00f4631c5aed4f4669adac1a0a9b25"
)
EXPECTED_Z2_SHA256 = {
    -1: "71c92ddb3fab7fe6319665ad1b1c2343dd45b791b2ba3ed98c52ae05198c0de3",
    1: "97706d838541e87a2d4dacfa8188627b8160e59fbd6aa83167db0e0b91dab585",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _mod_rank(rows: list[list[int]]) -> int:
    work = [[value % P for value in row] for row in rows]
    rank = 0
    if not work:
        return rank
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, P)
        work[rank] = [value * inverse % P for value in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            scalar = work[row][column]
            work[row] = [
                (left - scalar * right) % P
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == len(work[0]):
            break
    return rank


def _catalog_digest(rows: object) -> str:
    payload = json.dumps(rows, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@lru_cache(maxsize=1)
def literal_root_and_hard_alphabet_certificate() -> dict[str, object]:
    """Check literal moments, the hard alphabet, and the five-root rule."""
    degrees = (2, 4, 6)
    literal_star_sums = {
        degree: {
            center: sum(
                pow(center - other, degree, P)
                for other in range(P)
                if other != center
            )
            % P
            for center in range(P)
        }
        for degree in degrees
    }
    every_literal_moment_zero = all(
        value == 0
        for by_center in literal_star_sums.values()
        for value in by_center.values()
    )

    alphabet = tuple(sorted(hard_moment_alphabet()))
    alphabet_digest = _catalog_digest(alphabet)
    fourth_values = sorted({n4 for _n2, n4, _n6 in alphabet})
    fourth_given_n2_zero = sorted(
        {n4 for n2, n4, _n6 in alphabet if n2 == 0}
    )

    five_root_ranks = []
    for roots in combinations(POINTS, 5):
        evaluation = [
            [
                pow(x, 4 - monomial, P) * pow(y, monomial, P) % P
                for monomial in range(5)
            ]
            for x, y in roots
        ]
        five_root_ranks.append(_mod_rank(evaluation))

    proved = bool(
        len(POINTS) == 14
        and every_literal_moment_zero
        and len(alphabet) == 69
        and alphabet_digest == EXPECTED_ALPHABET_SHA256
        and fourth_values == list(range(1, P))
        and fourth_given_n2_zero == sorted(FOURTH_ALPHABET) == [7, 8, 11]
        and len(five_root_ranks) == comb(14, 5) == 2002
        and set(five_root_ranks) == {5}
    )
    _require(proved, "the literal-root or hard-alphabet premise changed")
    return {
        "field": "F_13",
        "projective_direction_count": len(POINTS),
        "literal_star_even_degrees": list(degrees),
        "literal_star_power_sum_residue_sets": {
            str(degree): sorted(set(by_center.values()))
            for degree, by_center in literal_star_sums.items()
        },
        "every_exact_literal_direction_is_a_common_M2_M4_M6_root": (
            every_literal_moment_zero
        ),
        "hard_moment_alphabet_size": len(alphabet),
        "hard_moment_alphabet_sha256": alphabet_digest,
        "hard_fourth_moment_value_set": fourth_values,
        "hard_fourth_moment_never_zero": 0 not in fourth_values,
        "hard_fourth_values_when_N2_zero": fourth_given_n2_zero,
        "five_direction_subsets_checked": len(five_root_ranks),
        "five_root_quartic_evaluation_rank_set": sorted(set(five_root_ranks)),
        "five_roots_force_M4_identically_zero": True,
        "z_at_least_five_impossible": True,
        "proved": proved,
    }


def _audit_z2_survivor(
    survivor: dict[str, object],
    hard: tuple[int, ...],
    opposite: tuple[int, ...],
    alphabet: frozenset[tuple[int, int, int]],
) -> bool:
    roots = tuple(int(value) for value in survivor["roots"])
    scalar = int(survivor["M2_scalar"])
    quadratic = tuple(
        int(value) for value in survivor["M4_quotient_coefficients"]
    )
    n2 = tuple(int(value) for value in survivor["hard_N2"])
    n4 = tuple(int(value) for value in survivor["hard_N4"])
    n6 = tuple(int(value) for value in survivor["hard_N6"])
    if not (
        len(roots) == 2
        and set(roots) <= set(opposite)
        and len(quadratic) == 3
        and len(n2) == len(n4) == len(n6) == len(hard) == 7
    ):
        return False
    root_values = root_product_values(roots)
    expected_n2 = tuple(scalar * root_values[index] % P for index in hard)
    expected_n4 = tuple(
        root_values[index]
        * sum(
            coefficient
            * pow(POINTS[index][0], 2 - monomial, P)
            * pow(POINTS[index][1], monomial, P)
            for monomial, coefficient in enumerate(quadratic)
        )
        % P
        for index in hard
    )
    if n2 != expected_n2 or n4 != expected_n4:
        return False
    if not all(row in alphabet for row in zip(n2, n4, n6)):
        return False
    quotient_values = tuple(
        value * pow(root_values[index], -1, P) % P
        for value, index in zip(n6, hard)
    )
    evaluation, inverse = quartic_code_membership_setup(hard)
    return quartic_code_contains(quotient_values, evaluation, inverse)


@lru_cache(maxsize=1)
def exact_literal_interpolation_certificate() -> dict[str, object]:
    """Package and independently replay the exact z=4,3,2 enumeration."""
    raw = audit()
    alphabet = hard_moment_alphabet()
    sign_rows: dict[str, dict[str, object]] = {}
    rows_proved = []
    for hard_sign in (-1, 1):
        source = raw["sign_rows"][str(hard_sign)]
        hard = tuple(int(value) for value in source["hard_direction_indices"])
        opposite = tuple(
            int(value) for value in source["opposite_direction_indices"]
        )
        z2 = list(source["z2_survivors"])
        z2_digest = _catalog_digest(z2)
        replayed = all(
            _audit_z2_survivor(survivor, hard, opposite, alphabet)
            for survivor in z2
        )
        unique_count = len(
            {json.dumps(row, sort_keys=True) for row in z2}
        )
        row_proved = bool(
            len(hard) == len(opposite) == 7
            and set(hard).isdisjoint(opposite)
            and set(hard) | set(opposite) == set(range(14))
            and source["z4_survivor_count"] == 0
            and source["z3_survivor_count"] == 0
            and source["z4_survivors"] == []
            and source["z3_survivors"] == []
            and source["z2_M2_M4_candidate_count"] == 1554
            and source["z2_N6_vectors_checked"] == 2688
            and source["z2_survivor_count"] == 336
            and len(z2) == unique_count == 336
            and z2_digest == EXPECTED_Z2_SHA256[hard_sign]
            and replayed
        )
        _require(row_proved, f"the hard-sign {hard_sign} interpolation changed")
        sign_rows[str(hard_sign)] = {
            "hard_sign": hard_sign,
            "hard_direction_indices": list(hard),
            "opposite_direction_indices": list(opposite),
            "z4_parameter_cases_before_alphabet_filter": comb(7, 4) * (P - 1),
            "z4_survivor_count": 0,
            "z3_parameter_cases_before_alphabet_filter": (
                comb(7, 3) * (P**2 - 1)
            ),
            "z3_survivor_count": 0,
            "z2_M2_M4_candidate_count_after_alphabet_filter": int(
                source["z2_M2_M4_candidate_count"]
            ),
            "z2_N6_vectors_checked": int(source["z2_N6_vectors_checked"]),
            "z2_survivor_count": len(z2),
            "z2_unique_survivor_count": unique_count,
            "z2_survivor_catalog_sha256": z2_digest,
            "representative_z2_survivor": json.loads(
                json.dumps(z2[0], sort_keys=True)
            ),
            "every_z2_survivor_independently_replayed": replayed,
            "proved": row_proved,
        }
        rows_proved.append(row_proved)

    proved = bool(
        raw["hard_moment_alphabet_size"] == 69
        and raw["fourth_moment_alphabet_under_M2_zero"] == [7, 8, 11]
        and raw["payload_sha256"] == EXPECTED_PAYLOAD_SHA256
        and raw["all_z2_z3_z4_interpolation_cases_empty"] is False
        and raw["proved"] is False
        and all(rows_proved)
    )
    _require(proved, "the exact literal interpolation certificate changed")
    return {
        "field": "F_13",
        "runner": "scripts/p13_p5_literal_interpolation.py",
        "hard_moment_alphabet_size": 69,
        "fourth_moment_alphabet_under_M2_zero": [7, 8, 11],
        "raw_payload_sha256": EXPECTED_PAYLOAD_SHA256,
        "sign_rows": sign_rows,
        "z4_empty_for_both_hard_signs": True,
        "z3_empty_for_both_hard_signs": True,
        "z2_survivors_per_hard_sign": 336,
        "all_z2_z3_z4_cases_empty": False,
        "raw_runner_proved_flag_is_false_because_z2_survives": True,
        "result_status": "exhaustive finite interpolation certificate",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p5_excess_partition_reduction() -> dict[str, object]:
    """Deduce that only two literal minima and excess partition 1^5 remain."""
    prior = proposition_15747()
    p5_ledger = t4_u4_catalog_consequence()["family_ledgers"][
        "all_equal_triple"
    ]
    roots = literal_root_and_hard_alphabet_certificate()
    interpolation = exact_literal_interpolation_certificate()
    prior_z_lower_bound = int(prior["P5_minimum_literal_count_at_least"])
    z_minimum = prior_z_lower_bound
    positive_excess_count = OPPOSITE_DIRECTION_COUNT - z_minimum
    surviving_partition = [1] * positive_excess_count
    proved = bool(
        prior["proved"]
        and prior["P5_Q3_minimum_cells_forced_literal"]
        and p5_ledger["common_hard_parallel_count_P"] == 5
        and p5_ledger["minimum_opposite_Q"] == 3
        and p5_ledger["opposite_excess_sum"] == OPPOSITE_EXCESS_SUM
        and p5_ledger["directions_at_minimum_at_least"]
        == prior_z_lower_bound
        == 2
        and roots["proved"]
        and roots["z_at_least_five_impossible"]
        and interpolation["proved"]
        and interpolation["z4_empty_for_both_hard_signs"]
        and interpolation["z3_empty_for_both_hard_signs"]
        and interpolation["z2_survivors_per_hard_sign"] == 336
        and OPPOSITE_EXCESS_SUM == 5
        and positive_excess_count == 5
        and sum(surviving_partition) == OPPOSITE_EXCESS_SUM
    )
    _require(proved, "the P=5 excess-partition reduction failed")
    return {
        "p": 13,
        "t": 4,
        "k": 60,
        "u": 4,
        "hard_family": "all_equal_triple",
        "hard_parallel_count_P": 5,
        "minimum_opposite_parallel_count_Q": 3,
        "z_definition": "number of minimum Q=3 exact-literal opposite directions",
        "prior_prop_15747_forces_every_minimum_Q3_cell_literal": True,
        "opposite_direction_count": OPPOSITE_DIRECTION_COUNT,
        "opposite_excess_sum": OPPOSITE_EXCESS_SUM,
        "prior_lower_bound_on_z": prior_z_lower_bound,
        "z_at_least_5_excluded_by_M4_root_count": True,
        "z4_excluded_by_exact_interpolation": True,
        "z3_excluded_by_exact_interpolation": True,
        "z2_moment_level_survivors_per_hard_sign": 336,
        "forced_z": z_minimum,
        "positive_opposite_excess_count": positive_excess_count,
        "only_remaining_opposite_excess_partition": surviving_partition,
        "partition_identity": "five positive integers sum to five, hence 1^5",
        "moment_level_survivors_are_not_common_graph_realizations": True,
        "P5_branch_closed": False,
        "result_status": "proved open reduction",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15748() -> dict[str, object]:
    """Package the exact interpolation reduction without claiming closure."""
    prior = proposition_15747()
    roots = literal_root_and_hard_alphabet_certificate()
    interpolation = exact_literal_interpolation_certificate()
    reduction = p5_excess_partition_reduction()
    proved = bool(
        prior["proved"]
        and prior["P5_Q3_minimum_cells_forced_literal"]
        and roots["proved"]
        and interpolation["proved"]
        and reduction["proved"]
        and not reduction["P5_branch_closed"]
    )
    _require(proved, "Proposition 15.748 failed")
    return {
        "prop": "15.748",
        "title": "Exact literal-root interpolation in the p13 P=5 branch",
        "result_status": (
            "exhaustive finite interpolation certificate and proved open reduction"
        ),
        "statement": (
            "z>=3 minimum literal cells are impossible; z=2 has 336 "
            "moment-level survivors per hard sign, so only opposite excess "
            "partition 1^5 remains"
        ),
        "prior_prop_15747_dependency": {
            "P5_Q3_minimum_cells_forced_literal": True,
            "proved": bool(prior["proved"]),
        },
        "literal_root_and_hard_alphabet": roots,
        "exact_literal_interpolation": interpolation,
        "P5_excess_partition_reduction": reduction,
        "p13_t4_u4_P5_branch_closed": False,
        "p13_t4_u4_closed": False,
        "p13_k_eq_60_closed": False,
        "remaining_p13_t4_residues": [4, 6],
        "next_exact_gate": (
            "couple the 336 z=2 moment-level survivors per sign to the five "
            "excess-one opposite cells and one common 61-edge graph"
        ),
        "residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "proved_means_exact_reduction_verified_not_all_cases_empty": True,
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    """Atomically write the deterministic Proposition 15.748 artifact."""
    target = path or ROOT / "evidence" / "e1_gmin_m4_prop15748.json"
    write_json_atomic(target, proposition_15748())
    return target


def main() -> None:
    result = proposition_15748()
    path = write_evidence()
    print(
        json.dumps(
            {
                "prop": result["prop"],
                "result_status": result["result_status"],
                "P5_branch_closed": result["p13_t4_u4_P5_branch_closed"],
                "output": str(path),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
