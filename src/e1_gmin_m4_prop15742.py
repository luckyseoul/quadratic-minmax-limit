#!/usr/bin/env python3
r"""Prop. 15.742 -- six-dilate energy close of the generic p=13 row.

Proposition 15.741 reduces the last generic ``p=13,t=3`` partition to ten
integer cyclic-distance rows: three elevated hard rows and seven opposite
rows.  Four exact stars force the common quadratic moment ``M_2`` to vanish.
For a nonexact row ``q=(q_1,...,q_6)``, the already proved data are

* elevated: ``sum q=11``, ``||q||_1<=53``, ``||q||_2^2<=86`` and
  ``c_X.q<=91``;
* opposite: ``sum q=-20``, ``||q||_1<=56``, ``||q||_2^2<=106`` and
  ``c_X.q<=-130``;
* in both cases, ``sum a^2 q_a=0 (mod 13)``.

Only the six multiplicative dilates of the interval seven-set are needed.
An exhaustive six-integer enumeration gives sharp row-energy maxima 31 and
82.  Consequently the ten nonstar rows have energy at most

    3*31 + 7*82 = 667,

contradicting Proposition 15.741's exact Parseval value ``707+26*C>=707``.
The recursive enumeration is independently checked by deterministic
one-worker, nineteen-variable CP-SAT exclusions of energies at least 32 and
83.  This is an aggregate finite certificate, not a directional matrix or
midpoint census; no quartic value code or root-quartet split is used.
"""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from itertools import product
from pathlib import Path

from ortools import __version__ as ORTOOLS_VERSION
from ortools.sat.python import cp_model

from e1_gmin_m4_prop15740 import (
    translated_cut_vector,
    translated_cut_vector_catalog,
)
from e1_gmin_m4_prop15741 import (
    DISTANCES,
    H_EDGE_COUNT,
    P,
    difference_radon_gram_certificate,
    exact_star_moment_certificate,
    quartic_root_rank_certificate,
    six_dilate_cut_energy_certificate,
)


BASE_INTERVAL = tuple(range(7))
DILATE_MULTIPLIERS = tuple(range(1, 7))
ROOT = Path(__file__).resolve().parents[1]

EXPECTED_PRE_CUT_COUNTS = {"elevated": 5844, "opposite": 1704}
EXPECTED_SURVIVING_COUNTS = {"elevated": 30, "opposite": 24}
EXPECTED_ENERGY_MAXIMA = {"elevated": 31, "opposite": 82}
EXPECTED_ROW_DIGESTS = {
    "elevated": "7bff1ebb77ac362b5089b46588f603be812ea4a96fdeaa2f2b52881803b486b5",
    "opposite": "5226c7ee0c44d3cf7e460db2a309a08368667686dfcb2cd45ec24ce932081c1a",
}

ROW_PARAMETERS = {
    "elevated": {
        "parallel_count": 6,
        "total": 11,
        "l1_bound": 53,
        "energy_bound": 86,
        "cut_bound": 91,
    },
    "opposite": {
        "parallel_count": 3,
        "total": -20,
        "l1_bound": 56,
        "energy_bound": 106,
        "cut_bound": -130,
    },
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _row_parameters(kind: str) -> dict[str, int]:
    try:
        return ROW_PARAMETERS[kind]
    except KeyError as exc:
        raise ValueError(f"unknown row kind: {kind}") from exc


def _dilate_subset(multiplier: int) -> tuple[int, ...]:
    return tuple(sorted(multiplier * value % P for value in BASE_INTERVAL))


SIX_DILATE_SUBSETS = tuple(
    _dilate_subset(multiplier) for multiplier in DILATE_MULTIPLIERS
)
SIX_DILATE_CUTS = tuple(
    translated_cut_vector(subset) for subset in SIX_DILATE_SUBSETS
)


def _rows_digest(rows: tuple[tuple[int, ...], ...]) -> str:
    payload = ";".join(",".join(map(str, row)) for row in rows)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


@lru_cache(maxsize=1)
def aggregate_dependencies_certificate() -> dict[str, object]:
    """Audit every imported identity used by the six-row enumeration."""
    stars = exact_star_moment_certificate()
    root_ranks = quartic_root_rank_certificate()
    cut_catalog = translated_cut_vector_catalog()
    prior_energy = six_dilate_cut_energy_certificate()
    radon = difference_radon_gram_certificate()

    canonical_dilate_cuts = {
        translated_cut_vector(subset)
        for subset in prior_energy["dilated_interval_seven_sets"]
    }
    full_cut_catalog = {tuple(row) for row in cut_catalog["vectors"]}

    # If T/h=17, an elevated hard row has signed off-bin total 17-6,
    # while an opposite row has total -(17+3).  Triangle inequality on the
    # same nonparallel edge multiplicities gives ||q||_1 <= 59-P_L.
    signed_total = radon["branch_signed_total_T_over_h"]
    elevated_sum = signed_total - ROW_PARAMETERS["elevated"]["parallel_count"]
    opposite_sum = -(
        signed_total + ROW_PARAMETERS["opposite"]["parallel_count"]
    )
    elevated_l1 = H_EDGE_COUNT - ROW_PARAMETERS["elevated"]["parallel_count"]
    opposite_l1 = H_EDGE_COUNT - ROW_PARAMETERS["opposite"]["parallel_count"]

    # Translation-summing cut_W(X)<=7 or <=-10 over thirteen translates
    # gives the two aggregate cut bounds.
    elevated_cut_bound = P * 7
    opposite_cut_bound = P * -10

    # Every entry of every cut vector is even.  Thus r=Cq is even for
    # integer q.  The elevated slack 91-r is odd positive; the opposite
    # slack -130-r is even nonnegative.  Column sums are 42, so the six
    # slack sums are 84 and 60, exactly as in Proposition 15.741.
    cut_column_sums = tuple(
        sum(row[column] for row in SIX_DILATE_CUTS) for column in range(6)
    )
    elevated_slack_sum = 6 * 91 - 42 * elevated_sum
    opposite_slack_sum = 6 * -130 - 42 * opposite_sum

    proved = bool(
        stars["proved"]
        and stars["all_M2_T3_M4_U4_zero"]
        and root_ranks["proved"]
        and root_ranks["degree_2_four_roots_force_zero"]
        and cut_catalog["proved"]
        and len(SIX_DILATE_CUTS) == len(set(SIX_DILATE_CUTS)) == 6
        and set(SIX_DILATE_CUTS) == canonical_dilate_cuts
        and set(SIX_DILATE_CUTS) <= full_cut_catalog
        and all(value % 2 == 0 for row in SIX_DILATE_CUTS for value in row)
        and cut_column_sums == (42,) * 6
        and elevated_sum == 11
        and opposite_sum == -20
        and elevated_l1 == 53
        and opposite_l1 == 56
        and elevated_cut_bound == 91
        and opposite_cut_bound == -130
        and elevated_slack_sum == 84
        and opposite_slack_sum == 60
        and prior_energy["proved"]
        and prior_energy["elevated_row"]["integer_q_energy_bound"] == 86
        and prior_energy["opposite_row"]["integer_q_energy_bound"] == 106
        and ROW_PARAMETERS["elevated"]["total"] == elevated_sum
        and ROW_PARAMETERS["opposite"]["total"] == opposite_sum
        and ROW_PARAMETERS["elevated"]["l1_bound"] == elevated_l1
        and ROW_PARAMETERS["opposite"]["l1_bound"] == opposite_l1
        and ROW_PARAMETERS["elevated"]["cut_bound"] == elevated_cut_bound
        and ROW_PARAMETERS["opposite"]["cut_bound"] == opposite_cut_bound
        and ROW_PARAMETERS["elevated"]["energy_bound"]
        == prior_energy["elevated_row"]["integer_q_energy_bound"]
        and ROW_PARAMETERS["opposite"]["energy_bound"]
        == prior_energy["opposite_row"]["integer_q_energy_bound"]
        and prior_energy["collision_parameter_upper_bound"] == 11
        and radon["proved"]
        and radon["three_elevated_plus_seven_opposite_off_bin_energy"]
        == "707+26*C"
    )
    _require(proved, "the Proposition 15.742 dependency ledger changed")
    return {
        "p": P,
        "four_exact_stars_force_global_M2_zero": True,
        "M2_dependency": {
            "exact_star_certificate": stars["proved"],
            "degree_2_four_roots_force_zero": root_ranks[
                "degree_2_four_roots_force_zero"
            ],
        },
        "signed_total_T_over_h": signed_total,
        "row_sums": {"elevated": elevated_sum, "opposite": opposite_sum},
        "l1_derivation": "sum_a |q_L(a)| <= number of nonparallel edges = 59-P_L",
        "l1_bounds": {"elevated": elevated_l1, "opposite": opposite_l1},
        "cut_bound_derivation": (
            "sum_t cut_W(X+t)=c_X.q; multiply the cellwise bounds 7 and -10 by p=13"
        ),
        "cut_bounds": {
            "elevated": elevated_cut_bound,
            "opposite": opposite_cut_bound,
        },
        "six_interval_dilate_subsets": [list(row) for row in SIX_DILATE_SUBSETS],
        "six_interval_dilate_cuts": [list(row) for row in SIX_DILATE_CUTS],
        "six_cut_set_matches_15_741": True,
        "six_cuts_in_full_74_catalog": True,
        "cut_image_is_even": True,
        "cut_column_sums": list(cut_column_sums),
        "slacks": {
            "elevated": {
                "definition": "y=91*1-Cq",
                "parity_and_sign": "odd positive",
                "sum": elevated_slack_sum,
            },
            "opposite": {
                "definition": "y=-130*1-Cq",
                "parity_and_sign": "even nonnegative",
                "sum": opposite_slack_sum,
            },
        },
        "imported_integer_energy_bounds": {
            "elevated": prior_energy["elevated_row"]["integer_q_energy_bound"],
            "opposite": prior_energy["opposite_row"]["integer_q_energy_bound"],
        },
        "row_parameters_equal_imported_live_values": True,
        "imported_nonstar_parseval": "707+26*C",
        "collision_parameter_nonnegative": True,
        "proved": proved,
    }


def _remaining_energy_floor(total: int, slots: int) -> int:
    """Minimum square energy of ``slots`` integers with the given sum."""
    magnitude, residue = divmod(abs(total), slots)
    return residue * (magnitude + 1) ** 2 + (slots - residue) * magnitude**2


def _component_range(kind: str) -> range:
    """Exact coordinate range implied only by row sum and energy."""
    parameters = _row_parameters(kind)
    total = parameters["total"]
    energy_bound = parameters["energy_bound"]
    l1_bound = parameters["l1_bound"]
    feasible = [
        value
        for value in range(-l1_bound, l1_bound + 1)
        if value * value + _remaining_energy_floor(total - value, 5)
        <= energy_bound
    ]
    _require(bool(feasible), f"empty {kind} coordinate range")
    return range(min(feasible), max(feasible) + 1)


@lru_cache(maxsize=None)
def enumerate_six_dilate_rows(
    kind: str,
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    """Exhaust all integer rows before and after the six cut inequalities."""
    parameters = _row_parameters(kind)
    total = parameters["total"]
    l1_bound = parameters["l1_bound"]
    energy_bound = parameters["energy_bound"]
    cut_bound = parameters["cut_bound"]
    values = _component_range(kind)
    pre_cut: list[tuple[int, ...]] = []
    surviving: list[tuple[int, ...]] = []
    for prefix in product(values, repeat=5):
        final = total - sum(prefix)
        if final not in values:
            continue
        row = prefix + (final,)
        if sum(abs(value) for value in row) > l1_bound:
            continue
        if sum(value * value for value in row) > energy_bound:
            continue
        if sum(
            distance * distance * value
            for distance, value in zip(DISTANCES, row)
        ) % P:
            continue
        pre_cut.append(row)
        if all(
            sum(coefficient * value for coefficient, value in zip(cut, row))
            <= cut_bound
            for cut in SIX_DILATE_CUTS
        ):
            surviving.append(row)
    return tuple(pre_cut), tuple(surviving)


@lru_cache(maxsize=None)
def independent_cpsat_energy_exclusion(kind: str) -> dict[str, object]:
    """Independently exclude one more unit than the recursive sharp maximum."""
    parameters = _row_parameters(kind)
    l1_bound = parameters["l1_bound"]
    forbidden_floor = EXPECTED_ENERGY_MAXIMA[kind] + 1
    model = cp_model.CpModel()
    # The broad domains follow directly from ||q||_1<=l1_bound; the solver
    # does not consume the recursive catalog or its derived coordinate range.
    q = [
        model.NewIntVar(-l1_bound, l1_bound, f"q_{index}")
        for index in range(6)
    ]
    q_abs = [
        model.NewIntVar(0, l1_bound, f"qabs_{index}") for index in range(6)
    ]
    q_square = [
        model.NewIntVar(0, l1_bound * l1_bound, f"qsq_{index}")
        for index in range(6)
    ]
    for value, absolute, square in zip(q, q_abs, q_square):
        model.AddAbsEquality(absolute, value)
        model.AddMultiplicationEquality(square, [value, value])
    model.Add(sum(q) == parameters["total"])
    model.Add(sum(q_abs) <= l1_bound)
    model.Add(sum(q_square) >= forbidden_floor)
    # |sum a^2*q_a| <= 36*||q||_1, so [-200,200] exhausts both rows
    # (36*56/13 < 156) without consuming the prior energy cap.
    quotient = model.NewIntVar(-200, 200, "M2_quotient")
    model.Add(
        sum(
            distance * distance * value
            for distance, value in zip(DISTANCES, q)
        )
        == P * quotient
    )
    for index, cut in enumerate(SIX_DILATE_CUTS):
        model.Add(
            sum(coefficient * value for coefficient, value in zip(cut, q))
            <= parameters["cut_bound"]
        ).WithName(f"interval_dilate_cut_{index}")

    validation = model.Validate()
    _require(not validation, f"invalid {kind} energy model: {validation}")
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    _require(
        status == cp_model.INFEASIBLE,
        f"the {kind} forbidden-energy model is no longer infeasible",
    )
    return {
        "forbidden_energy_floor": forbidden_floor,
        "status": solver.StatusName(status),
        "infeasible": True,
        "solver": "OR-Tools CP-SAT",
        "solver_version": ORTOOLS_VERSION,
        "workers": 1,
        "seed": 0,
        "model_validation": validation,
        "model_proto_sha256": hashlib.sha256(
            str(model.Proto()).encode("utf-8")
        ).hexdigest(),
        "variables": len(model.Proto().variables),
        "constraints": len(model.Proto().constraints),
        "M2_quotient_domain": [-200, 200],
        "M2_quotient_absolute_bound_from_l1": 36 * l1_bound // P + 1,
        "prior_energy_upper_constraint_used": False,
        "proved": True,
    }


@lru_cache(maxsize=None)
def row_catalog_certificate(kind: str) -> dict[str, object]:
    """Return the exact recursive catalog and its independent upper-bound check."""
    pre_cut, rows = enumerate_six_dilate_rows(kind)
    energies = tuple(sum(value * value for value in row) for row in rows)
    maximum = max(energies)
    maximizers = tuple(row for row, energy in zip(rows, energies) if energy == maximum)
    digest = _rows_digest(rows)
    expected_pre_cut = EXPECTED_PRE_CUT_COUNTS[kind]
    expected_count = EXPECTED_SURVIVING_COUNTS[kind]
    expected_maximum = EXPECTED_ENERGY_MAXIMA[kind]
    expected_digest = EXPECTED_ROW_DIGESTS[kind]
    independent = independent_cpsat_energy_exclusion(kind)
    proved = bool(
        len(pre_cut) == expected_pre_cut
        and len(rows) == expected_count
        and maximum == expected_maximum
        and digest == expected_digest
        and len(maximizers) == 6
        and independent["infeasible"]
        and independent["forbidden_energy_floor"] == maximum + 1
        and independent["prior_energy_upper_constraint_used"] is False
    )
    _require(proved, f"the {kind} exact row catalog changed")
    values = _component_range(kind)
    return {
        "kind": kind,
        "coordinate_range_from_sum_and_energy": [values.start, values.stop - 1],
        "pre_cut_row_count": len(pre_cut),
        "surviving_row_count": len(rows),
        "surviving_rows_sha256": digest,
        "catalog_scope": (
            "exact catalog of the stated necessary six-bin relaxation; "
            "a superset of realizable directional rows"
        ),
        "prior_energy_cap_redundant_by_independent_cpsat": True,
        "surviving_rows": [list(row) for row in rows],
        "sharp_energy_maximum": maximum,
        "maximizer_count": len(maximizers),
        "maximizers": [list(row) for row in maximizers],
        "independent_cpsat_exclusion": independent,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15742() -> dict[str, object]:
    """Package the exhaustive aggregate contradiction and its exact scope."""
    dependencies = aggregate_dependencies_certificate()
    elevated = row_catalog_certificate("elevated")
    opposite = row_catalog_certificate("opposite")
    three_elevated_upper = 3 * elevated["sharp_energy_maximum"]
    seven_opposite_upper = 7 * opposite["sharp_energy_maximum"]
    nonstar_upper = three_elevated_upper + seven_opposite_upper
    parseval_lower = 707
    gap = parseval_lower - nonstar_upper
    proved = bool(
        dependencies["proved"]
        and elevated["proved"]
        and opposite["proved"]
        and nonstar_upper == 667
        and gap == 40
        and nonstar_upper < parseval_lower
    )
    _require(proved, "Proposition 15.742 aggregate contradiction failed")
    return {
        "prop": "15.742",
        "title": "six-dilate energy close of the generic p=13 fourth shell",
        "result_status": "exhaustive finite certificate",
        "p": P,
        "layer_index_t": 3,
        "original_k": 4 * P + 6,
        "remaining_hard_quotient_partition": [1, 1, 1, 1, 2, 2, 2],
        "dependencies": dependencies,
        "dependency_chain": {
            "15.739": "closes the exceptional p=13,t=3,u=3 branch",
            "15.740": "leaves only the generic four-exact partition 1^4 2^3",
            "15.741": (
                "supplies M2=0, the six-cut row energy bounds, and exact "
                "nonstar Parseval energy 707+26*C"
            ),
        },
        "model_scope": "six integer cyclic-distance bins per nonexact direction",
        "elevated_row_catalog": elevated,
        "opposite_row_catalog": opposite,
        "global_energy": {
            "three_elevated_upper": three_elevated_upper,
            "seven_opposite_upper": seven_opposite_upper,
            "nonstar_upper": nonstar_upper,
            "exact_parseval": "707+26*C",
            "collision_parameter_lower_bound": 0,
            "parseval_lower": parseval_lower,
            "gap": gap,
            "contradiction": True,
        },
        "quartic_code_used": False,
        "root_quartet_case_split_used": False,
        "hard_sign_normalization_used": False,
        "directional_coefficient_matrix_census_used": False,
        "binary_midpoint_lift_used": False,
        "common_graph_exists": False,
        "p13_generic_four_exact_partition_closed": True,
        "p13_generic_t3_branch_closed": True,
        "p13_t3_exceptional_u3_closed_by_15_739": True,
        "p13_k_eq_58_closed": True,
        "generic_p_ge_17_t3_branch_closed": False,
        "k_eq_4p_plus_6_shell_closed": False,
        "residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "remaining_scope": (
            "critical p=5,7; p=11 at k>=50; p=13 at k>=60; generic p>=17 "
            "fourth and later residual layers; multi-level Type I; and the limit"
        ),
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    """Write the canonical exhaustive-certificate payload."""
    target = path or ROOT / "evidence" / "e1_gmin_m4_prop15742.json"
    target.write_text(json.dumps(proposition_15742(), indent=2, sort_keys=True) + "\n")
    return target


def main() -> None:
    theorem = proposition_15742()
    target = write_evidence()
    print(
        "Prop. 15.742: six-dilate row energy 667 contradicts "
        "common-graph energy >=707"
    )
    print(f"result status={theorem['result_status']}")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
