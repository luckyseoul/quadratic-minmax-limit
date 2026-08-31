#!/usr/bin/env python3
r"""Prop. 15.747 -- remove the mass-12 lifts in the p13 u=4 gate.

Proposition 15.746 leaves minimum opposite cells of scaled mean twelve at
``Q=5`` in the omitted-pair ``P=3`` branch and at ``Q=3`` in the
all-equal-triple ``P=5`` branch.  A phase-zero ``b=0`` cell has ``A=2C``
with ``4p E[C]=12`` and height one or four.

Height one is impossible for every integral parallel count by one exact cut
second moment.  Height four at Q=3 and Q=5 is excluded by deterministic
one-worker CP-SAT models using only necessary coefficient, cut, parity, and
edge-budget constraints.  No field-moment identity is needed.  Consequently
the P=3 branch is empty, while every minimum Q=3 cell in the P=5 branch is
the exact b=12 literal.
"""
from __future__ import annotations

import hashlib
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path

from ortools import __version__ as ORTOOLS_VERSION
from ortools.sat.python import cp_model

from e1_gmin_m4_prop15738 import middle_slice_points, pair_coordinates
from e1_gmin_m4_prop15746 import (
    H_EDGE_COUNT,
    P,
    SUPPORT396,
    mass12_phase_zero_dichotomy,
    t4_u4_catalog_consequence,
)
from io_atomic import write_json_atomic


ROOT = Path(__file__).resolve().parents[1]
SCALED_MASS = 12
HEIGHT = 4
DOMAIN_SIZE = 1716
PAIR_COUNT = 78
VALUE_SUM = SUPPORT396
PARALLEL_COUNTS = (3, 5)
PAIR_SEPARATION_MULTIPLICITY = 2 * comb(P - 2, 6)
ANCHOR_INDEX = 0
SEARCH_WORKERS = 1
EXPECTED_MODEL_SHA256 = {
    3: "e8404a5684e033b73750b1f36a338aa13038861d6dbfc614cc99b6f0666423d9",
    5: "8f992368fac869f29c23e6ecd20400228c2c10d5bda4d1001b291242dd6e3941",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _model_textproto_sha256(model: cp_model.CpModel) -> str:
    return hashlib.sha256(str(model.Proto()).encode("utf-8")).hexdigest()


def coefficient_sum(parallel_count: int) -> int:
    """Return the forced normalized coefficient sum for mass twelve."""
    return P * (parallel_count - 3) - SCALED_MASS


@lru_cache(maxsize=1)
def mass12_boolean_second_moment_exclusion() -> dict[str, object]:
    """Exclude every Boolean phase-zero b=0 mass-12 lift at p=13."""
    edge_cut_probability = Fraction(7, 13)
    adjacent_pair_probability = Fraction(7, 26)
    disjoint_pair_probability = Fraction(42, 143)
    rows = {}
    for q in PARALLEL_COUNTS:
        total_w = coefficient_sum(q)
        cut_without_c = 7 * q - 27
        required_second_moment = Fraction(
            (DOMAIN_SIZE - VALUE_SUM) * cut_without_c**2
            + VALUE_SUM * (cut_without_c - 2) ** 2,
            DOMAIN_SIZE,
        )
        residual = 182 * q * q - 1428 * q + 2598
        rows[str(q)] = {
            "parallel_count_Q": q,
            "coefficient_sum_S": total_w,
            "cut_values": [cut_without_c - 2, cut_without_c],
            "low_cut_multiplicity": VALUE_SUM,
            "high_cut_multiplicity": DOMAIN_SIZE - VALUE_SUM,
            "required_cut_second_moment": str(required_second_moment),
            "integer_energy_equation": (
                f"-7*D2+84*E2{residual:+d}=0"
            ),
            "residual_mod_7": residual % 7,
            "impossible_mod_7": residual % 7 != 0,
        }

    # For a uniform six-set, write A for the sum of products of adjacent
    # distinct edge weights and B for the analogous disjoint sum.  Then
    # D2=2E2+2A and S^2=E2+2A+2B.  Substitution of the three probabilities
    # gives (-7D2+84E2+84S^2)/286 exactly.
    general_residues = {
        (182 * q * q - 1428 * q + 2598) % 7 for q in range(7)
    }
    proved = bool(
        edge_cut_probability == Fraction(7, 13)
        and adjacent_pair_probability == Fraction(7, 26)
        and disjoint_pair_probability == Fraction(42, 143)
        and VALUE_SUM == 396
        and general_residues == {1}
        and all(row["impossible_mod_7"] for row in rows.values())
        and rows["3"]["required_cut_second_moment"] == "552/13"
        and rows["5"]["required_cut_second_moment"] == "748/13"
    )
    _require(proved, "the mass-12 Boolean cut obstruction changed")
    return {
        "p": P,
        "uniform_cut_side_size": 6,
        "same_cuts_as_slice": "six-sets are complements of J(13,7)",
        "one_edge_cut_probability": str(edge_cut_probability),
        "two_adjacent_edges_both_cut_probability": str(
            adjacent_pair_probability
        ),
        "two_disjoint_edges_both_cut_probability": str(
            disjoint_pair_probability
        ),
        "notation": {
            "S": "sum_e W_e",
            "E2": "sum_e W_e^2",
            "D2": "sum_v (sum_(e incident v) W_e)^2",
        },
        "weighted_cut_second_moment_identity": (
            "E[cut_W^2]=(-7*D2+84*E2+84*S^2)/286"
        ),
        "coefficient_sum_formula": "S=13*(Q-3)-12=13Q-51",
        "boolean_cut_formula": "cut_W=7Q-27-2C",
        "boolean_support_size": VALUE_SUM,
        "general_integer_energy_equation": (
            "-7*D2+84*E2+182Q^2-1428Q+2598=0"
        ),
        "general_residual_is_one_mod_7": True,
        "parallel_rows": rows,
        "l1_row_parity_third_differences_and_field_moments_needed": False,
        "all_integral_Q_boolean_mass12_lifts_excluded": True,
        "result_status": "proved theorem",
        "proved": proved,
    }


@lru_cache(maxsize=None)
def mass12_height_four_arithmetic(parallel_count: int) -> dict[str, object]:
    """Derive the projected height-four model constants at Q=3 or Q=5."""
    if parallel_count not in PARALLEL_COUNTS:
        raise ValueError("parallel_count must be 3 or 5")
    total_w = coefficient_sum(parallel_count)
    l1_budget = H_EDGE_COUNT - parallel_count
    offset = parallel_count - 3 + total_w
    cut_lower = (offset - 16) // 2
    cut_upper = offset // 2
    cut_total = PAIR_SEPARATION_MULTIPLICITY * total_w
    derived_value_sum = (DOMAIN_SIZE * offset - 2 * cut_total) // 4
    proved = bool(
        offset % 2 == 0
        and cut_lower * 2 == offset - 16
        and cut_upper * 2 == offset
        and derived_value_sum == VALUE_SUM
        and PAIR_SEPARATION_MULTIPLICITY == 924
    )
    _require(proved, f"Q={parallel_count} height-four arithmetic changed")
    return {
        "parallel_count_Q": parallel_count,
        "coefficient_sum": total_w,
        "coefficient_sum_formula": "sum W=13*(Q-3)-12",
        "l1_budget": l1_budget,
        "l1_formula": "sum |W|<=61-Q",
        "cut_identity": (
            f"4C(X)={offset}-2*cut_W(X)"
        ),
        "cut_lower_at_C4": cut_lower,
        "cut_upper_at_C0": cut_upper,
        "even_row_sums_make_every_cut_even": True,
        "pair_separation_multiplicity": PAIR_SEPARATION_MULTIPLICITY,
        "derived_value_sum": derived_value_sum,
        "height_four_anchor_cut": cut_lower,
        "height_four_anchor_is_wlog_before_field_moments": True,
        "anchor_symmetry": "S_13 transitive on J(13,7)",
        "field_moment_constraint_used": False,
        "proved": proved,
    }


def build_mass12_height_four_model(
    parallel_count: int,
) -> tuple[cp_model.CpModel, dict[tuple[int, int], cp_model.IntVar], dict[str, object]]:
    """Build the exact necessary coefficient/cut relaxation."""
    arithmetic = mass12_height_four_arithmetic(parallel_count)
    points = middle_slice_points()
    pairs = pair_coordinates()
    l1_budget = int(arithmetic["l1_budget"])
    model = cp_model.CpModel()
    weights = {
        pair: model.NewIntVar(-l1_budget, l1_budget, f"W_{pair[0]}_{pair[1]}")
        for pair in pairs
    }
    absolutes = {
        pair: model.NewIntVar(0, l1_budget, f"absW_{pair[0]}_{pair[1]}")
        for pair in pairs
    }
    model.Add(sum(weights.values()) == arithmetic["coefficient_sum"])
    for pair in pairs:
        model.AddAbsEquality(absolutes[pair], weights[pair])
    model.Add(sum(absolutes.values()) <= l1_budget)
    row_halves = [
        model.NewIntVar(-l1_budget, l1_budget, f"row_half_{vertex}")
        for vertex in range(P)
    ]
    for vertex in range(P):
        model.Add(
            sum(weight for pair, weight in weights.items() if vertex in pair)
            == 2 * row_halves[vertex]
        )
    for index, point in enumerate(points):
        point_set = set(point)
        cut = sum(
            weight
            for (left, right), weight in weights.items()
            if (left in point_set) != (right in point_set)
        )
        model.Add(cut >= arithmetic["cut_lower_at_C4"])
        model.Add(cut <= arithmetic["cut_upper_at_C0"])
        if index == ANCHOR_INDEX:
            model.Add(cut == arithmetic["height_four_anchor_cut"])
    validation = model.Validate()
    _require(not validation, f"invalid Q={parallel_count} model: {validation}")
    metadata = {
        "arithmetic": arithmetic,
        "coefficient_variable_count": len(weights),
        "absolute_value_variable_count": len(absolutes),
        "row_half_variable_count": len(row_halves),
        "integer_variable_count": len(model.Proto().variables),
        "constraint_count": len(model.Proto().constraints),
        "cut_lower_inequality_count": len(points),
        "cut_upper_inequality_count": len(points),
        "height_anchor_constraint_count": 1,
        "model_validation": validation,
        "model_textproto_sha256": _model_textproto_sha256(model),
        "model_is_exact_projection_of_value_formulation": True,
    }
    return model, weights, metadata


@lru_cache(maxsize=None)
def mass12_height_four_exclusion(parallel_count: int) -> dict[str, object]:
    """Run one deterministic exact height-four exclusion."""
    model, _weights, metadata = build_mass12_height_four_model(parallel_count)
    _require(
        metadata["model_textproto_sha256"]
        == EXPECTED_MODEL_SHA256[parallel_count],
        f"Q={parallel_count} model hash changed",
    )
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = SEARCH_WORKERS
    solver.parameters.random_seed = 0
    solver.parameters.randomize_search = False
    solver.parameters.symmetry_level = 3
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    infeasible = status == cp_model.INFEASIBLE
    _require(infeasible, f"Q={parallel_count} model is {solver.StatusName(status)}")
    return {
        "parallel_count_Q": parallel_count,
        "height_under_test": HEIGHT,
        "arithmetic": metadata["arithmetic"],
        "model": {key: value for key, value in metadata.items() if key != "arithmetic"},
        "solver": {
            "name": "OR-Tools CP-SAT",
            "version": ORTOOLS_VERSION,
            "status": solver.StatusName(status),
            "status_code": int(status),
            "num_search_workers": SEARCH_WORKERS,
            "deterministic_one_worker_replay": True,
            "no_timeout": True,
            "branches": solver.NumBranches(),
            "conflicts": solver.NumConflicts(),
        },
        "height_four_model_infeasible": infeasible,
        "field_moment_constraint_needed": False,
        "result_status": "exhaustive finite necessary-relaxation certificate",
        "proved": infeasible,
    }


@lru_cache(maxsize=1)
def proposition_15747() -> dict[str, object]:
    """Close P=3 and rigidify the minimum cells in P=5."""
    dependency = t4_u4_catalog_consequence()
    dichotomy = mass12_phase_zero_dichotomy()
    boolean = mass12_boolean_second_moment_exclusion()
    height_four = {
        str(q): mass12_height_four_exclusion(q) for q in PARALLEL_COUNTS
    }
    ledgers = dependency["family_ledgers"]
    p3 = ledgers["omitted_pair"]
    p5 = ledgers["all_equal_triple"]
    proved = bool(
        dependency["proved"]
        and dichotomy["proved"]
        and boolean["proved"]
        and all(row["proved"] for row in height_four.values())
        and p3["minimum_opposite_Q"] == 5
        and p3["directions_at_minimum_at_least"] >= 2
        and not p3["b12_literal_compatible_at_minimum_Q"]
        and p5["minimum_opposite_Q"] == 3
        and p5["directions_at_minimum_at_least"] >= 2
        and p5["b12_literal_compatible_at_minimum_Q"]
    )
    _require(proved, "Proposition 15.747 failed")
    return {
        "prop": "15.747",
        "title": "Mass-12 cut obstruction closes the P3 half of p13 u4",
        "result_status": "proved branch exclusion with exact finite certificates",
        "statement": (
            "the p13,t4,u4 omitted-pair P=3 branch is empty; in the "
            "all-equal-triple P=5 branch every minimum Q=3 cell is b=12 literal"
        ),
        "boolean_mass12_exclusion": boolean,
        "height_four_exclusions": height_four,
        "p13_t4_u4_P3_branch_closed": True,
        "P5_Q3_minimum_cells_forced_literal": True,
        "P5_minimum_literal_count_at_least": p5["directions_at_minimum_at_least"],
        "p13_t4_u4_closed": False,
        "remaining_u4_hard_family": "all_equal_triple P=5",
        "remaining_p13_t4_residues": [4, 6],
        "next_exact_gate": (
            "use the common projective roots supplied by the forced Q=3 literals "
            "against the P=5 hard moment alphabet"
        ),
        "residual_ii_closed": False,
        "quadratic_minmax_limit_closed": False,
        "proved": proved,
    }


def write_evidence() -> Path:
    target = ROOT / "evidence" / "e1_gmin_m4_prop15747.json"
    write_json_atomic(target, proposition_15747())
    return target


if __name__ == "__main__":
    result = proposition_15747()
    print("Prop 15.747: P=3 branch CLOSED; P=5 minimum cells are literals")
    print(f"wrote {write_evidence()}")
