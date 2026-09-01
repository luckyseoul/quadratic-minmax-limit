#!/usr/bin/env python3
r"""Prop. 15.753 -- the p=17,19 fifth-shell endpoints are empty.

This is an exact aggregate-row certificate for the two endpoints not covered
by Proposition 15.752.  At ``t=4`` the isolated-chart arithmetic of
Propositions 15.734--15.735 leaves two branches at each prime.  Full
translation-averaged middle-slice cuts, common difference-Radon Parseval,
and the global even moments exclude all four branches.

The only sign-sensitive step occurs in the last partition of the p=17 XNOR
branch.  If ``S_j=epsilon_L M_j`` with ``epsilon_L=h`` on hard directions
and ``epsilon_L=-h`` on opposite directions, the binary quartic

    G(L) = h M_4(L) - M_2(L)^2

vanishes on eight exact hard XNOR directions and hence identically.  Thus
``S_4=S_2^2`` on hard rows but ``S_4=-S_2^2`` on opposite rows.  The latter
minus sign is enforced explicitly below and covered by a regression test.

Every finite row bound is certified by an explicit integral maximizer and
an exact one-worker CP-SAT infeasibility replay at energy one larger.  The
two empty row systems are replayed directly.  This is a finite-prime
aggregate-row census, but not a graph, orbit, coefficient-cell, or common-
realization census.
"""
from __future__ import annotations

import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from typing import Iterable

from ortools import __version__ as ORTOOLS_VERSION
from ortools.sat.python import cp_model

from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    BRANCH_P1_LAST,
    BRANCH_P3_LAST,
    baseline_coefficient_rules,
    residual_even_floor_table,
)
from e1_gmin_m4_prop15743 import translated_cut_vectors as p17_cut_vectors
from io_atomic import write_json_atomic


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE_PATH = ROOT / "evidence" / "e1_gmin_m4_prop15753.json"

EXPECTED_CUT_CATALOGS = {
    17: (
        698,
        72,
        "a8ac7349cb601db5163ef1526949587c766914d774fe26858fe93eac1d940708",
    ),
    19: (
        2338,
        90,
        "5f07e9ced107e6dc1551b806043a92147c00d80eb009b70d0cbfd3ce9631c5b7",
    ),
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical_sha256(value: object) -> str:
    return _sha256_bytes(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )


def _rows_digest(rows: Iterable[tuple[int, ...]]) -> str:
    payload = ";".join(",".join(map(str, row)) for row in rows)
    return _sha256_bytes(payload.encode("ascii"))


def translated_cut_vector(p: int, subset: Iterable[int]) -> tuple[int, ...]:
    """Return the nonzero distance-bin cut vector on the middle slice."""
    if p not in (17, 19):
        raise ValueError("this certificate uses only p=17 or p=19")
    m = (p + 1) // 2
    values = tuple(sorted(int(value) for value in subset))
    if len(values) != m or len(set(values)) != m or not all(
        0 <= value < p for value in values
    ):
        raise ValueError(f"need {m} distinct elements of F_{p}")
    chosen = set(values)
    return tuple(
        sum(
            (value in chosen) != ((value + distance) % p in chosen)
            for value in range(p)
        )
        for distance in range(1, (p + 1) // 2)
    )


@lru_cache(maxsize=2)
def translated_cut_vectors(p: int) -> tuple[tuple[int, ...], ...]:
    """Generate the exact distinct translated-cut catalog for p=17 or 19."""
    if p == 17:
        return p17_cut_vectors()
    if p != 19:
        raise ValueError("this certificate uses only p=17 or p=19")
    m = (p + 1) // 2
    return tuple(
        sorted(
            {
                translated_cut_vector(p, subset)
                for subset in combinations(range(p), m)
            }
        )
    )


@lru_cache(maxsize=2)
def cut_catalog_certificate(p: int) -> dict[str, object]:
    vectors = translated_cut_vectors(p)
    expected_count, expected_sum, expected_hash = EXPECTED_CUT_CATALOGS[p]
    digest = _rows_digest(vectors)
    proved = bool(
        len(vectors) == expected_count
        and digest == expected_hash
        and all(
            len(row) == (p - 1) // 2
            and sum(row) == expected_sum
            and all(value % 2 == 0 and 0 <= value <= p - 1 for value in row)
            for row in vectors
        )
    )
    _require(proved, f"the p={p} translated-cut catalog changed")
    return {
        "p": p,
        "slice_size": (p + 1) // 2,
        "middle_slice_point_count": len(tuple(combinations(range(p), (p + 1) // 2))),
        "distinct_translated_cut_vectors": len(vectors),
        "every_vector_sum": expected_sum,
        "catalog_sha256": digest,
        "proved": proved,
    }


@dataclass(frozen=True)
class RowSpec:
    name: str
    p: int
    branch: str
    direction: str
    row_kind: str
    total: int
    l1_bound: int
    cut_upper: int
    moment_relation: str
    expected_energy: int | None
    witness: tuple[int, ...] | None


ROW_SPECS = (
    # p=17, branch A, raw rows used in every partition except (5).
    RowSpec("p17_A_hard_e1_raw", 17, "A_XNOR", "hard", "excess_1", 0, 72, 17, "none", 28, (-1, 0, 3, -2, 0, 2, -3, 1)),
    RowSpec("p17_A_hard_e2_raw", 17, "A_XNOR", "hard", "excess_2", -1, 71, 17, "none", 81, (-4, 5, -4, 3, -3, 2, -1, 1)),
    RowSpec("p17_A_hard_e3_raw", 17, "A_XNOR", "hard", "excess_3", -2, 70, 17, "none", 200, (1, -7, -2, 7, 3, -6, -4, 6)),
    RowSpec("p17_A_hard_e4_raw", 17, "A_XNOR", "hard", "excess_4", -3, 69, 17, "none", 289, (-5, -2, 8, -7, 1, 7, -9, 4)),
    RowSpec("p17_A_opposite_Q4_raw", 17, "A_XNOR", "opposite", "Q4", -9, 73, -68, "none", 23, (-2, -2, 0, 0, -3, -1, 1, -2)),
    # The last p=17 A partition has eight exact hard XNOR roots.
    RowSpec("p17_A_hard_e5_quartic", 17, "A_XNOR", "hard", "excess_5", -4, 68, 17, "hard_quartic", 384, (0, -10, -2, 10, 4, -8, -6, 8)),
    RowSpec("p17_A_opposite_Q4_quartic", 17, "A_XNOR", "opposite", "Q4", -9, 73, -68, "opposite_quartic", 11, (-1, -1, -2, -1, -1, -1, -1, -1)),
    # p=17, branch B.  At least five exact literal stars give M2=M4=0.
    RowSpec("p17_B_hard_e1_zero", 17, "B_LITERAL", "hard", "excess_1", 15, 71, 153, "zero", None, None),
    RowSpec("p17_B_hard_e2_zero", 17, "B_LITERAL", "hard", "excess_2", 14, 70, 153, "zero", 70, (2, -1, 1, 5, 3, -2, 1, 5)),
    RowSpec("p17_B_hard_e4_zero", 17, "B_LITERAL", "hard", "excess_4", 12, 68, 153, "zero", 218, (3, -5, 1, 8, 6, -5, -3, 7)),
    RowSpec("p17_B_opposite_Q3_zero", 17, "B_LITERAL", "opposite", "Q3", -24, 74, -204, "zero", 72, (-3, -3, -3, -3, -3, -3, -3, -3)),
    RowSpec("p17_B_opposite_Q4_zero", 17, "B_LITERAL", "opposite", "Q4", -25, 73, -204, "zero", 101, (-6, -1, -5, -1, -3, -3, -4, -2)),
    # p=19, branch A.  Raw translated-cut energy already suffices.
    RowSpec("p19_A_hard_e1_raw", 19, "A_XNOR", "hard", "excess_1", 0, 80, 19, "none", 36, (0, -3, 1, 3, -1, -2, 2, 2, -2)),
    RowSpec("p19_A_hard_e2_raw", 19, "A_XNOR", "hard", "excess_2", -1, 79, 19, "none", 97, (4, 1, -2, -4, -4, -3, -1, 3, 5)),
    RowSpec("p19_A_hard_e3_raw", 19, "A_XNOR", "hard", "excess_3", -2, 78, 19, "none", 194, (-5, 4, -1, -3, 5, -7, 7, -4, 2)),
    RowSpec("p19_A_hard_e4_raw", 19, "A_XNOR", "hard", "excess_4", -3, 77, 19, "none", 325, (9, 7, 5, 2, -1, -4, -6, -7, -8)),
    RowSpec("p19_A_hard_e5_raw", 19, "A_XNOR", "hard", "excess_5", -4, 76, 19, "none", 494, (-1, -11, 3, 10, -4, -9, 6, 9, -7)),
    RowSpec("p19_A_opposite_Q4_raw", 19, "A_XNOR", "opposite", "Q4", -9, 81, -76, "none", 23, (-2, 0, -1, -2, 0, -3, 1, -2, 0)),
    # p=19, branch C.  Five exact complement literals give M2=M4=0.
    RowSpec("p19_C_opposite_Q5_zero", 19, "C_COMPLEMENT_LITERAL", "opposite", "Q5", 10, 80, 114, "zero", None, None),
)
ROW_SPEC_BY_NAME = {spec.name: spec for spec in ROW_SPECS}


def moment_residues(p: int, row: Iterable[int]) -> tuple[int, int]:
    values = tuple(int(value) for value in row)
    return tuple(
        sum(pow(distance, degree, p) * value for distance, value in enumerate(values, 1)) % p
        for degree in (2, 4)
    )  # type: ignore[return-value]


def moment_relation_holds(spec: RowSpec, row: Iterable[int]) -> bool:
    m2, m4 = moment_residues(spec.p, row)
    if spec.moment_relation == "none":
        return True
    if spec.moment_relation == "zero":
        return m2 == 0 and m4 == 0
    if spec.moment_relation == "hard_quartic":
        return m4 == m2 * m2 % spec.p
    if spec.moment_relation == "opposite_quartic":
        return m4 == -(m2 * m2) % spec.p
    raise ValueError(f"unknown moment relation {spec.moment_relation!r}")


def validate_witness(spec: RowSpec, witness: Iterable[int]) -> dict[str, object]:
    row = tuple(int(value) for value in witness)
    if len(row) != (spec.p - 1) // 2:
        raise ArithmeticError(f"{spec.name}: wrong witness length")
    energy = sum(value * value for value in row)
    cut_values = tuple(
        sum(coefficient * value for coefficient, value in zip(cut, row))
        for cut in translated_cut_vectors(spec.p)
    )
    m2, m4 = moment_residues(spec.p, row)
    proved = bool(
        spec.expected_energy is not None
        and sum(row) == spec.total
        and sum(abs(value) for value in row) <= spec.l1_bound
        and max(cut_values) <= spec.cut_upper
        and moment_relation_holds(spec, row)
        and energy == spec.expected_energy
    )
    _require(proved, f"{spec.name}: explicit maximizer failed exact replay")
    return {
        "row": list(row),
        "sum": sum(row),
        "l1": sum(abs(value) for value in row),
        "energy": energy,
        "maximum_translated_cut": max(cut_values),
        "M2_mod_p": m2,
        "M4_mod_p": m4,
        "moment_relation": spec.moment_relation,
        "proved": proved,
    }


def _add_moment_constraints(
    model: cp_model.CpModel, spec: RowSpec, values: list[cp_model.IntVar]
) -> None:
    if spec.moment_relation == "none":
        return
    residues: list[cp_model.IntVar] = []
    for degree in (2, 4):
        coefficients = [
            pow(distance, degree, spec.p)
            for distance in range(1, (spec.p + 1) // 2)
        ]
        residue = model.NewIntVar(0, spec.p - 1, f"M{degree}_residue")
        # The crude quotient bounds are far wider than the exact l1 bound and
        # therefore do not remove any integral row.
        quotient = model.NewIntVar(-10_000, 10_000, f"M{degree}_quotient")
        model.Add(
            sum(coefficient * value for coefficient, value in zip(coefficients, values))
            == spec.p * quotient + residue
        )
        residues.append(residue)
    m2, m4 = residues
    if spec.moment_relation == "zero":
        model.Add(m2 == 0)
        model.Add(m4 == 0)
        return
    allowed = [
        (
            residue,
            (residue * residue if spec.moment_relation == "hard_quartic" else -residue * residue)
            % spec.p,
        )
        for residue in range(spec.p)
    ]
    model.AddAllowedAssignments([m2, m4], allowed)


def build_row_model(spec: RowSpec) -> tuple[cp_model.CpModel, list[cp_model.IntVar]]:
    """Build the exact infeasibility model proving one row bound."""
    model = cp_model.CpModel()
    q = (spec.p - 1) // 2
    values = [
        model.NewIntVar(-spec.l1_bound, spec.l1_bound, f"q_{distance}")
        for distance in range(1, q + 1)
    ]
    absolutes = [
        model.NewIntVar(0, spec.l1_bound, f"abs_{distance}")
        for distance in range(1, q + 1)
    ]
    squares = [
        model.NewIntVar(0, spec.l1_bound * spec.l1_bound, f"sq_{distance}")
        for distance in range(1, q + 1)
    ]
    for value, absolute, square in zip(values, absolutes, squares):
        model.AddAbsEquality(absolute, value)
        model.AddMultiplicationEquality(square, [value, value])
    model.Add(sum(values) == spec.total)
    model.Add(sum(absolutes) <= spec.l1_bound)
    _add_moment_constraints(model, spec, values)
    for index, cut in enumerate(translated_cut_vectors(spec.p)):
        model.Add(
            sum(coefficient * value for coefficient, value in zip(cut, values))
            <= spec.cut_upper
        ).WithName(f"translated_cut_{index}")
    if spec.expected_energy is not None:
        model.Add(sum(squares) >= spec.expected_energy + 1).WithName(
            "energy_strictly_above_claimed_maximum"
        )
    return model, values


def row_model_sha256(spec: RowSpec) -> str:
    model, _ = build_row_model(spec)
    return _sha256_bytes(str(model.Proto()).encode("utf-8"))


@lru_cache(maxsize=None)
def replay_row_certificate(name: str) -> dict[str, object]:
    """Replay one bound/emptiness claim exactly with one CP-SAT worker."""
    spec = ROW_SPEC_BY_NAME[name]
    witness = (
        validate_witness(spec, spec.witness) if spec.witness is not None else None
    )
    model, _ = build_row_model(spec)
    model_hash = _sha256_bytes(str(model.Proto()).encode("utf-8"))
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    status_name = solver.StatusName(status)
    proved = status == cp_model.INFEASIBLE
    _require(proved, f"{name}: exact one-worker model returned {status_name}")
    return {
        "spec": asdict(spec),
        "claim": (
            "row system infeasible"
            if spec.expected_energy is None
            else f"maximum energy exactly {spec.expected_energy}"
        ),
        "explicit_maximizer": witness,
        "proof_model": (
            "base constraints"
            if spec.expected_energy is None
            else f"base constraints plus energy>={spec.expected_energy + 1}"
        ),
        "model_proto_sha256": model_hash,
        "solver": "OR-Tools CP-SAT",
        "ortools_version": ORTOOLS_VERSION,
        "num_search_workers": 1,
        "random_seed": 0,
        "status": status_name,
        "proved": proved,
    }


def replay_all_row_certificates() -> dict[str, dict[str, object]]:
    """Replay independent one-worker models concurrently, preserving order."""
    # Warm the two immutable catalogs once.  Without this, simultaneous first
    # calls can redundantly enumerate the same middle slice in several threads.
    for p in (17, 19):
        translated_cut_vectors(p)
    names = tuple(spec.name for spec in ROW_SPECS)
    with ThreadPoolExecutor(
        max_workers=len(names), thread_name_prefix="prop15753-row"
    ) as executor:
        replays = tuple(executor.map(replay_row_certificate, names))
    rows = dict(zip(names, replays))
    _require(
        set(rows) == set(ROW_SPEC_BY_NAME)
        and all(row["proved"] for row in rows.values()),
        "not every exact one-worker row model replayed",
    )
    return rows


def p17_opposite_sign_regression() -> dict[str, object]:
    """Pin the hard/opposite quartic sign distinction at p=17."""
    spec = ROW_SPEC_BY_NAME["p17_A_opposite_Q4_quartic"]
    correct = spec.witness
    assert correct is not None
    old_wrong_sign_witness = (-2, -2, 0, -1, -2, 0, -1, -1)
    correct_m2, correct_m4 = moment_residues(17, correct)
    wrong_m2, wrong_m4 = moment_residues(17, old_wrong_sign_witness)
    proved = bool(
        correct_m4 == -(correct_m2 * correct_m2) % 17
        and correct_m4 != correct_m2 * correct_m2 % 17
        and wrong_m4 == wrong_m2 * wrong_m2 % 17
        and wrong_m4 != -(wrong_m2 * wrong_m2) % 17
        and sum(value * value for value in correct) == 11
        and sum(value * value for value in old_wrong_sign_witness) == 15
    )
    _require(proved, "the p17 opposite quartic sign regression failed")
    return {
        "global_quartic": "G=h*M4-M2^2",
        "normalization": "S_j=h*M_j on hard and S_j=-h*M_j on opposite",
        "hard_relation": "S4=S2^2 mod 17",
        "opposite_relation": "S4=-S2^2 mod 17",
        "correct_opposite_maximizer": list(correct),
        "correct_opposite_energy": 11,
        "wrong_plus_sign_witness": list(old_wrong_sign_witness),
        "wrong_plus_sign_energy": 15,
        "wrong_plus_sign_witness_rejected_by_correct_relation": True,
        "proved": proved,
    }


def endpoint_arithmetic(p: int) -> dict[str, object]:
    if p not in (17, 19):
        raise ValueError("this proposition closes only p=17 and p=19")
    t = 4
    k = 4 * p + 2 * t
    H_edges = k + 1
    q = (p - 1) // 2
    m = q + 1
    lift = sharp_integral_quadratic_lift_floor(p)
    floors = residual_even_floor_table(p)
    maximum_low_mean = p + 1 + 2 * t
    next_nonbaseline_phase_one_floor = 2 * p - 6 if p % 4 == 1 else 2 * p
    isolated_gap = p * p + 1 - 2 * H_edges
    proved = bool(
        H_edges == (77 if p == 17 else 85)
        and isolated_gap > 0
        and maximum_low_mean < next_nonbaseline_phase_one_floor
        and int(lift["sharp_scaled_floor"]) == p - 3
        and int(floors["least_nonzero_phase_zero_floor"]) == (16 if p == 17 else 20)
        and 2 * t + 2 < p - 3
        and t + 1 < m
    )
    _require(proved, f"p={p} fifth-shell arithmetic changed")
    return {
        "p": p,
        "layer_t": t,
        "original_k": k,
        "H_edge_count": H_edges,
        "ambient_vertex_count": p * p + 1,
        "guaranteed_isolated_vertices": isolated_gap,
        "q": q,
        "m": m,
        "distance_bin_count": q,
        "hard_direction_count": m,
        "opposite_direction_count": m,
        "type_budget": 2 * m * (m + t),
        "phase_one_mean_form": f"a_d=2u+{p + 1}*k_d",
        "phase_one_quotient_sum": "sum k_d=m+t-u",
        "maximum_low_phase_one_mean": maximum_low_mean,
        "next_nonbaseline_phase_one_floor": next_nonbaseline_phase_one_floor,
        "maximum_baseline_lift_excess": 2 * t + 2,
        "sharp_integral_lift_floor": int(lift["sharp_scaled_floor"]),
        "least_nonzero_phase_zero_floor": int(floors["least_nonzero_phase_zero_floor"]),
        "proved": proved,
    }


def hard_residue_ledger(p: int) -> dict[str, object]:
    """Derive, rather than assume, the exhaustive fifth-shell branches."""
    arithmetic = endpoint_arithmetic(p)
    t = int(arithmetic["layer_t"])
    q = int(arithmetic["q"])
    m = int(arithmetic["m"])
    phase_one = residual_even_floor_table(p)["phase_one_floors"]
    least_phase_one_floor = min(int(value) for value in phase_one.values())
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    low_rows: list[dict[str, object]] = []
    for u in range(t + 1):
        low_mean = p + 1 + 2 * u
        available = [
            (int(b), int(floor), low_mean - int(floor))
            for b, floor in phase_one.items()
            if int(floor) <= low_mean
        ]
        _require(
            all(
                b in {2, p - 1} and 0 <= excess < lift_floor
                for b, _, excess in available
            ),
            f"p={p}: a nonbaseline low phase-one cell entered the fifth shell",
        )
        exact = [(b, floor) for b, floor, excess in available if excess == 0]
        expected_exact = [(p - 1, p + 1)] if p == 17 and u == 0 else []
        _require(exact == expected_exact, f"p={p},u={u}: exact survivor list changed")
        k_zero_mean = 2 * u
        k_zero_forbidden = k_zero_mean < least_phase_one_floor
        minimum_k_after_lift_sieve = 1 if exact else 2
        quotient_sum = m + t - u
        excluded_by_sum = quotient_sum < m * minimum_k_after_lift_sieve
        low_rows.append(
            {
                "u": u,
                "low_mean": low_mean,
                "hypothetical_k_zero_mean": k_zero_mean,
                "least_phase_one_floor": least_phase_one_floor,
                "k_zero_forbidden": k_zero_forbidden,
                "available_endpoint_cells_b_floor_excess": [list(row) for row in available],
                "exact_survivors_b_floor": [list(row) for row in exact],
                "minimum_k_after_lift_sieve": minimum_k_after_lift_sieve,
                "quotient_sum": quotient_sum,
                "minimum_possible_quotient_sum": m * minimum_k_after_lift_sieve,
                "excluded_by_quotient_sum": excluded_by_sum,
                "surviving_branch": "B_LITERAL" if exact else None,
            }
        )

    endpoint_candidates = [
        int(value)
        for value in residual_even_floor_table(p)["phase_one_cells_at_mean_p_minus_one"]
    ]
    rules = baseline_coefficient_rules(p)
    if p == 17:
        possible = ["A_XNOR", "B_LITERAL"]
        expected_candidates = [2]
        endpoint_branches = {2: "A_XNOR"}
        no_mixing = True
        no_mixing_reason = "only b=2 attains the phase-one floor p-1"
    else:
        possible = ["A_XNOR", "C_COMPLEMENT_LITERAL"]
        expected_candidates = [2, p - 1]
        endpoint_branches = {2: "A_XNOR", p - 1: "C_COMPLEMENT_LITERAL"}
        b2_offset = int(rules[BRANCH_B2]["offset"])
        complement_offset = int(rules[BRANCH_P3_LAST]["offset"])
        no_mixing = (b2_offset - complement_offset) % q != 0
        no_mixing_reason = (
            "the common parallel count cannot satisfy both coefficient "
            f"congruences P={b2_offset} and P={complement_offset} mod {q}"
        )
    _require(endpoint_candidates == expected_candidates, f"p={p}: endpoint candidates changed")

    # For t<u<=m-2, sum k_d=m+t-u<m while every one of the m quotients is
    # at least one.  Thus only the low rows and u=m-1 endpoint remain.
    intermediate_u = list(range(t + 1, m - 1))
    intermediate_rows = [
        {
            "u": u,
            "hypothetical_k_zero_mean": 2 * u,
            "least_phase_one_floor": least_phase_one_floor,
            "k_zero_forbidden": 2 * u < least_phase_one_floor,
            "quotient_sum": m + t - u,
            "quotient_sum_less_than_direction_count": m + t - u < m,
        }
        for u in intermediate_u
    ]
    intermediate_impossible = all(
        row["k_zero_forbidden"]
        and row["quotient_sum_less_than_direction_count"]
        for row in intermediate_rows
    )
    endpoint_exact_count_lower = m - (t + 1)
    expected_p1_offset = int(rules[BRANCH_P1_LAST]["offset"])
    low_literal_offset_ok = p != 17 or expected_p1_offset == 5
    proved = bool(
        arithmetic["proved"]
        and intermediate_impossible
        and endpoint_exact_count_lower > 0
        and no_mixing
        and low_literal_offset_ok
        and all(
            row["k_zero_forbidden"]
            and (
                row["excluded_by_quotient_sum"]
                if row["surviving_branch"] is None
                else (p == 17 and row["u"] == 0)
            )
            for row in low_rows
        )
    )
    _require(proved, f"p={p}: hard-residue branch exhaustiveness failed")
    return {
        "p": p,
        "layer_t": t,
        "u_0_through_t_rows": low_rows,
        "intermediate_u_values": intermediate_u,
        "intermediate_u_rows": intermediate_rows,
        "intermediate_u_impossible_reason": (
            "k_d=0 would have mean 2u below the least phase-one floor; "
            "therefore all k_d>=1, but sum k_d=m+t-u<m"
        ),
        "u_equals_m_minus_1_exact_cell_count_at_least": endpoint_exact_count_lower,
        "u_equals_m_minus_1_endpoint_b_candidates": endpoint_candidates,
        "endpoint_candidate_branch_map": {str(key): value for key, value in endpoint_branches.items()},
        "equal_mean_endpoint_cells_cannot_mix": no_mixing,
        "no_mixing_reason": no_mixing_reason,
        "possible_branches": possible,
        "proved": proved,
    }


def _partitions(total: int) -> tuple[tuple[int, ...], ...]:
    def rec(remaining: int, maximum: int) -> list[tuple[int, ...]]:
        if remaining == 0:
            return [()]
        rows: list[tuple[int, ...]] = []
        for first in range(min(remaining, maximum), 0, -1):
            rows.extend((first,) + tail for tail in rec(remaining - first, first))
        return rows

    return tuple(sorted(rec(total, total), key=lambda row: (max(row), row)))


def _parseval_lower(
    *, p: int, H_edges: int, hT: int, parallel_squares: int, exact_energy: int
) -> int:
    return p * H_edges + 2 * hT * hT - 2 * parallel_squares - exact_energy


def branch_normalization_certificate(p: int, branch: str) -> dict[str, object]:
    """Glue local row sums to one global signed edge total.

    This is the load-bearing step that prevents a locally admissible row
    from choosing its own parallel-bin normalization.
    """
    if (p, branch) not in {
        (17, "A_XNOR"),
        (17, "B_LITERAL"),
        (19, "A_XNOR"),
        (19, "C_COMPLEMENT_LITERAL"),
    }:
        raise ValueError("unknown fifth-shell endpoint branch")
    arithmetic = endpoint_arithmetic(p)
    H_edges = int(arithmetic["H_edge_count"])
    q = int(arithmetic["q"])
    m = int(arithmetic["m"])
    rules = baseline_coefficient_rules(p)
    if branch == "A_XNOR":
        rule_name = BRANCH_B2
        opposite_delta = 4
        local_constant = p - 1
        expected_baseline_P = 4
        first_k = 0
    elif branch == "B_LITERAL":
        rule_name = BRANCH_P1_LAST
        opposite_delta = 4
        local_constant = 0
        expected_baseline_P = 5
        first_k = 1
    else:
        rule_name = BRANCH_P3_LAST
        opposite_delta = 5
        local_constant = p - 1
        expected_baseline_P = 3
        first_k = 0
    offset = int(rules[rule_name]["offset"])

    parameter_rows = []
    for parallel in range(H_edges + 1):
        numerator = parallel - offset
        if numerator % q:
            continue
        rho = numerator // q
        if rho < 0:
            continue
        s = parallel + rho
        opposite_edges = q * (8 - s) + opposite_delta
        if opposite_edges >= 0:
            parameter_rows.append((parallel, rho, s, opposite_edges))
    _require(
        len(parameter_rows) == 1
        and parameter_rows[0][0] == expected_baseline_P
        and parameter_rows[0][1] == 0,
        f"p={p} {branch}: coefficient congruence no longer fixes the edge split",
    )
    _, rho, s, opposite_edges = parameter_rows[0]
    hard_edges = H_edges - opposite_edges
    hT = hard_edges - opposite_edges

    # The unspecialized local row sum is
    #   p(P_L-3)-c-(p+1)k_L,
    # while every row of the one common graph has sum hT-P_L.
    # Solve this equality over the entire integral P range for several k;
    # it has the asserted affine solution and no alternative normalization.
    normalization_rows = []
    for k in range(first_k, first_k + 6):
        matches = []
        for parallel in range(H_edges + 1):
            local_sum = p * (parallel - 3) - local_constant - (p + 1) * k
            common_sum = hT - parallel
            if local_sum == common_sum:
                matches.append(parallel)
        expected = (expected_baseline_P - first_k) + k
        normalization_rows.append(
            {
                "k": k,
                "matching_parallel_counts_in_full_range": matches,
                "forced_parallel_count": expected,
                "local_sum": p * (expected - 3) - local_constant - (p + 1) * k,
                "common_sum": hT - expected,
            }
        )
    affine_intercept = expected_baseline_P - first_k
    if branch == "A_XNOR":
        exact_baseline_row = (1,) + (0,) * (q - 1)
    elif branch == "B_LITERAL":
        exact_baseline_row = (2,) * q
    else:
        exact_baseline_row = (-2,) * q
    exact_baseline_sum = sum(exact_baseline_row)
    exact_baseline_energy = sum(value * value for value in exact_baseline_row)
    proved = bool(
        rho == 0
        and s == expected_baseline_P
        and hard_edges + opposite_edges == H_edges
        and hT == hard_edges - opposite_edges
        and exact_baseline_sum == hT - expected_baseline_P
        and all(
            row["matching_parallel_counts_in_full_range"]
            == [row["forced_parallel_count"]]
            and row["local_sum"] == row["common_sum"]
            for row in normalization_rows
        )
    )
    _require(proved, f"p={p} {branch}: common normalization glue failed")
    return {
        "p": p,
        "branch": branch,
        "baseline_rule": rule_name,
        "coefficient_offset": offset,
        "feasible_baseline_P_rho_s_opposite_edges": [list(row) for row in parameter_rows],
        "hard_edge_count": hard_edges,
        "opposite_edge_count": opposite_edges,
        "hT_from_edge_split": hT,
        "local_unspecialized_sum_formula": (
            f"{p}*(P_L-3)-{local_constant}-{p + 1}*k_L"
        ),
        "common_difference_sum_formula": "hT-P_L",
        "opposite_difference_sum_formula": "-hT-Q_L",
        "exact_baseline_row": list(exact_baseline_row),
        "exact_baseline_row_sum": exact_baseline_sum,
        "exact_baseline_row_energy": exact_baseline_energy,
        "exact_baseline_agrees_with_independent_edge_split": True,
        "normalization_rows": normalization_rows,
        "forced_affine_normalization": f"P_L={affine_intercept}+k_L",
        "proved": proved,
    }


def p17_branch_A_certificate(
    row_replays: dict[str, dict[str, object]],
) -> dict[str, object]:
    normalization = branch_normalization_certificate(17, "A_XNOR")
    floor_gate = endpoint_arithmetic(17)
    opposite_q3_excluded = bool(
        8 < int(floor_gate["least_nonzero_phase_zero_floor"])
        and 8 < int(floor_gate["sharp_integral_lift_floor"])
    )
    raw_hard = {excess: ROW_SPEC_BY_NAME[f"p17_A_hard_e{excess}_raw"].expected_energy for excess in range(1, 5)}
    raw_opposite = ROW_SPEC_BY_NAME["p17_A_opposite_Q4_raw"].expected_energy
    quartic_hard = ROW_SPEC_BY_NAME["p17_A_hard_e5_quartic"].expected_energy
    quartic_opposite = ROW_SPEC_BY_NAME["p17_A_opposite_Q4_quartic"].expected_energy
    assert all(value is not None for value in raw_hard.values())
    assert raw_opposite is not None and quartic_hard is not None and quartic_opposite is not None
    rows: list[dict[str, object]] = []
    for partition in _partitions(5):
        exact_count = 9 - len(partition)
        hard_parallel = [4 + excess for excess in partition] + [4] * exact_count
        parallel_squares = sum(value * value for value in hard_parallel) + 9 * 4 * 4
        lower = _parseval_lower(
            p=17,
            H_edges=77,
            hT=5,
            parallel_squares=parallel_squares,
            exact_energy=exact_count,
        )
        if partition == (5,):
            upper = quartic_hard + 9 * quartic_opposite
            mechanism = "eight exact XNOR roots force G identically zero"
        else:
            upper = sum(int(raw_hard[value]) for value in partition) + 9 * raw_opposite
            mechanism = "raw translated-cut row maxima"
        rows.append(
            {
                "excess_partition": list(partition),
                "hard_parallel_counts": sorted(hard_parallel),
                "opposite_parallel_counts": [4] * 9,
                "parallel_square_sum": parallel_squares,
                "exact_XNOR_row_count": exact_count,
                "exact_XNOR_row_energy_each": 1,
                "nonexact_Parseval_lower": lower,
                "collision_increment": "34*C, C>=0",
                "row_energy_upper": upper,
                "strict_gap_at_C_zero": lower - upper,
                "mechanism": mechanism,
                "excluded": lower > upper,
            }
        )
    expected_gaps = [342, 312, 282, 212, 182, 138, 162]
    proved = bool(
        len(rows) == 7
        and [row["strict_gap_at_C_zero"] for row in rows] == expected_gaps
        and all(row["excluded"] for row in rows)
        and p17_opposite_sign_regression()["proved"]
        and normalization["hT_from_edge_split"] == 5
        and normalization["forced_affine_normalization"] == "P_L=4+k_L"
        and normalization["exact_baseline_row_energy"] == 1
        and opposite_q3_excluded
        and all(
            row_replays[name]["proved"]
            for name in (
                "p17_A_hard_e1_raw",
                "p17_A_hard_e2_raw",
                "p17_A_hard_e3_raw",
                "p17_A_hard_e4_raw",
                "p17_A_opposite_Q4_raw",
                "p17_A_hard_e5_quartic",
                "p17_A_opposite_Q4_quartic",
            )
        )
    )
    _require(proved, "p17 branch A ledger failed")
    return {
        "branch": "A_XNOR",
        "common_normalization": normalization,
        "hard_edge_count": 41,
        "opposite_edge_count": 36,
        "hT": 5,
        "hard_parallel_normalization": "P_L=4+k_L",
        "hard_excess_sum": 5,
        "opposite_Q3_mean_8_exclusion": {
            "nonzero_b_floor": int(floor_gate["least_nonzero_phase_zero_floor"]),
            "b_zero_integral_lift_floor": int(floor_gate["sharp_integral_lift_floor"]),
            "proved": opposite_q3_excluded,
        },
        "opposite_parallel_counts": [4] * 9,
        "difference_Radon_Parseval": "17*77+2*5^2-2*sum_L(P_L^2)+34*C",
        "quartic_root_argument": {
            "last_partition": [5],
            "exact_hard_XNOR_roots": 8,
            "binary_quartic_degree": 4,
            "root_count_forces_zero_polynomial": True,
            "sign_regression": p17_opposite_sign_regression(),
        },
        "partition_ledgers": rows,
        "excluded": proved,
        "proved": proved,
    }


def p17_branch_B_certificate(
    row_replays: dict[str, dict[str, object]],
) -> dict[str, object]:
    normalization = branch_normalization_certificate(17, "B_LITERAL")
    h1 = ROW_SPEC_BY_NAME["p17_B_hard_e1_zero"]
    h2 = ROW_SPEC_BY_NAME["p17_B_hard_e2_zero"].expected_energy
    h4 = ROW_SPEC_BY_NAME["p17_B_hard_e4_zero"].expected_energy
    o3 = ROW_SPEC_BY_NAME["p17_B_opposite_Q3_zero"].expected_energy
    o4 = ROW_SPEC_BY_NAME["p17_B_opposite_Q4_zero"].expected_energy
    assert h1.expected_energy is None
    assert h2 is not None and h4 is not None and o3 is not None and o4 is not None
    ledgers = []
    for partition, exact_count, hard_parallel, hard_upper in (
        ((2, 2), 7, [5] * 7 + [7] * 2, 2 * h2),
        ((4,), 8, [5] * 8 + [9], h4),
    ):
        opposite_parallel = [3] * 8 + [4]
        parallel_squares = sum(value * value for value in hard_parallel + opposite_parallel)
        lower = _parseval_lower(
            p=17,
            H_edges=77,
            hT=21,
            parallel_squares=parallel_squares,
            exact_energy=exact_count * 32,
        )
        upper = hard_upper + 8 * o3 + o4
        ledgers.append(
            {
                "excess_partition": list(partition),
                "hard_parallel_counts": sorted(hard_parallel),
                "opposite_parallel_counts": sorted(opposite_parallel),
                "parallel_square_sum": parallel_squares,
                "exact_literal_star_count": exact_count,
                "exact_literal_star_energy_each": 32,
                "nonexact_Parseval_lower": lower,
                "collision_increment": "34*C, C>=0",
                "row_energy_upper": upper,
                "strict_gap_at_C_zero": lower - upper,
                "excluded": lower > upper,
            }
        )
    proved = bool(
        h1.expected_energy is None
        and [row["strict_gap_at_C_zero"] for row in ledgers] == [428, 302]
        and all(row["excluded"] for row in ledgers)
        and normalization["hT_from_edge_split"] == 21
        and normalization["forced_affine_normalization"] == "P_L=4+k_L"
        and normalization["exact_baseline_row_energy"] == 32
        and all(
            row_replays[name]["proved"]
            for name in (
                "p17_B_hard_e1_zero",
                "p17_B_hard_e2_zero",
                "p17_B_hard_e4_zero",
                "p17_B_opposite_Q3_zero",
                "p17_B_opposite_Q4_zero",
            )
        )
    )
    _require(proved, "p17 branch B ledger failed")
    return {
        "branch": "B_LITERAL",
        "common_normalization": normalization,
        "hard_edge_count": 49,
        "opposite_edge_count": 28,
        "hT": 21,
        "hard_parallel_normalization": "P_L=4+k_L",
        "hard_excess_sum": 4,
        "opposite_parallel_counts": [3] * 8 + [4],
        "exact_literal_stars_at_least": 5,
        "global_moments_forced": ["M2=0", "M4=0"],
        "partitions_containing_excess_1_excluded_by_empty_row_model": True,
        "remaining_partition_ledgers": ledgers,
        "excluded": proved,
        "proved": proved,
    }


def p19_branch_A_certificate(
    row_replays: dict[str, dict[str, object]],
) -> dict[str, object]:
    normalization = branch_normalization_certificate(19, "A_XNOR")
    floor_gate = endpoint_arithmetic(19)
    opposite_q3_excluded = bool(
        8 < int(floor_gate["least_nonzero_phase_zero_floor"])
        and 8 < int(floor_gate["sharp_integral_lift_floor"])
    )
    hard = {
        excess: int(ROW_SPEC_BY_NAME[f"p19_A_hard_e{excess}_raw"].expected_energy)  # type: ignore[arg-type]
        for excess in range(1, 6)
    }
    opposite = int(ROW_SPEC_BY_NAME["p19_A_opposite_Q4_raw"].expected_energy)  # type: ignore[arg-type]
    rows: list[dict[str, object]] = []
    for partition in _partitions(5):
        exact_count = 10 - len(partition)
        hard_parallel = [4 + excess for excess in partition] + [4] * exact_count
        parallel_squares = sum(value * value for value in hard_parallel) + 10 * 4 * 4
        lower = _parseval_lower(
            p=19,
            H_edges=85,
            hT=5,
            parallel_squares=parallel_squares,
            exact_energy=exact_count,
        )
        upper = sum(hard[value] for value in partition) + 10 * opposite
        rows.append(
            {
                "excess_partition": list(partition),
                "hard_parallel_counts": sorted(hard_parallel),
                "opposite_parallel_counts": [4] * 10,
                "parallel_square_sum": parallel_squares,
                "exact_XNOR_row_count": exact_count,
                "exact_XNOR_row_energy_each": 1,
                "nonexact_Parseval_lower": lower,
                "collision_increment": "38*C, C>=0",
                "row_energy_upper": upper,
                "strict_gap_at_C_zero": lower - upper,
                "excluded": lower > upper,
            }
        )
    expected_gaps = [520, 490, 460, 420, 390, 312, 162]
    proved = bool(
        [row["strict_gap_at_C_zero"] for row in rows] == expected_gaps
        and all(row["excluded"] for row in rows)
        and normalization["hT_from_edge_split"] == 5
        and normalization["forced_affine_normalization"] == "P_L=4+k_L"
        and normalization["exact_baseline_row_energy"] == 1
        and opposite_q3_excluded
        and all(
            row_replays[name]["proved"]
            for name in (
                "p19_A_hard_e1_raw",
                "p19_A_hard_e2_raw",
                "p19_A_hard_e3_raw",
                "p19_A_hard_e4_raw",
                "p19_A_hard_e5_raw",
                "p19_A_opposite_Q4_raw",
            )
        )
    )
    _require(proved, "p19 branch A ledger failed")
    return {
        "branch": "A_XNOR",
        "common_normalization": normalization,
        "hard_edge_count": 45,
        "opposite_edge_count": 40,
        "hT": 5,
        "hard_parallel_normalization": "P_L=4+k_L",
        "hard_excess_sum": 5,
        "opposite_Q3_mean_8_exclusion": {
            "nonzero_b_floor": int(floor_gate["least_nonzero_phase_zero_floor"]),
            "b_zero_integral_lift_floor": int(floor_gate["sharp_integral_lift_floor"]),
            "proved": opposite_q3_excluded,
        },
        "opposite_parallel_counts": [4] * 10,
        "difference_Radon_Parseval": "19*85+2*5^2-2*sum_L(P_L^2)+38*C",
        "partition_ledgers": rows,
        "excluded": proved,
        "proved": proved,
    }


def p19_branch_C_certificate(
    row_replays: dict[str, dict[str, object]],
) -> dict[str, object]:
    normalization = branch_normalization_certificate(19, "C_COMPLEMENT_LITERAL")
    spec = ROW_SPEC_BY_NAME["p19_C_opposite_Q5_zero"]
    floor_gate = endpoint_arithmetic(19)
    opposite_q4_excluded = bool(
        8 < int(floor_gate["least_nonzero_phase_zero_floor"])
        and 8 < int(floor_gate["sharp_integral_lift_floor"])
    )
    proved = bool(
        spec.expected_energy is None
        and spec.moment_relation == "zero"
        and normalization["hT_from_edge_split"] == -15
        and normalization["forced_affine_normalization"] == "P_L=3+k_L"
        and normalization["exact_baseline_row_energy"] == 36
        and opposite_q4_excluded
        and row_replays[spec.name]["proved"]
    )
    _require(proved, "p19 branch C arithmetic failed")
    return {
        "branch": "C_COMPLEMENT_LITERAL",
        "common_normalization": normalization,
        "hard_edge_count": 35,
        "opposite_edge_count": 50,
        "hT": -15,
        "hard_parallel_normalization": "P_L=3+k_L",
        "hard_excess_sum": 5,
        "opposite_Q4_mean_8_exclusion": {
            "nonzero_b_floor": int(floor_gate["least_nonzero_phase_zero_floor"]),
            "b_zero_integral_lift_floor": int(floor_gate["sharp_integral_lift_floor"]),
            "proved": opposite_q4_excluded,
        },
        "opposite_parallel_counts": [5] * 10,
        "exact_complement_literal_stars_at_least": 5,
        "global_moments_forced": ["M2=0", "M4=0"],
        "opposite_Q5_row_constraints": {
            "sum": spec.total,
            "l1_bound": spec.l1_bound,
            "translated_cut_upper": spec.cut_upper,
            "translated_cut_count": EXPECTED_CUT_CATALOGS[19][0],
            "moment_relation": spec.moment_relation,
        },
        "opposite_Q5_row_system_infeasible": row_replays[spec.name]["proved"],
        "excluded": proved,
        "proved": proved,
    }


def branch_ledger_certificate() -> dict[str, object]:
    row_replays = replay_all_row_certificates()
    _require(
        set(row_replays) == set(ROW_SPEC_BY_NAME)
        and all(row["proved"] for row in row_replays.values()),
        "the branch ledger requires all nineteen live row-model replays",
    )
    p17_residues = hard_residue_ledger(17)
    p19_residues = hard_residue_ledger(19)
    p17_A = p17_branch_A_certificate(row_replays)
    p17_B = p17_branch_B_certificate(row_replays)
    p19_A = p19_branch_A_certificate(row_replays)
    p19_C = p19_branch_C_certificate(row_replays)
    proved = bool(
        p17_residues["possible_branches"] == ["A_XNOR", "B_LITERAL"]
        and p19_residues["possible_branches"] == ["A_XNOR", "C_COMPLEMENT_LITERAL"]
        and p17_A["proved"]
        and p17_B["proved"]
        and p19_A["proved"]
        and p19_C["proved"]
    )
    _require(proved, "a fifth-shell endpoint branch remained")
    return {
        "p17": {
            "hard_residue_exhaustiveness": p17_residues,
            "possible_branches": p17_residues["possible_branches"],
            "branch_A": p17_A,
            "branch_B": p17_B,
            "p17_k76_closed": True,
        },
        "p19": {
            "hard_residue_exhaustiveness": p19_residues,
            "possible_branches": p19_residues["possible_branches"],
            "branch_A": p19_A,
            "branch_C": p19_C,
            "p19_k84_closed": True,
        },
        "proved": proved,
    }


def _supporting_artifact_hashes() -> dict[str, str]:
    paths = (
        ROOT / "src" / "e1_gmin_m4_prop15753.py",
        ROOT / "tests" / "test_prop15753.py",
        ROOT / "evidence" / "NOTE_2026-09-01_P17_P19_FIFTH_SHELL_CLOSE.md",
    )
    return {
        str(path.relative_to(ROOT)): _sha256_bytes(path.read_bytes()) for path in paths
    }


def proposition_15753() -> dict[str, object]:
    cuts = {str(p): cut_catalog_certificate(p) for p in (17, 19)}
    arithmetic = {str(p): endpoint_arithmetic(p) for p in (17, 19)}
    ledgers = branch_ledger_certificate()
    rows = replay_all_row_certificates()
    model_hashes = {spec.name: row_model_sha256(spec) for spec in ROW_SPECS}
    spec_manifest = [asdict(spec) for spec in ROW_SPECS]
    certificate_manifest = {
        "cut_catalog_sha256": {
            p: row["catalog_sha256"] for p, row in cuts.items()
        },
        "row_spec_sha256": _canonical_sha256(spec_manifest),
        "row_model_proto_sha256": model_hashes,
        "supporting_artifact_sha256": _supporting_artifact_hashes(),
    }
    global_gate = residual_ii_k_ge_4p_ND_closed()
    _require(global_gate is False, "a finite endpoint certificate must not flip the global gate")
    proved = bool(
        ledgers["proved"]
        and all(row["proved"] for row in cuts.values())
        and all(row["proved"] for row in arithmetic.values())
        and set(rows) == set(ROW_SPEC_BY_NAME)
        and all(row["proved"] for row in rows.values())
        and not global_gate
    )
    _require(proved, "Proposition 15.753 certificate failed")
    return {
        "proposition": "15.753",
        "title": "translated-cut energy closes the p=17,19 fifth-shell endpoints",
        "result_class": "exhaustive finite aggregate certificate and proved endpoint theorem",
        "changed_premise": (
            "Proposition 15.752 closes the fifth shell for p>=23; this exact "
            "certificate handles its sharp exceptional endpoints p=17,19"
        ),
        "endpoint_arithmetic": arithmetic,
        "translated_cut_catalogs": cuts,
        "branch_ledgers": ledgers,
        "exact_one_worker_row_replays": rows,
        "all_row_models_replayed": True,
        "independent_row_models_run_concurrently": True,
        "maximum_concurrent_row_models": len(ROW_SPECS),
        "certificate_manifest": certificate_manifest,
        "certificate_manifest_sha256": _canonical_sha256(certificate_manifest),
        "p17_k76_closed": True,
        "p19_k84_closed": True,
        "finite_prime_aggregate_census_used": True,
        "graph_or_configuration_census_used": False,
        "residual_ii_k_ge_4p_ND_closed": global_gate,
        "E1_closed": False,
        "quadratic_minmax_limit_closed": False,
        "proved": proved,
    }


def write_evidence(path: Path = EVIDENCE_PATH) -> Path:
    write_json_atomic(path, proposition_15753())
    return path


def main() -> None:
    path = write_evidence()
    print(path)


if __name__ == "__main__":
    main()
