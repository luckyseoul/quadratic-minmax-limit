#!/usr/bin/env python3
r"""Six-bin common-energy close of the high-root p=13,t=4,u=6 strata.

This module treats only the excess partitions ``(3,2)``, ``(4,1)``, and
``(5)``.  Their five or six exact XNOR directions force

    G(L) = h M_4(L) - M_2(L)^2

to vanish identically.  The remaining rows are bounded by all 74 pinned
translated cuts and by the global displacement-collision budget.  No graph,
coefficient-cell, orbit, or common-realization census is used.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable

from ortools import __version__ as ORTOOLS_VERSION
from ortools.sat.python import cp_model

from e1_gmin_m4_prop15740 import translated_cut_vector_catalog
from e1_gmin_m4_prop15745 import finite_field_bucket_sign_audit


P = 13
DISTANCES = tuple(range(1, 7))
H_EDGE_COUNT = 61
HARD_SIGN_TIMES_T = 5
OPPOSITE_DIRECTION_COUNT = 7
OPPOSITE_PARALLEL = 4


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _proto_sha256(model: cp_model.CpModel) -> str:
    return hashlib.sha256(str(model.Proto()).encode("utf-8")).hexdigest()


@lru_cache(maxsize=1)
def translated_cuts() -> tuple[tuple[int, ...], ...]:
    payload = translated_cut_vector_catalog()
    rows = tuple(tuple(int(value) for value in row) for row in payload["vectors"])
    _require(
        len(rows) == 74
        and all(len(row) == 6 and sum(row) == 42 for row in rows),
        "the pinned p=13 translated-cut catalog changed",
    )
    return rows


def balanced_collision(total: int, boxes: int) -> int:
    """Minimum ``sum binom(m_i,2)`` for ``total`` objects in ``boxes``."""
    quotient, remainder = divmod(int(total), int(boxes))
    return (
        remainder * quotient * (quotient + 1) // 2
        + (boxes - remainder) * quotient * (quotient - 1) // 2
    )


def signed_collision_floor(value: int) -> int:
    """Minimum collisions needed for one normalized signed bucket value."""
    value = int(value)
    return (
        balanced_collision(value, 6)
        if value >= 0
        else balanced_collision(-value, 7)
    )


@lru_cache(maxsize=1)
def collision_floor_certificate() -> dict[str, object]:
    incidence = finite_field_bucket_sign_audit()
    sample_rows = []
    for value in range(-20, 21):
        floor = signed_collision_floor(value)
        brute = min(
            balanced_collision(positive_total, 6)
            + balanced_collision(positive_total - value, 7)
            for positive_total in range(max(value, 0), 41)
            if 0 <= positive_total - value <= 40
        )
        _require(floor == brute, f"signed collision floor failed at {value}")
        sample_rows.append((value, floor))
    proved = bool(
        incidence["proved"]
        and incidence["normalized_hard_sign_classes_per_bucket"] == 6
        and incidence["normalized_opposite_sign_classes_per_bucket"] == 7
    )
    _require(proved, "the p=13 signed bucket incidence changed")
    return {
        "positive_classes_per_nonzero_bucket": 6,
        "negative_classes_per_nonzero_bucket": 7,
        "checked_values": [list(row) for row in sample_rows],
        "inequality": (
            "sum_a signed_collision_floor(q_a) plus the zero-bin floor "
            "is at most the global collision C"
        ),
        "proved": proved,
    }


@dataclass(frozen=True)
class RowSpec:
    name: str
    kind: str
    excess: int | None
    collision_upper: int | None
    expected_energy: int
    witness: tuple[int, ...]

    @property
    def total(self) -> int:
        return 1 - int(self.excess) if self.kind == "hard" else -9

    @property
    def l1_bound(self) -> int:
        return 57 - int(self.excess) if self.kind == "hard" else 57

    @property
    def cut_upper(self) -> int:
        return 13 if self.kind == "hard" else -52

    @property
    def zero_bin_collision_floor(self) -> int:
        if self.kind == "opposite":
            return 0
        return balanced_collision(4 + int(self.excess), 6)

    @property
    def moment_relation(self) -> str:
        return "hard_quartic" if self.kind == "hard" else "opposite_quartic"


ROW_SPECS = (
    RowSpec("hard_e1_quartic", "hard", 1, None, 8, (0, 0, 2, -2, 0, 0)),
    RowSpec("hard_e2_quartic", "hard", 2, None, 51, (3, -1, -4, -3, 0, 4)),
    RowSpec("hard_e3_quartic_C1", "hard", 3, 1, 104, (4, -1, -6, -5, 1, 5)),
    RowSpec("opposite_q4_quartic", "opposite", None, None, 19, (-1, -3, -2, -1, 0, -2)),
    RowSpec("hard_e4_quartic_raw", "hard", 4, None, 229, (-7, 1, 5, -8, 9, -3)),
    RowSpec("hard_e4_quartic_C2", "hard", 4, 2, 159, (-6, 2, 5, -7, 6, -3)),
    RowSpec("hard_e4_quartic_C3", "hard", 4, 3, 163, (-5, 2, 5, -8, 6, -3)),
    RowSpec("hard_e5_quartic_raw", "hard", 5, None, 310, (-10, 10, -8, 6, -3, 1)),
    RowSpec("hard_e5_quartic_C3", "hard", 5, 3, 166, (5, 5, 3, -3, -7, -7)),
    RowSpec("hard_e5_quartic_C4", "hard", 5, 4, 188, (6, -2, -7, -7, -1, 7)),
    RowSpec("hard_e5_quartic_C5", "hard", 5, 5, 200, (-7, -1, 7, -8, 6, -1)),
    RowSpec("hard_e5_quartic_C6", "hard", 5, 6, 224, (-7, -2, 7, -7, 8, -3)),
    RowSpec("hard_e5_quartic_C7", "hard", 5, 7, 226, (0, -9, -3, 8, 6, -6)),
)
ROW_SPEC_BY_NAME = {spec.name: spec for spec in ROW_SPECS}

EXPECTED_MODEL_HASHES = {
    "hard_e1_quartic": ("9313c66a89a832f442abe3fce840eb9a8523fe9df3c94196fbba94e81b3169e7", "068bb75320071401a0f65018902f0a815be8b5e5f88fd900f594c2b1d52f1421"),
    "hard_e2_quartic": ("69bc08c9d47985803e4f24e5a14166ab57ef907f4fc1b1780c60bbdd99a08e09", "2d1da854399ff3b940dad82cd95a5ad223bde3807dab88240090c852e83c383c"),
    "hard_e3_quartic_C1": ("4d72383a807c84c59c052659e847b39624d933fb688bc0544ba4c73b568111fb", "654e4f74f914c0ecbe8b0d94b705d0ff70322b3b720486da5a4db87cfa38e322"),
    "opposite_q4_quartic": ("c4cbd6bf46f504d02c44ade1fd879f7e9d518e730883e480a6cb9278367daadb", "f728775847e07e543a6e39dd8eb7eca690cb1629be53b0b49185ff071b270cff"),
    "hard_e4_quartic_raw": ("e1182942d157f74973a31a6a5eecb80bbe46663b8e6d66f971fa604e291cb795", "041b9ad9e47ba36bfbb0c84cb88531b84858b4b1be997d0cc04b6e10b7dbb0e0"),
    "hard_e4_quartic_C2": ("6be278fe4a387917e8b1ca251fea0b84cddd106d9c5d98c0705e196957a09966", "ea36467ef0733e954cda62ca32cb26ad671ed95fdf2ebfed674815e6f3cf191c"),
    "hard_e4_quartic_C3": ("ca0b2786e68b68a97b97c8f89faa678a68ea056d2f63959a59eda3e2e210a148", "4ea22c6c722c7897b4119261c39cf6b81ea49ab8935dba76a7a67aa1ba450c4b"),
    "hard_e5_quartic_raw": ("6311109100f1e381c57a077434fed12288d5d4ac20c7993a0d67326d3cc79f0a", "568595418e47ec02fa172410a88e5febbc5e5a1b777cc576fdd3a236f026cd4f"),
    "hard_e5_quartic_C3": ("6b4fec0efec99ef4de9bf56ab18f94a19e5803ce329a72005bc5f9c30be4e617", "b5b06f10f15c6b10844ca461fcea2c59cdfdde68f92ed0b08b97bd6811fafe2a"),
    "hard_e5_quartic_C4": ("277f6510efe875a09e07bac3bc0b4873f60d194e124bdcb39f22e1e78da6a0b8", "062d000a0aad4ef331cb6912c85fba7980b0463a3bbd8840d202401593f28fe7"),
    "hard_e5_quartic_C5": ("506708abbc5e25124f33cf65398158d14868343dc76d0deb064158f7e8d5de0e", "4aac4898ae13647cf5c0836a7004d646e345d0e230d79c90b7005e92f7acbd99"),
    "hard_e5_quartic_C6": ("10ec598f3e5dd33fe7b3badaa0a4042e909db18653f8af753c388a5e2ee4aa3f", "5859b397d6139926dca610b0d353966ee43c99865eab5c16d9909dbcee41cb95"),
    "hard_e5_quartic_C7": ("7366a1f98a82f37d94a798386eb322b46006bbd04e4a0885e5aef677c516d58a", "133a6e1e7e2087eec4a575d6079b2c757045a746a207b71fcd861b3637ae987c"),
}


def moment_residues(row: Iterable[int]) -> tuple[int, int]:
    values = tuple(int(value) for value in row)
    return tuple(
        sum(pow(distance, degree, P) * value for distance, value in zip(DISTANCES, values)) % P
        for degree in (2, 4)
    )  # type: ignore[return-value]


def moment_relation_holds(spec: RowSpec, row: Iterable[int]) -> bool:
    m2, m4 = moment_residues(row)
    if spec.kind == "hard":
        return m4 == m2 * m2 % P
    return m4 == -(m2 * m2) % P


def row_satisfies(spec: RowSpec, row: Iterable[int]) -> bool:
    values = tuple(int(value) for value in row)
    collision_ok = bool(
        spec.collision_upper is None
        or sum(signed_collision_floor(value) for value in values)
        + spec.zero_bin_collision_floor
        <= spec.collision_upper
    )
    return bool(
        len(values) == 6
        and sum(values) == spec.total
        and sum(abs(value) for value in values) <= spec.l1_bound
        and all(
            sum(coefficient * value for coefficient, value in zip(cut, values))
            <= spec.cut_upper
            for cut in translated_cuts()
        )
        and collision_ok
        and moment_relation_holds(spec, values)
    )


def _build_row_model(
    spec: RowSpec,
    *,
    threshold: int | None,
    table_encoding: bool,
) -> tuple[cp_model.CpModel, list[cp_model.IntVar], cp_model.LinearExpr]:
    lower, upper = -spec.l1_bound, spec.l1_bound
    model = cp_model.CpModel()
    values = [model.NewIntVar(lower, upper, f"q{distance}") for distance in DISTANCES]
    absolute = [model.NewIntVar(0, spec.l1_bound, f"a{distance}") for distance in DISTANCES]
    squares = [
        model.NewIntVar(0, spec.l1_bound * spec.l1_bound, f"s{distance}")
        for distance in DISTANCES
    ]
    collisions = [
        model.NewIntVar(0, balanced_collision(spec.l1_bound, 6), f"c{distance}")
        for distance in DISTANCES
    ]
    full_table = [
        (value, abs(value), value * value, signed_collision_floor(value))
        for value in range(lower, upper + 1)
    ]
    collision_table = [(value, signed_collision_floor(value)) for value in range(lower, upper + 1)]
    for value, magnitude, square, collision in zip(values, absolute, squares, collisions):
        if table_encoding:
            model.AddAllowedAssignments([value, magnitude, square, collision], full_table)
        else:
            model.AddAbsEquality(magnitude, value)
            model.AddMultiplicationEquality(square, [value, value])
            model.AddAllowedAssignments([value, collision], collision_table)

    model.Add(sum(values) == spec.total)
    model.Add(sum(absolute) <= spec.l1_bound)
    if spec.collision_upper is not None:
        model.Add(
            sum(collisions) + spec.zero_bin_collision_floor
            <= spec.collision_upper
        )
    for cut in translated_cuts():
        model.Add(
            sum(coefficient * value for coefficient, value in zip(cut, values))
            <= spec.cut_upper
        )

    residues = []
    for degree in (2, 4):
        residue = model.NewIntVar(0, P - 1, f"M{degree}")
        quotient = model.NewIntVar(-1000, 1000, f"z{degree}")
        model.Add(
            sum(
                pow(distance, degree, P) * value
                for distance, value in zip(DISTANCES, values)
            )
            == P * quotient + residue
        )
        residues.append(residue)
    m2, m4 = residues
    allowed = [
        (
            residue,
            (
                residue * residue
                if spec.kind == "hard"
                else -(residue * residue)
            )
            % P,
        )
        for residue in range(P)
    ]
    model.AddAllowedAssignments([m2, m4], allowed)

    energy = sum(squares)
    if threshold is None:
        model.Maximize(energy)
    else:
        model.Add(energy >= threshold)
    return model, values, energy


def _solve(model: cp_model.CpModel) -> tuple[cp_model.CpSolver, int]:
    validation = model.Validate()
    _require(not validation, f"invalid CP-SAT model: {validation}")
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    return solver, status


@lru_cache(maxsize=None)
def row_energy_certificate(name: str) -> dict[str, object]:
    try:
        spec = ROW_SPEC_BY_NAME[name]
    except KeyError as exc:
        raise ValueError(f"unknown p13 u6 row spec {name!r}") from exc
    _require(
        row_satisfies(spec, spec.witness)
        and sum(value * value for value in spec.witness) == spec.expected_energy,
        f"{name}: explicit maximizer failed",
    )

    optimization, values, energy = _build_row_model(
        spec, threshold=None, table_encoding=True
    )
    optimization_solver, optimization_status = _solve(optimization)
    _require(
        optimization_status == cp_model.OPTIMAL
        and int(round(optimization_solver.ObjectiveValue())) == spec.expected_energy,
        f"{name}: sharp maximum changed",
    )
    replay, _values, _energy = _build_row_model(
        spec, threshold=spec.expected_energy + 1, table_encoding=False
    )
    replay_solver, replay_status = _solve(replay)
    _require(replay_status == cp_model.INFEASIBLE, f"{name}: upper replay became feasible")
    optimization_hash = _proto_sha256(optimization)
    replay_hash = _proto_sha256(replay)
    _require(
        (optimization_hash, replay_hash) == EXPECTED_MODEL_HASHES[name],
        f"{name}: pinned model hash changed",
    )
    return {
        "name": name,
        "kind": spec.kind,
        "excess": spec.excess,
        "collision_upper": spec.collision_upper,
        "zero_bin_collision_floor": spec.zero_bin_collision_floor,
        "moment_relation": spec.moment_relation,
        "sharp_energy": spec.expected_energy,
        "explicit_maximizer": list(spec.witness),
        "optimization": {
            "status": optimization_solver.StatusName(optimization_status),
            "workers": 1,
            "model_proto_sha256": optimization_hash,
            "solver_version": ORTOOLS_VERSION,
        },
        "upper_replay": {
            "status": replay_solver.StatusName(replay_status),
            "workers": 1,
            "forbidden_energy_floor": spec.expected_energy + 1,
            "model_proto_sha256": replay_hash,
            "solver_version": ORTOOLS_VERSION,
        },
        "proved": True,
    }


def nonexact_parseval_base(partition: tuple[int, ...]) -> int:
    exact_count = 7 - len(partition)
    hard_parallel = [4 + excess for excess in partition] + [4] * exact_count
    all_parallel_squares = sum(value * value for value in hard_parallel) + 7 * 4 * 4
    return (
        P * H_EDGE_COUNT
        + 2 * HARD_SIGN_TIMES_T * HARD_SIGN_TIMES_T
        - 2 * all_parallel_squares
        - exact_count
    )


@lru_cache(maxsize=1)
def high_root_partition_certificate() -> dict[str, object]:
    collision_floor_certificate()
    rows = {name: row_energy_certificate(name) for name in ROW_SPEC_BY_NAME}
    h1 = int(rows["hard_e1_quartic"]["sharp_energy"])
    h2 = int(rows["hard_e2_quartic"]["sharp_energy"])
    h3 = int(rows["hard_e3_quartic_C1"]["sharp_energy"])
    opposite = int(rows["opposite_q4_quartic"]["sharp_energy"])

    partition_32_base = nonexact_parseval_base((3, 2))
    partition_32_upper = h3 + h2 + 7 * opposite
    partition_32_lower = partition_32_base + 26

    partition_41_base = nonexact_parseval_base((4, 1))
    partition_41_raw_upper = int(rows["hard_e4_quartic_raw"]["sharp_energy"]) + h1 + 7 * opposite
    partition_41_collision_upper = (partition_41_raw_upper - partition_41_base) // 26
    partition_41_cases = []
    for collision in (2, 3):
        upper = int(rows[f"hard_e4_quartic_C{collision}"]["sharp_energy"]) + h1 + 7 * opposite
        lower = partition_41_base + 26 * collision
        partition_41_cases.append(
            {"C": collision, "lower": lower, "upper": upper, "gap": lower - upper}
        )

    partition_5_base = nonexact_parseval_base((5,))
    partition_5_raw_upper = int(rows["hard_e5_quartic_raw"]["sharp_energy"]) + 7 * opposite
    partition_5_collision_upper = (partition_5_raw_upper - partition_5_base) // 26
    partition_5_cases = []
    for collision in range(3, 8):
        upper = int(rows[f"hard_e5_quartic_C{collision}"]["sharp_energy"]) + 7 * opposite
        lower = partition_5_base + 26 * collision
        partition_5_cases.append(
            {"C": collision, "lower": lower, "upper": upper, "gap": lower - upper}
        )

    proved = bool(
        partition_32_base == 284
        and partition_32_upper == 288
        and partition_32_lower == 310
        and partition_41_base == 276
        and partition_41_raw_upper == 370
        and partition_41_collision_upper == 3
        and [row["gap"] for row in partition_41_cases] == [28, 50]
        and partition_5_base == 259
        and partition_5_raw_upper == 443
        and partition_5_collision_upper == 7
        and [row["gap"] for row in partition_5_cases] == [38, 42, 56, 58, 82]
        and all(row["proved"] for row in rows.values())
    )
    _require(proved, "a high-root p13 u6 partition remained")
    return {
        "scope": "only p=13,t=4,u=6 excess partitions (3,2), (4,1), and (5)",
        "common_normalization": {
            "hard_edge_count": 33,
            "opposite_edge_count": 28,
            "hT": 5,
            "hard_parallel_count": "P_L=4+e_L",
            "opposite_parallel_counts": [4] * 7,
            "hard_excess_sum": 5,
        },
        "global_quartic": {
            "polynomial": "G=h*M4-M2^2",
            "at_least_five_exact_XNOR_roots": True,
            "therefore_G_identically_zero": True,
            "hard_relation": "N4=N2^2",
            "opposite_relation": "N4=-N2^2",
        },
        "collision_floor": collision_floor_certificate(),
        "row_certificates": rows,
        "partitions": [
            {
                "partition": [3, 2],
                "exact_XNOR_roots": 5,
                "collision_minimum": 1,
                "nonexact_Parseval_base": partition_32_base,
                "lower": partition_32_lower,
                "upper": partition_32_upper,
                "gap": partition_32_lower - partition_32_upper,
                "excluded": partition_32_lower > partition_32_upper,
            },
            {
                "partition": [4, 1],
                "exact_XNOR_roots": 5,
                "collision_minimum": 2,
                "raw_collision_upper": partition_41_collision_upper,
                "collision_cases": partition_41_cases,
                "excluded": all(row["gap"] > 0 for row in partition_41_cases),
            },
            {
                "partition": [5],
                "exact_XNOR_roots": 6,
                "collision_minimum": 3,
                "raw_collision_upper": partition_5_collision_upper,
                "collision_cases": partition_5_cases,
                "excluded": all(row["gap"] > 0 for row in partition_5_cases),
            },
        ],
        "closed_partitions": [[3, 2], [4, 1], [5]],
        "p13_t4_u6_fully_closed": False,
        "graph_or_configuration_census_used": False,
        "finite_six_bin_aggregate_models_used": True,
        "proved": proved,
    }


__all__ = [
    "ROW_SPECS",
    "ROW_SPEC_BY_NAME",
    "EXPECTED_MODEL_HASHES",
    "balanced_collision",
    "collision_floor_certificate",
    "high_root_partition_certificate",
    "moment_relation_holds",
    "nonexact_parseval_base",
    "row_energy_certificate",
    "row_satisfies",
    "signed_collision_floor",
    "translated_cuts",
]
