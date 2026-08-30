#!/usr/bin/env python3
"""Bounded joint mod-3/mod-7 supports for the positive p=7 high catalogs.

The supplied Normaliz ``.gen`` file is interpreted as a Hilbert basis for

    S = {L in Z^35_{≥0} : K L = 0}.

Every z=7 slack has the form ``parity_floor + 2 L``.  If the excess grade is
``g = sum(L)/5``, the projected support therefore obeys the exact recurrence

    T[0] = {0},
    T[g] = union_(h in H, deg(h) <= g) (T[g-deg(h)] - 2 D h).

Here ``D`` is a block of the derived torsion-dependency quotient of one of
the existing pointed systems.  A Hilbert-basis row is projected modulo 3 and
7 first and the two components are concatenated before deduplication.  Thus
one generator row, not independently chosen marginal rows, supplies both
prime components.

The computation is deliberately bounded.  A state-cap hit makes that grade
and every dependent later grade an explicit skip; partial supports are never
used for rejection.  Exact finite-group saturation is reported as saturation
and has no exclusion force.  Grades at most six automatically satisfy the
slack box because the audited basis is binary.  Grades seven and eight are a
safe outer support with the coordinate box dropped.  The script makes no z=7
closure claim and currently permits smoke runs only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z7_compact_symmetry_audit as symmetry  # noqa: E402
import p7_infinity7_positive_z7_mod7_projection as parent  # noqa: E402
import p7_infinity7_positive_z7_pointed_affine_hull_multimod as affine  # noqa: E402
import p7_infinity7_positive_z7_torsion_support_projection as torsion  # noqa: E402
import p7_size_four_slack_classify as johnson  # noqa: E402


P = 7
AMBIENT_DIMENSION = 35
GRADING_DENOMINATOR = 5
MAX_EXCESS_GRADE = 8
MODULI = (3, 7)
EXPECTED_GENERATOR_ROWS = 896
EXPECTED_DEGREE_HISTOGRAM = {1: 56, 2: 168, 3: 672}
CURRENT_REFERENCE_SHA256 = "3b582d6a0e7c83cb8ed41a421e4950be2645ce2d3aa18dab17432628b787b789"
DEFAULT_HILBERT_BASIS = Path("/tmp/p7_johnson_semigroup.gen")
DEFAULT_STATE_CAP = 200_000
DEFAULT_PAIR_CHUNK_CAP = 250_000
MAX_SMOKE_STATE_CAP = 2_000_000
MAX_SMOKE_MOD7_COORDINATES = 6


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def sorted_unique_integer_rows(rows: np.ndarray) -> np.ndarray:
    source = np.ascontiguousarray(rows, dtype=np.int16)
    require(source.ndim == 2 and source.shape[1] == AMBIENT_DIMENSION, "integer row shape changed")
    return np.ascontiguousarray(np.unique(source, axis=0), dtype=np.int16)


def load_hilbert_basis(path: Path, kernel_rows: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    """Parse and structurally audit a Normaliz ``.gen`` Hilbert-basis file."""
    require(path.is_file(), f"Normaliz Hilbert-basis file does not exist: {path}")
    raw = path.read_bytes()
    lines = [line.strip() for line in raw.decode("ascii").splitlines() if line.strip()]
    require(len(lines) >= 2, "Normaliz .gen file is truncated")
    try:
        row_count = int(lines[0])
        dimension = int(lines[1])
        parsed = [tuple(int(value) for value in line.split()) for line in lines[2:]]
    except ValueError as error:
        raise AssertionError("Normaliz .gen file contains a non-integer field") from error
    require(row_count == EXPECTED_GENERATOR_ROWS, "Normaliz generator-row count changed")
    require(dimension == AMBIENT_DIMENSION, "Normaliz ambient dimension changed")
    require(len(parsed) == row_count, "Normaliz .gen row count disagrees with its header")
    require(all(len(row) == dimension for row in parsed), "Normaliz generator width changed")

    basis = np.ascontiguousarray(parsed, dtype=np.int64)
    require(np.all(basis >= 0), "Normaliz basis contains a negative entry")
    require(len(np.unique(basis, axis=0)) == row_count, "Normaliz basis contains duplicate rows")
    require(np.all((basis == 0) | (basis == 1)), "reference Hilbert basis is no longer binary")
    require(not np.any(np.asarray(kernel_rows, dtype=np.int64) @ basis.T), "basis escaped the exact Johnson kernel")

    masses = basis.sum(axis=1, dtype=np.int64)
    require(np.all(masses > 0), "Normaliz basis contains the zero row")
    require(not np.any(masses % GRADING_DENOMINATOR), "basis mass is not divisible by the grading denominator")
    degrees = np.ascontiguousarray(masses // GRADING_DENOMINATOR, dtype=np.int64)
    histogram = dict(sorted(Counter(int(value) for value in degrees).items()))
    require(histogram == EXPECTED_DEGREE_HISTOGRAM, "Hilbert generator-degree histogram changed")

    # In low grades the claimed Hilbert basis can be checked against complete
    # repository catalogs without trusting a projected calculation.
    grade_one = sorted_unique_integer_rows(basis[degrees == 1])
    direct_grade_one = sorted_unique_integer_rows(
        affine.canonical_catalog(7, 4).astype(np.int64) // 2
    )
    require(np.array_equal(grade_one, direct_grade_one), "degree-one basis disagrees with S56")

    decomposable_grade_two = (
        grade_one[:, None, :].astype(np.int16)
        + grade_one[None, :, :].astype(np.int16)
    ).reshape(-1, AMBIENT_DIMENSION)
    grade_two_slice = sorted_unique_integer_rows(
        np.vstack((basis[degrees == 2], decomposable_grade_two))
    )
    direct_grade_two = sorted_unique_integer_rows(
        affine.canonical_catalog(7, 8).astype(np.int64) // 2
    )
    require(len(grade_two_slice) == 1_764, "degree-two semigroup slice census changed")
    require(np.array_equal(grade_two_slice, direct_grade_two), "degree-two recurrence disagrees with M1764")

    digest = hashlib.sha256(raw).hexdigest()
    return basis, degrees, {
        "path": str(path.resolve()),
        "file_bytes": len(raw),
        "file_sha256": digest,
        "matches_current_reference_sha256": digest == CURRENT_REFERENCE_SHA256,
        "header_generator_rows": row_count,
        "header_ambient_dimension": dimension,
        "grading_denominator": GRADING_DENOMINATOR,
        "generator_degree_histogram": {str(key): value for key, value in histogram.items()},
        "all_rows_nonnegative_distinct_binary_and_nonzero": True,
        "all_rows_in_exact_Johnson_kernel": True,
        "degree_one_rows_equal_complete_S56_lift_catalog": True,
        "degree_two_semigroup_slice_rows": len(grade_two_slice),
        "degree_two_semigroup_slice_equals_complete_M1764_lift_catalog": True,
        "basis_sha256_int64": array_sha256(basis),
        "degrees_sha256_int64": array_sha256(degrees),
        "trust_boundary": (
            "The .gen file is consumed as a complete Normaliz Hilbert basis. "
            "Its structure and complete grade-one/two slices are independently audited here; "
            "Normaliz's all-grade Hilbert-basis completeness is not independently reproved."
        ),
        "all_grade_completeness_independently_reproved": False,
    }


def derive_torsion_quotients(
    matrix: np.ndarray, common: np.ndarray
) -> tuple[dict[int, dict], dict]:
    """Derive the varying modular quotient for either pointed branch.

    The earlier torsion prototype pins branch A's rational rank to 168.  This
    version proves the rank from the selected system itself, so branch B's
    seven additional independent fixed-edge rows (rational rank 175) are
    handled without changing the quotient argument.
    """
    modular_data = {}
    rank_rows = []
    for modulus in (3, 7, 11):
        matrix_rank, dependencies = torsion.left_dependencies(matrix, modulus)
        dependencies = np.ascontiguousarray(dependencies, dtype=np.int64)
        require(
            not np.any(dependencies @ (matrix % modulus) % modulus),
            f"complete mod-{modulus} left kernel failed direct audit",
        )
        require(
            affine.modular_rank(dependencies, modulus) == len(dependencies),
            f"complete mod-{modulus} dependency basis is not independent",
        )
        modular_data[modulus] = {
            "matrix_rank": int(matrix_rank),
            "dependencies": dependencies,
        }
        rank_rows.append(
            {
                "modulus": modulus,
                "matrix_rank": int(matrix_rank),
                "left_dependency_dimension": len(dependencies),
                "dependency_sha256_uint8": array_sha256(dependencies.astype(np.uint8)),
                "left_null_and_full_basis_audit": True,
            }
        )

    common_rank_mod11 = affine.modular_rank(common, 11)
    require(common_rank_mod11 == len(common) == 114, "common exact dependency rank changed")
    rational_rank_upper = matrix.shape[0] - common_rank_mod11
    rational_rank_lower = modular_data[11]["matrix_rank"]
    require(
        rational_rank_upper == rational_rank_lower,
        "selected pointed system's exact rational rank was not pinned",
    )

    quotient_data: dict[int, dict] = {}
    quotient_rows = []
    for modulus in MODULI:
        complete = modular_data[modulus]["dependencies"]
        require(
            affine.modular_rank(np.vstack((complete, common)), modulus) == len(complete),
            f"common dependencies escaped the complete mod-{modulus} kernel",
        )
        complement, indices, audit = torsion.extend_common_to_full(
            common, complete, modulus
        )
        quotient_data[modulus] = {
            "complete": complete,
            "complement": complement,
            "selected_indices": indices,
        }
        quotient_rows.append(audit)

    dimensions = {
        modulus: len(quotient_data[modulus]["complement"]) for modulus in MODULI
    }
    return quotient_data, {
        "integer_pointed_matrix_shape": list(matrix.shape),
        "integer_pointed_matrix_sha256_int64": array_sha256(matrix),
        "exact_common_dependency_dimension": len(common),
        "rational_rank_upper_from_complete_exact_dependencies": rational_rank_upper,
        "rational_rank_lower_from_nonzero_mod11_minor": rational_rank_lower,
        "exact_rational_rank": rational_rank_upper,
        "exact_rational_left_dependency_dimension": matrix.shape[0] - rational_rank_upper,
        "complete_modular_rank_audits": rank_rows,
        "quotient_audits": quotient_rows,
        "actual_effective_dimensions": {
            f"F{modulus}": dimension for modulus, dimension in dimensions.items()
        },
        "dimensions_F3_6_F7_21_reproduced": dimensions == {3: 6, 7: 21},
        "rational_rank_derived_from_selected_system_not_branch_A_hardcoded": True,
    }


class MixedRadixCodec:
    """Collision-free encoding for a bounded product of prime fields."""

    def __init__(self, moduli: tuple[int, ...]):
        require(moduli and all(value > 1 for value in moduli), "mixed-radix moduli are invalid")
        self.moduli = tuple(int(value) for value in moduli)
        self.group_size = math.prod(self.moduli)
        require(self.group_size <= np.iinfo(np.uint64).max, "projected group does not fit uint64")
        multipliers = [1]
        for modulus in self.moduli[:-1]:
            multipliers.append(multipliers[-1] * modulus)
        self.multipliers = np.asarray(multipliers, dtype=np.uint64)
        self.moduli_array = np.asarray(self.moduli, dtype=np.uint16)

    def encode(self, rows: np.ndarray) -> np.ndarray:
        source = np.ascontiguousarray(rows, dtype=np.uint8)
        require(
            source.ndim == 2 and source.shape[1] == len(self.moduli),
            "mixed-radix row width changed",
        )
        for column, modulus in enumerate(self.moduli):
            require(np.all(source[:, column] < modulus), "mixed-radix digit escaped its field")
        return np.ascontiguousarray(
            np.sum(source.astype(np.uint64) * self.multipliers[None, :], axis=1),
            dtype=np.uint64,
        )

    def decode(self, codes: np.ndarray) -> np.ndarray:
        source = np.ascontiguousarray(codes, dtype=np.uint64)
        require(source.ndim == 1, "mixed-radix codes are not a vector")
        require(not len(source) or int(source.max()) < self.group_size, "mixed-radix code escaped group")
        return np.ascontiguousarray(
            np.column_stack(
                [
                    (source // self.multipliers[index] % modulus).astype(np.uint8)
                    for index, modulus in enumerate(self.moduli)
                ]
            ),
            dtype=np.uint8,
        )

    def unique_codes(self, rows: np.ndarray) -> np.ndarray:
        return np.ascontiguousarray(np.unique(self.encode(rows)), dtype=np.uint64)

    def translate(self, codes: np.ndarray, offset: np.ndarray) -> np.ndarray:
        offset = np.ascontiguousarray(offset, dtype=np.uint8)
        require(offset.shape == (len(self.moduli),), "translation offset width changed")
        rows = self.decode(codes).astype(np.uint16)
        translated = (rows + offset[None, :].astype(np.uint16)) % self.moduli_array[None, :]
        return self.unique_codes(translated.astype(np.uint8))


def bounded_minkowski_union(
    left_codes: np.ndarray,
    right_codes: np.ndarray,
    accumulated: np.ndarray,
    codec: MixedRadixCodec,
    state_cap: int,
    pair_chunk_cap: int,
) -> tuple[np.ndarray, dict]:
    """Add one exact support product, stopping only at saturation or a hard cap."""
    left = codec.decode(left_codes)
    right = codec.decode(right_codes)
    require(len(left) and len(right), "Minkowski factor support is empty")
    union = np.ascontiguousarray(accumulated, dtype=np.uint64)
    blocks = 0
    pairs = 0
    for right_start in range(0, len(right), pair_chunk_cap):
        right_block = right[right_start : right_start + pair_chunk_cap]
        left_chunk = max(1, pair_chunk_cap // len(right_block))
        for left_start in range(0, len(left), left_chunk):
            left_block = left[left_start : left_start + left_chunk]
            sums = (
                left_block[:, None, :].astype(np.uint16)
                + right_block[None, :, :].astype(np.uint16)
            ) % codec.moduli_array[None, None, :]
            block_codes = np.unique(
                codec.encode(sums.reshape(-1, len(codec.moduli)).astype(np.uint8))
            )
            union = np.ascontiguousarray(np.union1d(union, block_codes), dtype=np.uint64)
            blocks += 1
            pairs += len(left_block) * len(right_block)
            if len(union) == codec.group_size:
                return union, {
                    "status": "exact_full_group_saturation",
                    "pair_blocks": blocks,
                    "candidate_pairs_examined": pairs,
                }
            if len(union) > state_cap:
                return union, {
                    "status": "state_cap_exceeded",
                    "pair_blocks": blocks,
                    "candidate_pairs_examined": pairs,
                }
    return union, {
        "status": "complete",
        "pair_blocks": blocks,
        "candidate_pairs_examined": pairs,
    }


def support_recurrence(
    generator_codes: dict[int, np.ndarray],
    codec: MixedRadixCodec,
    state_cap: int,
    pair_chunk_cap: int,
    max_grade: int = MAX_EXCESS_GRADE,
) -> tuple[dict[int, np.ndarray], list[dict]]:
    """Compute exact unbounded-semigroup image supports until an explicit cap."""
    require(set(generator_codes) == {1, 2, 3}, "generator degree support changed")
    require(all(len(generator_codes[degree]) for degree in generator_codes), "empty generator degree")
    supports: dict[int, np.ndarray] = {0: np.asarray([0], dtype=np.uint64)}
    records = [
        {
            "excess_grade": 0,
            "decision_status": "complete_exact_unbounded_semigroup_projection",
            "completed": True,
            "skipped": False,
            "state_count": 1,
            "support_sha256_uint64": array_sha256(supports[0].astype("<u8", copy=False)),
            "full_projected_group_saturated": codec.group_size == 1,
            "coordinate_box_status": "exact_for_catalog",
        }
    ]
    first_incomplete_grade: int | None = None
    for grade in range(1, max_grade + 1):
        if first_incomplete_grade is not None:
            records.append(
                {
                    "excess_grade": grade,
                    "decision_status": "skipped_incomplete_predecessor_grade",
                    "completed": False,
                    "skipped": True,
                    "skip_reason": f"grade {first_incomplete_grade} exceeded the state cap",
                    "state_count": None,
                    "support_sha256_uint64": None,
                    "partial_support_used": False,
                    "rigorous_rejection_allowed": False,
                }
            )
            continue

        accumulated = np.empty(0, dtype=np.uint64)
        degree_terms = []
        saturated = False
        capped = False
        for degree in (1, 2, 3):
            if degree > grade:
                continue
            predecessor = supports[grade - degree]
            generators = generator_codes[degree]
            if len(predecessor) == codec.group_size:
                accumulated = np.arange(codec.group_size, dtype=np.uint64)
                term = {
                    "generator_degree": degree,
                    "predecessor_grade": grade - degree,
                    "predecessor_states": len(predecessor),
                    "generator_states": len(generators),
                    "status": "exact_full_group_translation_shortcut",
                    "pair_blocks": 0,
                    "candidate_pairs_examined": 0,
                }
                degree_terms.append(term)
                saturated = True
                break
            accumulated, term_audit = bounded_minkowski_union(
                predecessor,
                generators,
                accumulated,
                codec,
                state_cap,
                pair_chunk_cap,
            )
            degree_terms.append(
                {
                    "generator_degree": degree,
                    "predecessor_grade": grade - degree,
                    "predecessor_states": len(predecessor),
                    "generator_states": len(generators),
                    **term_audit,
                }
            )
            if term_audit["status"] == "exact_full_group_saturation":
                saturated = True
                break
            if term_audit["status"] == "state_cap_exceeded":
                capped = True
                break

        if capped:
            first_incomplete_grade = grade
            records.append(
                {
                    "excess_grade": grade,
                    "decision_status": "skipped_state_cap",
                    "completed": False,
                    "skipped": True,
                    "skip_reason": "exact accumulated support exceeded --state-cap",
                    "state_cap": state_cap,
                    "distinct_states_lower_bound_at_skip": len(accumulated),
                    "support_sha256_uint64": None,
                    "degree_terms": degree_terms,
                    "partial_support_used": False,
                    "skip_is_explicit_not_approximation": True,
                    "rigorous_rejection_allowed": False,
                }
            )
            continue

        require(len(accumulated) <= state_cap, "completed support escaped state cap")
        require(len(accumulated) > 0, "completed recurrence support is empty")
        supports[grade] = accumulated
        equal_previous = np.array_equal(accumulated, supports[grade - 1])
        box_exact = grade <= 6
        records.append(
            {
                "excess_grade": grade,
                "decision_status": (
                    "complete_exact_full_projected_group_saturation"
                    if saturated
                    else "complete_exact_unbounded_semigroup_projection"
                ),
                "completed": True,
                "skipped": False,
                "state_count": len(accumulated),
                "state_cap": state_cap,
                "finite_projected_group_size": codec.group_size,
                "support_sha256_uint64": array_sha256(accumulated.astype("<u8", copy=False)),
                "full_projected_group_saturated": saturated,
                "equal_to_previous_computed_grade_support": equal_previous,
                "equality_is_not_extrapolated_to_uncomputed_grades": True,
                "coordinate_box_status": (
                    "exact_for_catalog"
                    if box_exact
                    else "safe_outer_support_coordinate_upper_bound_dropped"
                ),
                "exact_bounded_catalog_projection_conditional_on_supplied_Hilbert_basis": box_exact,
                "rigorous_catalog_outer_support": True,
                "saturation_or_presence_does_not_prove_catalog_or_edge_feasibility": True,
                "degree_terms": degree_terms,
            }
        )
    require(len(records) == max_grade + 1, "recurrence grade coverage changed")
    return supports, records


def project_rows(
    rows: np.ndarray,
    direction: int,
    quotient_data: dict,
    mod7_coordinates: tuple[int, ...],
) -> np.ndarray:
    """Project the same integer row into both prime components."""
    source = np.ascontiguousarray(rows, dtype=np.int64)
    require(source.ndim == 2 and source.shape[1] == AMBIENT_DIMENSION, "projected source shape changed")
    block = slice(1 + AMBIENT_DIMENSION * direction, 1 + AMBIENT_DIMENSION * (direction + 1))
    components = []
    for modulus in MODULI:
        complement = np.asarray(quotient_data[modulus]["complement"], dtype=np.int64)
        selected = tuple(range(len(complement))) if modulus == 3 else mod7_coordinates
        block_rows = complement[list(selected), block]
        component = block_rows @ (source.T % modulus) % modulus
        components.append(np.ascontiguousarray(component.T, dtype=np.uint8))
    joint = np.ascontiguousarray(np.concatenate(components, axis=1), dtype=np.uint8)
    require(len(joint) == len(source), "joint projection lost source-row identity")
    return joint


def excess_grade(mask: int, scaled_mean: int) -> tuple[int, int, np.ndarray]:
    floor = affine.parity_for_mask(mask).astype(np.int64)
    numerator = 2 * int(floor.sum())
    require(numerator % 5 == 0, "parity floor has nonintegral scaled mean")
    floor_mean = numerator // 5
    require(scaled_mean >= floor_mean and (scaled_mean - floor_mean) % 4 == 0, "bad excess grade")
    return (scaled_mean - floor_mean) // 4, floor_mean, floor


def high_leaf_audit(context: dict) -> tuple[list[list[dict]], dict]:
    leaves_by_orbit, leaf_source = parent.exact_mean_leaves(context["orbits"])
    grade_histogram: Counter[int] = Counter()
    mask_mean_grade: Counter[tuple[int, int, int]] = Counter()
    per_orbit_direction_grades = []
    for orbit_index, (orbit, leaves) in enumerate(zip(context["orbits"], leaves_by_orbit)):
        by_direction = {direction: set() for direction in range(P + 1)}
        for leaf in leaves:
            for direction in leaf["high_directions"]:
                mask = int(orbit["masks"][direction])
                mean = int(leaf["scaled_means"][direction])
                grade, _floor_mean, _floor = excess_grade(mask, mean)
                level = int(leaf["catalog_levels"][direction])
                expected_grade = 2 * level - int(mean % 8 == 4)
                require(grade == expected_grade, "leaf level/excess grade conversion changed")
                grade_histogram[grade] += 1
                mask_mean_grade[(mask.bit_count(), mean, grade)] += 1
                by_direction[int(direction)].add(grade)
        per_orbit_direction_grades.append(
            {
                "branch_orbit_index": orbit_index,
                "grades_by_direction": {
                    str(direction): sorted(values) for direction, values in by_direction.items()
                },
            }
        )
    require(grade_histogram and min(grade_histogram) >= 3, "high catalog grade floor changed")
    require(max(grade_histogram) == MAX_EXCESS_GRADE, "high catalog maximum grade changed")
    require(set(grade_histogram) <= set(range(MAX_EXCESS_GRADE + 1)), "high grade escaped recurrence bound")
    return leaves_by_orbit, {
        "mean_leaf_source": leaf_source,
        "high_direction_occurrences": sum(grade_histogram.values()),
        "high_excess_grade_histogram": {
            str(key): value for key, value in sorted(grade_histogram.items())
        },
        "high_mask_mean_grade_histogram": {
            f"b{b}_mean{mean}_grade{grade}": value
            for (b, mean, grade), value in sorted(mask_mean_grade.items())
        },
        "per_orbit_direction_required_grades": per_orbit_direction_grades,
        "all_actual_high_catalogs_have_excess_grade_at_most_eight": True,
    }


def direct_catalog_calibration(
    direction: int,
    mask: int,
    anchors: affine.AnchorFactory,
    quotient_data: dict,
    mod7_coordinates: tuple[int, ...],
    codec: MixedRadixCodec,
    supports: dict[int, np.ndarray],
) -> list[dict]:
    """Compare recurrence grades one/two to the complete S/M catalogs."""
    rows = []
    for grade, catalog_class, expected_rows in ((1, "S", 56), (2, "M", 1_764)):
        if grade not in supports:
            rows.append(
                {
                    "excess_grade": grade,
                    "catalog_class": catalog_class,
                    "calibration_status": "skipped_recurrence_grade_incomplete",
                    "recurrence_support_equals_direct_complete_catalog_support": None,
                    "skip_is_not_a_failed_calibration": True,
                }
            )
            continue
        _unused, floor_mean, floor = excess_grade(mask, 0 if mask.bit_count() == 7 else 8)
        mean = floor_mean + 4 * grade
        anchor = anchors.get(mask, mean)
        offset = project_rows(
            (anchor - floor)[None, :], direction, quotient_data, mod7_coordinates
        )[0]
        recurrence = codec.translate(supports[grade], offset)
        catalog = affine.mapped_catalog(mask, mean).astype(np.int64)
        require(len(catalog) == expected_rows, f"complete {catalog_class} catalog size changed")
        delta = np.ascontiguousarray(anchor[None, :] - catalog, dtype=np.int64)
        require(not np.any(delta.sum(axis=1)), "direct catalog delta changed mean")
        require(not np.any(anchors.kernel_rows @ delta.T), "direct catalog delta left degree two")
        require(not np.any(delta % 2), "direct catalog delta changed parity")
        direct = codec.unique_codes(
            project_rows(delta, direction, quotient_data, mod7_coordinates)
        )
        require(np.array_equal(recurrence, direct), f"grade-{grade} projected recurrence calibration failed")
        rows.append(
            {
                "excess_grade": grade,
                "catalog_class": catalog_class,
                "scaled_mean": mean,
                "complete_catalog_rows": len(catalog),
                "projected_unique_states": len(direct),
                "direct_support_sha256_uint64": array_sha256(direct.astype("<u8", copy=False)),
                "recurrence_support_equals_direct_complete_catalog_support": True,
                "same_direct_catalog_row_used_mod3_mod7": True,
            }
        )
    return rows


def direction_support(
    direction: int,
    mask: int,
    basis: np.ndarray,
    degrees: np.ndarray,
    anchors: affine.AnchorFactory,
    common: np.ndarray,
    quotient_data: dict,
    mod7_coordinates: tuple[int, ...],
    codec: MixedRadixCodec,
    state_cap: int,
    pair_chunk_cap: int,
) -> dict:
    mapped = np.ascontiguousarray(
        np.stack([affine.map_canonical_row(mask, row) for row in basis]),
        dtype=np.int64,
    )
    require(not np.any(anchors.kernel_rows @ mapped.T), "mapped Hilbert row left degree two")
    variation = np.ascontiguousarray(-2 * mapped, dtype=np.int64)
    joint = project_rows(variation, direction, quotient_data, mod7_coordinates)
    q3 = len(quotient_data[3]["complement"])
    component_hashes = {
        "3": array_sha256(joint[:, :q3]),
        "7": array_sha256(joint[:, q3:]),
    }
    generator_codes = {
        degree: codec.unique_codes(joint[degrees == degree]) for degree in (1, 2, 3)
    }
    supports, records = support_recurrence(
        generator_codes, codec, state_cap, pair_chunk_cap
    )

    # Translate each raw excess support to the anchor-relative convention used
    # by the existing z=7 global catalog join: delta = anchor - catalog.
    translated_records = []
    for record in records:
        grade = int(record["excess_grade"])
        if not record["completed"]:
            translated_records.append(record)
            continue
        _dummy_grade, floor_mean, floor = excess_grade(
            mask, 0 if mask.bit_count() == 7 else 8
        )
        mean = floor_mean + 4 * grade
        anchor = anchors.get(mask, mean)
        anchor_minus_floor = np.ascontiguousarray(anchor - floor, dtype=np.int64)
        require(not np.any(anchors.kernel_rows @ anchor_minus_floor), "anchor offset left degree two")
        require(int(anchor_minus_floor.sum()) == 10 * grade, "anchor offset mass changed")
        offset = project_rows(
            anchor_minus_floor[None, :], direction, quotient_data, mod7_coordinates
        )[0]
        translated = codec.translate(supports[grade], offset)
        translated_records.append(
            {
                **record,
                "scaled_mean_for_this_mask": mean,
                "parity_floor_scaled_mean": floor_mean,
                "anchor_relative_support_state_count": len(translated),
                "anchor_relative_support_sha256_uint64": array_sha256(
                    translated.astype("<u8", copy=False)
                ),
                "anchor_relative_delta_formula": "anchor - parity_floor - 2*L",
                "every_recurrence_path_has_exact_total_excess_grade": True,
                "common_rational_dependencies_are_fixed_by_kernel_and_grade": True,
            }
        )

    calibration = direct_catalog_calibration(
        direction,
        mask,
        anchors,
        quotient_data,
        mod7_coordinates,
        codec,
        supports,
    )
    first_incomplete = next(
        (int(row["excess_grade"]) for row in translated_records if not row["completed"]),
        None,
    )
    required_block = slice(
        1 + AMBIENT_DIMENSION * direction,
        1 + AMBIENT_DIMENSION * (direction + 1),
    )
    require(
        all(
            not np.any(common[:, required_block] @ (mapped[degrees == degree].T - mapped[degrees == degree][0][:, None]))
            for degree in (1, 2, 3)
        ),
        "equal-grade generators are visible to a common rational dependency",
    )
    return {
        "direction": direction,
        "mask": mask,
        "mask_weight_b": mask.bit_count(),
        "mapped_basis_sha256_int64": array_sha256(mapped),
        "joint_generator_signature_shape": list(joint.shape),
        "joint_generator_signature_sha256_uint8": array_sha256(joint),
        "component_signature_sha256_uint8": component_hashes,
        "same_Hilbert_generator_row_index_used_mod3_mod7_before_deduplication": True,
        "projected_unique_generator_states_by_degree": {
            str(degree): len(generator_codes[degree]) for degree in (1, 2, 3)
        },
        "projected_unique_generator_support_sha256_uint64_by_degree": {
            str(degree): array_sha256(generator_codes[degree].astype("<u8", copy=False))
            for degree in (1, 2, 3)
        },
        "known_complete_catalog_calibration": calibration,
        "recurrence_grades": translated_records,
        "first_incomplete_grade": first_incomplete,
        "all_grades_through_eight_completed": first_incomplete is None,
        "no_partial_support_used_after_cap": True,
    }


def manufactured_self_audit() -> dict:
    """Brute-force recurrence, cross-prime identity, cap, and saturation traps."""
    codec = MixedRadixCodec((3, 7))
    generator_rows = {
        1: np.asarray([[0, 1], [1, 0]], dtype=np.uint8),
        2: np.asarray([[2, 3]], dtype=np.uint8),
        3: np.asarray([[1, 6]], dtype=np.uint8),
    }
    generator_codes = {
        degree: codec.unique_codes(rows) for degree, rows in generator_rows.items()
    }
    supports, records = support_recurrence(generator_codes, codec, 100, 100, max_grade=5)
    brute: dict[int, set[tuple[int, int]]] = {0: {(0, 0)}}
    for grade in range(1, 6):
        brute[grade] = {
            ((left[0] + right[0]) % 3, (left[1] + right[1]) % 7)
            for degree, rows in generator_rows.items()
            if degree <= grade
            for left in brute[grade - degree]
            for right in map(tuple, rows.tolist())
        }
        observed = set(map(tuple, codec.decode(supports[grade]).tolist()))
        require(observed == brute[grade], "manufactured recurrence disagrees with brute force")

    trap = set(map(tuple, generator_rows[1].tolist()))
    require(
        (0, 0) not in trap
        and any(row[0] == 0 for row in trap)
        and any(row[1] == 0 for row in trap),
        "same-row cross-prime trap was not exercised",
    )

    cap_codec = MixedRadixCodec((5, 5))
    all_rows = np.asarray(list(np.ndindex(5, 5)), dtype=np.uint8)
    capped_codes = {
        1: cap_codec.unique_codes(all_rows),
        2: cap_codec.unique_codes(np.asarray([[0, 0]], dtype=np.uint8)),
        3: cap_codec.unique_codes(np.asarray([[0, 0]], dtype=np.uint8)),
    }
    _capped_supports, capped_records = support_recurrence(
        capped_codes, cap_codec, 10, 10, max_grade=3
    )
    require(capped_records[1]["decision_status"] == "skipped_state_cap", "state-cap trap failed")
    require(
        all(row["decision_status"] == "skipped_incomplete_predecessor_grade" for row in capped_records[2:]),
        "dependent grades did not skip after cap",
    )

    saturated_supports, saturated_records = support_recurrence(
        capped_codes, cap_codec, 25, 25, max_grade=2
    )
    require(
        saturated_records[1]["full_projected_group_saturated"] is True
        and len(saturated_supports[1]) == 25,
        "finite-group saturation trap failed",
    )
    return {
        "brute_force_recurrence_grades_checked": 5,
        "brute_force_support_sizes": {str(grade): len(rows) for grade, rows in brute.items()},
        "same_row_cross_prime_false_positive_trap_rejected": True,
        "cap_produces_explicit_skip_and_blocks_dependent_grades": True,
        "full_group_saturation_detected_exactly": True,
        "saturation_is_not_reported_as_feasibility": True,
        "audit_sha256": json_sha256(
            {
                "records": records,
                "capped_records": capped_records,
                "saturated_records": saturated_records,
            }
        ),
    }


def parse_integer_tuple(value: str, label: str) -> tuple[int, ...]:
    try:
        rows = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as error:
        raise AssertionError(f"{label} must be a comma-separated integer list") from error
    require(len(rows) == len(set(rows)), f"{label} repeats an index")
    return rows


def run(
    *,
    smoke: bool,
    hilbert_basis_path: Path,
    orbit_index: int,
    branch: str,
    directions: tuple[int, ...],
    mod7_coordinates: tuple[int, ...],
    state_cap: int,
    pair_chunk_cap: int,
) -> dict:
    started = time.time()
    require(smoke, "this first implementation deliberately supports --smoke only")
    require(orbit_index in (0, 1), "orbit index must be zero or one")
    require(branch in ("A", "B"), "pointed branch must be A or B")
    require(directions and len(directions) == len(set(directions)), "directions must be nonempty and unique")
    require(all(0 <= value <= P for value in directions), "direction index is outside 0..7")
    require(0 < state_cap <= MAX_SMOKE_STATE_CAP, "smoke state cap is outside its audited bound")
    require(0 < pair_chunk_cap <= MAX_SMOKE_STATE_CAP, "pair chunk cap is outside its audited bound")
    require(
        len(mod7_coordinates) <= MAX_SMOKE_MOD7_COORDINATES,
        "too many mod-7 coordinates for bounded smoke mode",
    )

    context = symmetry.construct_pointed_systems()
    system = context["systems"][orbit_index][branch]
    matrix = np.ascontiguousarray(system["matrix"], dtype=np.int64)
    require(matrix.shape[1] == 1_225, "pointed edge-variable count changed")
    kernel_rows = np.asarray(johnson._primitive_left_kernel_rows(), dtype=np.int64)  # noqa: SLF001
    require(kernel_rows.shape == (14, AMBIENT_DIMENSION), "primitive Johnson kernel shape changed")
    basis, degrees, basis_audit = load_hilbert_basis(hilbert_basis_path, kernel_rows)

    common, common_audit = torsion.exact_common_dependency_basis(matrix, kernel_rows)
    quotient_data, quotient_audit = derive_torsion_quotients(matrix, common)
    q3 = len(quotient_data[3]["complement"])
    q7 = len(quotient_data[7]["complement"])
    require(q3 == 6 and q7 == 21, "derived torsion quotient dimensions changed")
    require(
        all(0 <= value < q7 for value in mod7_coordinates),
        "selected mod-7 quotient coordinate is out of range",
    )
    moduli = (3,) * q3 + (7,) * len(mod7_coordinates)
    codec = MixedRadixCodec(moduli)

    leaves_by_orbit, leaf_audit = high_leaf_audit(context)
    kernel_from_hull, _hull_bases, hull_audit = affine.build_hull_audit()
    require(np.array_equal(kernel_from_hull, kernel_rows), "hull and semigroup kernels disagree")
    source_catalog = affine.canonical_catalog(7, 4).astype(np.int64)
    anchors = affine.AnchorFactory(kernel_rows, source_catalog[1:] - source_catalog[0])
    anchor_audit = anchors.audit_universe(context["orbits"], leaves_by_orbit)

    orbit = context["orbits"][orbit_index]
    results = [
        direction_support(
            direction,
            int(orbit["masks"][direction]),
            basis,
            degrees,
            anchors,
            common,
            quotient_data,
            mod7_coordinates,
            codec,
            state_cap,
            pair_chunk_cap,
        )
        for direction in directions
    ]
    all_complete = all(row["all_grades_through_eight_completed"] for row in results)
    return {
        "experiment": "p7_infinity7_positive_z7_high_semigroup_support",
        "status": (
            "complete_self_audited_bounded_smoke_high_semigroup_support"
            if all_complete
            else "bounded_smoke_high_semigroup_support_with_explicit_cap_skips"
        ),
        "p": P,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "smoke_test": True,
        "full_run": False,
        "configuration": {
            "selected_branch_orbit_index": orbit_index,
            "selected_pointed_star_branch": branch,
            "selected_directions": list(directions),
            "retained_all_mod3_torsion_coordinates": list(range(q3)),
            "retained_mod7_torsion_coordinates": list(mod7_coordinates),
            "projected_group": f"F3^{q3} x F7^{len(mod7_coordinates)}",
            "projected_group_size": codec.group_size,
            "state_cap": state_cap,
            "pair_chunk_cap": pair_chunk_cap,
            "maximum_excess_grade": MAX_EXCESS_GRADE,
        },
        "normaliz_Hilbert_basis": basis_audit,
        "pointed_system_source": {
            "all_four_pointed_systems_reconstructed_and_audited": True,
            "selected_system_matrix_shape": list(matrix.shape),
            "selected_system_matrix_sha256_int64": array_sha256(matrix),
            "selected_base_rhs_sha256_int64": array_sha256(
                np.asarray(system["base_rhs"], dtype=np.int64)
            ),
            "selected_fixed_edge_rows": system["fixed_edge_rows"],
            "orbit_source": context["orbit_source"],
        },
        "high_catalog_leaf_audit": leaf_audit,
        "exact_common_rational_dependency_audit": common_audit,
        "varying_torsion_quotient_audit": quotient_audit,
        "degree_two_zero_mean_hull_audit": hull_audit,
        "anchor_universe_summary": {
            "unique_mask_mean_anchors": anchor_audit["unique_mask_mean_anchors"],
            "direction_uses_across_all_2160_leaves": anchor_audit[
                "direction_uses_across_all_2160_leaves"
            ],
            "all_anchors_exact_degree_two_mean_and_parity": anchor_audit[
                "all_anchors_exact_degree_two_mean_and_parity"
            ],
        },
        "manufactured_self_audit": manufactured_self_audit(),
        "direction_results": results,
        "selected_directions_completed_through_grade_eight": sum(
            row["all_grades_through_eight_completed"] for row in results
        ),
        "selected_directions_with_explicit_cap_skip": sum(
            not row["all_grades_through_eight_completed"] for row in results
        ),
        "logical_semantics": {
            "generator_degree_recurrence_is_exact_for_the_supplied_Hilbert_basis": True,
            "same_Hilbert_generator_row_supplies_mod3_and_mod7_components": True,
            "deduplication_is_on_collision_free_joint_mixed_radix_codes": True,
            "grades_zero_through_six_automatically_obey_coordinate_upper_bound": True,
            "reason_grade_at_most_six_obeys_box": (
                "Every Hilbert generator is binary and every summand has positive degree, "
                "so a grade-g recurrence path uses at most g generators and each L coordinate is at most g<=6."
            ),
            "grades_seven_and_eight_are_safe_outer_supports_with_box_dropped": True,
            "absence_from_a_completed_outer_support_can_be_used_as_a_rigorous_necessary_rejection": True,
            "presence_or_saturation_is_only_necessary_and_proves_no_catalog_or_edge_lift": True,
            "state_cap_excess_is_an_explicit_skip_not_a_rejection": True,
            "partial_support_after_cap_is_discarded": True,
            "no_z7_exclusion_claimed": True,
            "no_binary_edge_feasibility_claimed": True,
        },
        "all_selected_supports_completed_without_cap": all_complete,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--hilbert-basis", type=Path, default=DEFAULT_HILBERT_BASIS)
    parser.add_argument("--orbit", type=int, default=0)
    parser.add_argument("--branch", choices=("A", "B"), default="A")
    parser.add_argument(
        "--directions",
        default=",".join(str(value) for value in range(P + 1)),
        help="comma-separated direction indices (default: all 0..7)",
    )
    parser.add_argument(
        "--mod7-coordinates",
        default="0,1",
        help="comma-separated derived F7 quotient coordinates; empty retains none",
    )
    parser.add_argument("--state-cap", type=int, default=DEFAULT_STATE_CAP)
    parser.add_argument("--pair-chunk-cap", type=int, default=DEFAULT_PAIR_CHUNK_CAP)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        smoke=args.smoke,
        hilbert_basis_path=args.hilbert_basis,
        orbit_index=args.orbit,
        branch=args.branch,
        directions=parse_integer_tuple(args.directions, "directions"),
        mod7_coordinates=parse_integer_tuple(args.mod7_coordinates, "mod7 coordinates"),
        state_cap=args.state_cap,
        pair_chunk_cap=args.pair_chunk_cap,
    )
    torsion.pointed.atomic_write(args.output, result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.output),
                "selected_directions": result["configuration"]["selected_directions"],
                "selected_directions_completed_through_grade_eight": result[
                    "selected_directions_completed_through_grade_eight"
                ],
                "selected_directions_with_explicit_cap_skip": result[
                    "selected_directions_with_explicit_cap_skip"
                ],
                "projected_group": result["configuration"]["projected_group"],
                "elapsed_seconds": result["elapsed_seconds"],
            },
            indent=2,
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
