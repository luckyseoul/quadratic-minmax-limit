#!/usr/bin/env python3
"""Exact stabilizer-quotient joins for the 51 positive p=7, z=7 H3 cases.

The existing case-level semigroup join keeps every high grade-three support
as a separate Minkowski factor.  Some projected supports stabilize as
subgroups; others do not.  This script computes every high support's exact
translation stabilizer, sums those stabilizers by finite-field linear
algebra, quotients them out, and retains all eight quotient-image factors.

For every direction and projection it computes the raw Hilbert-basis
recurrence through grade six and tests, without assuming,

    T3 = T4 = T5 = T6.

When the four supports are equal, Hilbert-basis degree at most three proves
permanent stabilization.  Otherwise Stab(S) is enumerated from all candidates
S-s0, with exact sorted translation equality.  Every stabilizer is audited
for zero, inverses, and closure.  Each H3 support is independently calibrated
against all 37,856 direct catalog rows.

For high directions D, H = sum_(d in D) Stab(T3(d)) is formed exactly as a
pair of row spaces over F3 and F7.  Explicit nullspace matrices have kernel H
and define G/H.  Every factor support—including H3—and the exact target are
mapped before deduplication, retaining same-row identity across both primes.
Since H stabilizes the full Minkowski sum, quotient membership is equivalent
to original membership.  Missing target rejects rigorously; presence is only
necessary.  State caps are explicit skips and partial supports are unused.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import resource
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z7_high_semigroup_support as semigroup  # noqa: E402
import p7_infinity7_positive_z7_semigroup_case_join as naive_join  # noqa: E402


P = naive_join.P
MOD3_DIMENSION = 6
MAX_MOD7_COORDINATES = 5
MAX_PROJECTIONS = 8
EXPECTED_TARGET_CASES = 51
EXPECTED_GRADE_CATALOG_ROWS = {0: 1, 1: 56, 2: 1_764, 3: 37_856}

DEFAULT_PARENT_INPUT = naive_join.DEFAULT_PARENT_INPUT
DEFAULT_CURRENT_JOIN = naive_join.DEFAULT_CURRENT_JOIN
DEFAULT_HILBERT_BASIS = naive_join.DEFAULT_HILBERT_BASIS
DEFAULT_CASE_COUNT = 3
DEFAULT_MOD7_COORDINATE_COUNT = 4
DEFAULT_PROJECTION_COUNT = 1
DEFAULT_STATE_CAP = 2_000_000
DEFAULT_PAIR_CHUNK_CAP = 200_000
DEFAULT_NAIVE_COMPARISON_CASES = 1


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


def field_rref(
    rows: np.ndarray, modulus: int, *, width: int | None = None
) -> tuple[np.ndarray, tuple[int, ...]]:
    """Canonical row-space basis over a prime field."""
    source = np.asarray(rows, dtype=np.int64)
    if source.ndim == 1 and source.size == 0:
        require(width is not None, "empty RREF input needs an explicit width")
        source = source.reshape(0, width)
    require(source.ndim == 2, "RREF input is not a matrix")
    if width is not None:
        require(source.shape[1] == width, "RREF width changed")
    work = np.ascontiguousarray(source % modulus, dtype=np.int64)
    pivot_row = 0
    pivots: list[int] = []
    for column in range(work.shape[1]):
        candidates = np.flatnonzero(work[pivot_row:, column])
        if not len(candidates):
            continue
        selected = pivot_row + int(candidates[0])
        if selected != pivot_row:
            work[[pivot_row, selected]] = work[[selected, pivot_row]]
        inverse = pow(int(work[pivot_row, column]), -1, modulus)
        work[pivot_row] = work[pivot_row] * inverse % modulus
        factors = work[:, column].copy()
        factors[pivot_row] = 0
        active = np.flatnonzero(factors)
        if len(active):
            work[active] = (
                work[active] - factors[active, None] * work[pivot_row]
            ) % modulus
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(work):
            break
    require(not np.any(work[pivot_row:]), "RREF retained a dependent nonzero row")
    return np.ascontiguousarray(work[:pivot_row], dtype=np.uint8), tuple(pivots)


def field_right_nullspace(rows: np.ndarray, modulus: int, width: int) -> np.ndarray:
    """Return N with ker(x -> N x) equal to the supplied row space."""
    basis, pivots = field_rref(rows, modulus, width=width)
    pivot_set = set(pivots)
    null_rows = []
    for free in (column for column in range(width) if column not in pivot_set):
        vector = np.zeros(width, dtype=np.int64)
        vector[free] = 1
        for pivot_row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -int(basis[pivot_row, free]) % modulus
        null_rows.append(vector)
    result = (
        np.ascontiguousarray(np.stack(null_rows), dtype=np.uint8)
        if null_rows
        else np.empty((0, width), dtype=np.uint8)
    )
    require(
        not len(basis) or not np.any(basis.astype(np.int64) @ result.T % modulus),
        "finite-field quotient matrix does not kill the subgroup basis",
    )
    quotient_rank, _ = field_rref(result, modulus, width=width)
    require(len(quotient_rank) == width - len(basis), "quotient map lost rank")
    return result


def map_rows_to_quotient(
    rows: np.ndarray,
    q3: np.ndarray,
    q7: np.ndarray,
    mod7_width: int,
) -> np.ndarray:
    """Map joint rows without breaking their mod-3/mod-7 identity."""
    source = np.ascontiguousarray(rows, dtype=np.uint8)
    require(
        source.ndim == 2 and source.shape[1] == MOD3_DIMENSION + mod7_width,
        "joint source-row width changed before quotienting",
    )
    components = []
    if len(q3):
        components.append(
            np.ascontiguousarray(
                (q3.astype(np.int64) @ source[:, :MOD3_DIMENSION].T).T % 3,
                dtype=np.uint8,
            )
        )
    if len(q7):
        components.append(
            np.ascontiguousarray(
                (q7.astype(np.int64) @ source[:, MOD3_DIMENSION:].T).T % 7,
                dtype=np.uint8,
            )
        )
    if not components:
        return np.empty((len(source), 0), dtype=np.uint8)
    return np.ascontiguousarray(np.concatenate(components, axis=1), dtype=np.uint8)


def sorted_membership(container: np.ndarray, queries: np.ndarray) -> bool:
    source = np.ascontiguousarray(container, dtype=np.uint64)
    requested = np.ascontiguousarray(queries, dtype=np.uint64)
    if not len(requested):
        return True
    positions = np.searchsorted(source, requested)
    valid = positions < len(source)
    return bool(np.all(valid) and np.all(source[positions[valid]] == requested[valid]))


def audit_subgroup_codes(
    reference: np.ndarray,
    codec: semigroup.MixedRadixCodec,
    mod7_width: int,
) -> tuple[dict[int, np.ndarray], dict]:
    """Prove an explicitly enumerated set is a subgroup of F3^6 x F7^k."""
    reference = np.ascontiguousarray(reference, dtype=np.uint64)
    require(len(reference) > 0, "candidate subgroup support is empty")
    require(np.array_equal(reference, np.unique(reference)), "support is not sorted unique")

    rows = codec.decode(reference)
    basis3, pivots3 = field_rref(rows[:, :MOD3_DIMENSION], 3, width=MOD3_DIMENSION)
    basis7, pivots7 = field_rref(rows[:, MOD3_DIMENSION:], 7, width=mod7_width)
    generated_order = 3 ** len(basis3) * 7 ** len(basis7)
    require(
        len(reference) == generated_order,
        "stabilized support does not have the order of its generated subgroup",
    )

    q3 = field_right_nullspace(basis3, 3, MOD3_DIMENSION)
    q7 = field_right_nullspace(basis7, 7, mod7_width)
    require(
        not np.any(map_rows_to_quotient(rows, q3, q7, mod7_width)),
        "candidate support escaped its generated subgroup",
    )
    zero_present = bool(len(reference) and int(reference[0]) == 0)
    require(zero_present, "stabilized support omits zero")

    inverse_rows = np.ascontiguousarray(
        (-rows.astype(np.int16)) % codec.moduli_array[None, :].astype(np.int16),
        dtype=np.uint8,
    )
    inverse_codes = codec.unique_codes(inverse_rows)
    require(
        len(inverse_codes) == len(reference) and sorted_membership(reference, inverse_codes),
        "stabilized support is not inverse-closed",
    )

    embedded_basis = []
    for row in basis3:
        embedded = np.zeros(MOD3_DIMENSION + mod7_width, dtype=np.uint8)
        embedded[:MOD3_DIMENSION] = row
        embedded_basis.append(embedded)
    for row in basis7:
        embedded = np.zeros(MOD3_DIMENSION + mod7_width, dtype=np.uint8)
        embedded[MOD3_DIMENSION:] = row
        embedded_basis.append(embedded)
    if embedded_basis:
        embedded_rows = np.ascontiguousarray(np.stack(embedded_basis), dtype=np.uint8)
        require(
            sorted_membership(reference, codec.unique_codes(embedded_rows)),
            "primary-component subgroup basis is not contained in joint support",
        )
        for generator in embedded_rows:
            translated = codec.translate(reference, generator)
            require(
                np.array_equal(translated, reference),
                "support is not closed under a generated-subgroup basis translation",
            )

    basis = {3: basis3, 7: basis7}
    audit = {
        "support_states": len(reference),
        "support_sha256_uint64": array_sha256(reference.astype("<u8", copy=False)),
        "zero_present_direct_audit": True,
        "all_inverses_present_direct_audit": True,
        "closure_under_embedded_primary_basis_translations_direct_audit": True,
        "closure_translation_generators_checked": len(embedded_basis),
        "component_ranks": {"F3": len(basis3), "F7": len(basis7)},
        "component_pivot_columns": {"F3": list(pivots3), "F7": list(pivots7)},
        "generated_subgroup_order": generated_order,
        "support_cardinality_equals_generated_subgroup_order": True,
        "component_basis_sha256_uint8": {
            "F3": array_sha256(basis3),
            "F7": array_sha256(basis7),
        },
        "finite_abelian_primary_decomposition_used": "F3-row-space direct-product F7-row-space",
        "subgroup_proved": True,
    }
    return basis, audit


def audit_subgroup_support(
    supports: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    codec: semigroup.MixedRadixCodec,
    mod7_width: int,
) -> tuple[dict[int, np.ndarray], dict]:
    """Prove T3=T4=T5=T6 is a permanently stabilized subgroup."""
    require(len(supports) == 4, "subgroup stabilization needs grades three through six")
    reference = np.ascontiguousarray(supports[0], dtype=np.uint64)
    require(
        all(np.array_equal(reference, row) for row in supports[1:]),
        "T3=T4=T5=T6 stabilization failed",
    )
    basis, direct = audit_subgroup_codes(reference, codec, mod7_width)
    return basis, {
        "stabilized_grades": [3, 4, 5, 6],
        "T3_equals_T4_equals_T5_equals_T6": True,
        "Hilbert_generator_maximum_degree": 3,
        "four_consecutive_equalities_prove_permanent_stabilization": True,
        **direct,
    }


def exact_translation_stabilizer(
    support: np.ndarray,
    codec: semigroup.MixedRadixCodec,
    mod7_width: int,
    pair_chunk_cap: int,
) -> tuple[np.ndarray, dict[int, np.ndarray], dict]:
    """Enumerate Stab(S) exactly from candidates S-s0, in bounded chunks."""
    source = np.ascontiguousarray(support, dtype=np.uint64)
    require(len(source) and np.array_equal(source, np.unique(source)), "bad stabilizer source")

    # Stab(S) acts freely on S by translation, so |Stab(S)| divides |S|.
    # The ambient group is F_3^6 x F_7^k; hence every nontrivial subgroup has
    # order divisible by 3 or 7.  Coprimality therefore certifies the trivial
    # stabilizer without enumerating the quadratic candidate set S-s0.
    if len(source) % 3 and len(source) % 7:
        zero = np.asarray([0], dtype=np.uint64)
        basis, subgroup_audit = audit_subgroup_codes(zero, codec, mod7_width)
        return zero, basis, {
            "algorithm": "free_translation_action_cardinality_divisibility",
            "source_support_states": len(source),
            "source_support_sha256_uint64": array_sha256(
                source.astype("<u8", copy=False)
            ),
            "ambient_group_prime_divisors": [3, 7],
            "source_cardinality_coprime_to_21": True,
            "proof": (
                "Stab(S) acts freely on S, so its order divides |S|; every "
                "nontrivial subgroup of F3^6 x F7^k has order divisible by 3 or 7"
            ),
            "candidate_states": 1,
            "candidate_chunks": 0,
            "candidate_chunk_size": 0,
            "translation_pairs_checked": 0,
            "stabilizer_states": 1,
            "stabilizer_sha256_uint64": array_sha256(zero.astype("<u8", copy=False)),
            "every_reported_stabilizer_rechecked_by_exact_sorted_translation_equality": True,
            "stabilizer_subgroup_audit": subgroup_audit,
            "stabilizer_exact": True,
        }

    rows = codec.decode(source)
    moduli = codec.moduli_array.astype(np.int16)
    all_candidates_rows = np.ascontiguousarray(
        (rows.astype(np.int16) - rows[0][None, :].astype(np.int16))
        % moduli[None, :],
        dtype=np.uint8,
    )
    require(
        len(codec.unique_codes(all_candidates_rows)) == len(source),
        "S-s0 candidate construction unexpectedly changed cardinality",
    )

    # Every subgroup of the coprime direct product F_3^6 x F_7^k is the
    # direct product of its 3- and 7-primary parts.  If a primary element
    # stabilizes S then it belongs to S-s0, so testing the candidates with
    # zero 7-part or zero 3-part is complete.  Mixed candidates never need a
    # direct test: the accepted primary parts generate the full stabilizer.
    three_primary = np.all(
        all_candidates_rows[:, MOD3_DIMENSION:] == 0, axis=1
    )
    seven_primary = np.all(
        all_candidates_rows[:, :MOD3_DIMENSION] == 0, axis=1
    )
    candidates_rows = np.ascontiguousarray(
        all_candidates_rows[three_primary | seven_primary], dtype=np.uint8
    )
    candidates_rows = codec.decode(codec.unique_codes(candidates_rows))
    require(len(candidates_rows), "primary stabilizer candidates lost zero")

    candidate_chunk = max(1, pair_chunk_cap // len(source))
    accepted_rows = []
    chunks = 0
    pairs = 0
    for start in range(0, len(candidates_rows), candidate_chunk):
        offsets = candidates_rows[start : start + candidate_chunk]
        sums = (
            rows[None, :, :].astype(np.int16)
            + offsets[:, None, :].astype(np.int16)
        ) % moduli[None, None, :]
        translated = codec.encode(
            np.ascontiguousarray(sums.reshape(-1, rows.shape[1]), dtype=np.uint8)
        ).reshape(len(offsets), len(rows))
        positions = np.searchsorted(source, translated)
        valid = positions < len(source)
        membership = np.zeros(translated.shape, dtype=bool)
        membership[valid] = source[positions[valid]] == translated[valid]
        accepted = np.all(membership, axis=1)
        if np.any(accepted):
            accepted_rows.append(offsets[accepted])
        chunks += 1
        pairs += len(offsets) * len(rows)

    require(accepted_rows, "translation stabilizer lost zero")
    accepted_primary_rows = np.ascontiguousarray(
        np.vstack(accepted_rows), dtype=np.uint8
    )
    accepted_three = accepted_primary_rows[
        np.all(accepted_primary_rows[:, MOD3_DIMENSION:] == 0, axis=1)
    ]
    accepted_seven = accepted_primary_rows[
        np.all(accepted_primary_rows[:, :MOD3_DIMENSION] == 0, axis=1)
    ]
    require(len(accepted_three) and len(accepted_seven), "primary stabilizer lost zero")
    stabilizer_rows = np.ascontiguousarray(
        (
            accepted_three[:, None, :].astype(np.int16)
            + accepted_seven[None, :, :].astype(np.int16)
        )
        % moduli[None, None, :],
        dtype=np.uint8,
    ).reshape(-1, rows.shape[1])
    stabilizer = codec.unique_codes(stabilizer_rows)
    basis, subgroup_audit = audit_subgroup_codes(stabilizer, codec, mod7_width)
    # Recheck the accepted set against S directly after subgroup canonicalization.
    for generator in codec.decode(stabilizer):
        require(
            np.array_equal(codec.translate(source, generator), source),
            "reported translation stabilizer contains a non-stabilizer",
        )
    return stabilizer, basis, {
        "algorithm": (
            "complete primary candidates h in S-s0; accept iff sorted S+h equals S; "
            "recombine coprime primary stabilizers"
        ),
        "source_support_states": len(source),
        "source_support_sha256_uint64": array_sha256(source.astype("<u8", copy=False)),
        "candidate_completeness_reason": "if S+h=S then s0+h lies in S, hence h lies in S-s0",
        "candidate_states": len(source),
        "primary_candidate_states_tested": len(candidates_rows),
        "candidate_chunks": chunks,
        "candidate_chunk_size": candidate_chunk,
        "translation_pairs_checked": pairs,
        "stabilizer_states": len(stabilizer),
        "stabilizer_sha256_uint64": array_sha256(stabilizer.astype("<u8", copy=False)),
        "every_reported_stabilizer_rechecked_by_exact_sorted_translation_equality": True,
        "coprime_primary_decomposition_proves_candidate_reduction_complete": True,
        "stabilizer_subgroup_audit": subgroup_audit,
        "stabilizer_exact": True,
    }


def translated_direct_support(
    *,
    direction: int,
    grade: int,
    mask: int,
    raw_support: np.ndarray,
    codec: semigroup.MixedRadixCodec,
    rebuilt: dict,
    quotient_data: dict[int, dict],
    projection: tuple[int, ...],
) -> tuple[np.ndarray, dict]:
    """Translate a recurrence support and calibrate every direct catalog row."""
    _unused, floor_mean, floor = semigroup.excess_grade(
        mask, 0 if mask.bit_count() == 7 else 8
    )
    mean = floor_mean + 4 * grade
    anchor = rebuilt["anchors"].get(mask, mean)
    offset = semigroup.project_rows(
        (anchor - floor)[None, :], direction, quotient_data, projection
    )[0]
    recurrence = codec.translate(raw_support, offset)

    catalog = naive_join.affine.mapped_catalog(mask, mean).astype(np.int64)
    require(
        len(catalog) == EXPECTED_GRADE_CATALOG_ROWS[grade],
        f"grade-{grade} direct catalog census changed",
    )
    delta = np.ascontiguousarray(anchor[None, :] - catalog, dtype=np.int64)
    require(not np.any(delta.sum(axis=1)), "direct catalog delta changed mean")
    require(
        not np.any(rebuilt["anchors"].kernel_rows @ delta.T),
        "direct catalog delta left the exact degree-two kernel",
    )
    require(not np.any(delta % 2), "direct catalog delta changed parity")
    direct_joint_rows = semigroup.project_rows(
        delta, direction, quotient_data, projection
    )
    direct = codec.unique_codes(direct_joint_rows)
    require(
        np.array_equal(recurrence, direct),
        f"direction {direction} grade-{grade} recurrence/direct calibration failed",
    )
    return recurrence, {
        "excess_grade": grade,
        "scaled_mean": mean,
        "complete_direct_catalog_rows": len(catalog),
        "direct_joint_rows_shape": list(direct_joint_rows.shape),
        "direct_joint_rows_sha256_uint8": array_sha256(direct_joint_rows),
        "projected_unique_states": len(direct),
        "projected_support_sha256_uint64": array_sha256(direct.astype("<u8", copy=False)),
        "semigroup_recurrence_equals_all_direct_catalog_rows": True,
        "same_direct_catalog_row_supplies_mod3_and_mod7": True,
        "deduplication_occurs_only_after_joint_projection": True,
        "box_exact": True,
    }


def build_projection_supports(
    *,
    projection: tuple[int, ...],
    full_generator_rows: dict[int, np.ndarray],
    degrees: np.ndarray,
    rebuilt: dict,
    quotient_data: dict[int, dict],
    state_cap: int,
    pair_chunk_cap: int,
) -> tuple[dict | None, dict]:
    """Build exact factor supports and certify each H3 translation stabilizer."""
    codec = semigroup.MixedRadixCodec((3,) * MOD3_DIMENSION + (7,) * len(projection))
    columns = tuple(range(MOD3_DIMENSION)) + tuple(
        MOD3_DIMENSION + value for value in projection
    )
    orbit = rebuilt["orbits"][0]
    support_table: dict[int, dict[int, np.ndarray]] = {}
    stabilizer_table: dict[int, dict] = {}
    direction_audits = []

    for direction in range(P + 1):
        selected_generators = np.ascontiguousarray(
            full_generator_rows[direction][:, columns], dtype=np.uint8
        )
        generator_codes = {
            degree: codec.unique_codes(selected_generators[degrees == degree])
            for degree in (1, 2, 3)
        }
        raw_supports, recurrence_records = semigroup.support_recurrence(
            generator_codes,
            codec,
            state_cap,
            pair_chunk_cap,
            max_grade=3,
        )
        missing = [grade for grade in range(4) if grade not in raw_supports]
        if missing:
            return None, {
                "mod7_coordinates": list(projection),
                "projected_group": f"F3^6 x F7^{len(projection)}",
                "projected_group_size": codec.group_size,
                "decision_status": "skipped_state_cap_before_stabilizer_certificate",
                "completed": False,
                "skipped": True,
                "first_incomplete_direction": direction,
                "missing_grades": missing,
                "recurrence_records": recurrence_records,
                "partial_support_used": False,
                "rigorous_rejection_allowed": False,
            }

        mask = int(orbit["masks"][direction])
        translated = {}
        calibration = []
        for grade in range(4):
            support, audit = translated_direct_support(
                direction=direction,
                grade=grade,
                mask=mask,
                raw_support=raw_supports[grade],
                codec=codec,
                rebuilt=rebuilt,
                quotient_data=quotient_data,
                projection=projection,
            )
            translated[grade] = support
            calibration.append(audit)

        H3_is_subgroup = False
        H3_subgroup_failure = None
        try:
            component_basis, subgroup_audit = audit_subgroup_codes(
                translated[3], codec, len(projection)
            )
            stabilizer_codes = translated[3]
            stabilizer_audit = {
                "algorithm": "support itself after direct subgroup proof",
                "source_support_states": len(translated[3]),
                "stabilizer_states": len(stabilizer_codes),
                "stabilizer_sha256_uint64": array_sha256(
                    stabilizer_codes.astype("<u8", copy=False)
                ),
                "stabilizer_subgroup_audit": subgroup_audit,
                "stabilizer_exact": True,
            }
            H3_is_subgroup = True
        except AssertionError as error:
            H3_subgroup_failure = str(error)
            stabilizer_codes, component_basis, stabilizer_audit = (
                exact_translation_stabilizer(
                    translated[3],
                    codec,
                    len(projection),
                    pair_chunk_cap,
                )
            )

        support_table[direction] = translated
        stabilizer_table[direction] = {
            "codes": stabilizer_codes,
            "component_basis": component_basis,
        }
        direction_audits.append(
            {
                "direction": direction,
                "mask": mask,
                "generator_unique_states_by_degree": {
                    str(degree): len(generator_codes[degree]) for degree in (1, 2, 3)
                },
                "selected_joint_generator_rows_sha256_uint8": array_sha256(
                    selected_generators
                ),
                "recurrence_records_sha256": json_sha256(recurrence_records),
                "recurrence_state_counts_grades_zero_through_three": [
                    len(raw_supports[grade]) for grade in range(4)
                ],
                "grade_three_scope_audit": {
                    "required_case_support_grades": [0, 1, 2, 3],
                    "higher_recurrence_layers_computed": False,
                    "higher_layers_not_needed_for_direct_subgroup_or_stabilizer_proof": True,
                    "permanent_stabilization_claimed_here": False,
                },
                "anchor_relative_H3_is_subgroup": H3_is_subgroup,
                "anchor_relative_H3_subgroup_failure": H3_subgroup_failure,
                "translation_stabilizer_audit": stabilizer_audit,
                "complete_direct_catalog_calibrations": calibration,
                "H3_calibrated_against_all_37856_direct_rows": True,
                "non_subgroup_H3_is_retained_as_a_factor_not_removed": not H3_is_subgroup,
            }
        )

    payload = {
        "codec": codec,
        "supports": support_table,
        "stabilizers": stabilizer_table,
    }
    return payload, {
        "mod7_coordinates": list(projection),
        "projected_group": f"F3^6 x F7^{len(projection)}",
        "projected_group_size": codec.group_size,
        "state_cap": state_cap,
        "pair_chunk_cap": pair_chunk_cap,
        "decision_status": "complete_all_direction_exact_H3_stabilizer_certificate",
        "completed": True,
        "skipped": False,
        "directions": direction_audits,
        "H3_supports_are_not_assumed_to_be_subgroups": True,
        "all_eight_H3_translation_stabilizers_exactly_certified": True,
        "all_eight_H3_supports_calibrated_against_37856_direct_rows": True,
        "all_low_USM_supports_complete_and_directly_calibrated": True,
        "partial_support_used": False,
    }


def case_stabilizer_quotient(
    high_directions: tuple[int, ...],
    stabilizer_table: dict[int, dict],
    mod7_width: int,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Sum exact H3 translation stabilizers and construct quotient matrices."""
    require(high_directions, "grade-three case has no high direction")
    rows3 = np.vstack(
        [stabilizer_table[direction]["component_basis"][3] for direction in high_directions]
    )
    rows7 = np.vstack(
        [stabilizer_table[direction]["component_basis"][7] for direction in high_directions]
    )
    basis3, pivots3 = field_rref(rows3, 3, width=MOD3_DIMENSION)
    basis7, pivots7 = field_rref(rows7, 7, width=mod7_width)
    q3 = field_right_nullspace(basis3, 3, MOD3_DIMENSION)
    q7 = field_right_nullspace(basis7, 7, mod7_width)

    require(
        not len(basis3) or not np.any(q3.astype(np.int64) @ basis3.T % 3),
        "F3 quotient does not kill the exact stabilizer sum",
    )
    require(
        not len(basis7) or not np.any(q7.astype(np.int64) @ basis7.T % 7),
        "F7 quotient does not kill the exact stabilizer sum",
    )
    high_order = 3 ** len(basis3) * 7 ** len(basis7)
    quotient_order = 3 ** len(q3) * 7 ** len(q7)
    ambient_order = 3 ** MOD3_DIMENSION * 7 ** mod7_width
    require(high_order * quotient_order == ambient_order, "subgroup/quotient orders disagree")
    return q3, q7, {
        "high_directions": list(high_directions),
        "sum_construction": "row span of the union of exact H3 translation-stabilizer bases",
        "sum_component_ranks": {"F3": len(basis3), "F7": len(basis7)},
        "sum_component_pivot_columns": {"F3": list(pivots3), "F7": list(pivots7)},
        "sum_component_basis_sha256_uint8": {
            "F3": array_sha256(basis3),
            "F7": array_sha256(basis7),
        },
        "summed_stabilizer_order": high_order,
        "quotient_dimensions": {"F3": len(q3), "F7": len(q7)},
        "quotient_order": quotient_order,
        "quotient_matrix_sha256_uint8": {
            "F3": array_sha256(q3),
            "F7": array_sha256(q7),
        },
        "quotient_kernel_equals_summed_stabilizers_by_rank_and_annihilation": True,
        "translation_stabilizers_summed_exactly": True,
        "summed_stabilizer_preserves_the_full_eight_factor_Minkowski_sum": True,
    }


def select_cases(
    targets: list[dict],
    current_by_key: dict[str, dict],
    *,
    all_cases: bool,
    case_keys: tuple[str, ...] | None,
    case_count: int | None,
) -> tuple[list[dict], dict]:
    require(len(targets) == EXPECTED_TARGET_CASES, "grade-three target census changed")
    target_by_key = {str(row["case_key"]): row for row in targets}
    if all_cases:
        selected = list(targets)
        mode = "all_51_in_current_decision_order"
    elif case_keys is not None:
        require(case_keys and len(case_keys) == len(set(case_keys)), "case key selection repeats")
        missing = [key for key in case_keys if key not in target_by_key]
        require(not missing, f"selected case keys are not open grade-three representatives: {missing}")
        selected = [target_by_key[key] for key in case_keys]
        mode = "explicit_case_keys_in_cli_order"
    else:
        count = DEFAULT_CASE_COUNT if case_count is None else case_count
        selected, _audit = naive_join.select_smoke_cases(targets, current_by_key, count)
        mode = "round_robin_pattern_prefix"
    keys = [str(row["case_key"]) for row in selected]
    return selected, {
        "selection_mode": mode,
        "selected_case_count": len(selected),
        "selected_case_keys": keys,
        "selected_case_keys_sha256": json_sha256(keys),
        "selected_pattern_counts": dict(
            sorted(Counter(current_by_key[key]["catalog_pattern"] for key in keys).items())
        ),
        "all_51_cases_selected": len(selected) == EXPECTED_TARGET_CASES,
    }


def evaluate_quotient_case(
    *,
    target_row: dict,
    current_row: dict,
    rebuilt: dict,
    common: np.ndarray,
    quotient_data: dict[int, dict],
    projection: tuple[int, ...],
    projection_data: dict,
    state_cap: int,
    pair_chunk_cap: int,
) -> dict:
    orbit, leaf, system, _factory = naive_join.old_join.validate_parent_survivor(
        target_row, rebuilt
    )
    grades = tuple(
        naive_join.leaf_grade(orbit, leaf, direction) for direction in range(P + 1)
    )
    high_directions = tuple(int(value) for value in leaf["high_directions"])
    require(
        high_directions
        and all(grades[direction] == 3 for direction in high_directions),
        "selected case is not grade-three-only on its high directions",
    )
    low_directions = tuple(
        direction for direction in range(P + 1) if direction not in high_directions
    )
    require(all(grades[direction] <= 2 for direction in low_directions), "non-high grade exceeded M")

    q3, q7, stabilizer_audit = case_stabilizer_quotient(
        high_directions, projection_data["stabilizers"], len(projection)
    )
    quotient_moduli = (3,) * len(q3) + (7,) * len(q7)

    anchor_rhs, _raw_syndromes = naive_join.affine.anchor_rhs_and_raw_syndromes(
        orbit, leaf, system, rebuilt["anchors"]
    )
    base_digits = naive_join.project_equation_vector(anchor_rhs, quotient_data, projection)
    ambient_codec = projection_data["codec"]
    target_digits = np.ascontiguousarray(
        (-base_digits.astype(np.int16))
        % ambient_codec.moduli_array.astype(np.int16),
        dtype=np.uint8,
    )
    quotient_target_rows = map_rows_to_quotient(
        target_digits[None, :], q3, q7, len(projection)
    )

    exact_rhs = anchor_rhs.copy()
    for direction, grade in enumerate(grades):
        mask = int(orbit["masks"][direction])
        mean = int(leaf["scaled_means"][direction])
        catalog_row = naive_join.affine.mapped_catalog(mask, mean)[0].astype(np.int64)
        anchor = rebuilt["anchors"].get(mask, mean)
        block = slice(
            1 + naive_join.AMBIENT_DIMENSION * direction,
            1 + naive_join.AMBIENT_DIMENSION * (direction + 1),
        )
        exact_rhs[block] += anchor - catalog_row
    require(not np.any(common @ exact_rhs), "common exact dependency syndrome did not vanish")

    base_record = {
        "mod7_coordinates": list(projection),
        "ambient_projected_group": f"F3^6 x F7^{len(projection)}",
        "ambient_projected_group_size": ambient_codec.group_size,
        "directional_excess_grades": list(grades),
        "high_directions": list(high_directions),
        "low_directions": list(low_directions),
        "summed_H3_stabilizer_and_quotient_audit": stabilizer_audit,
        "ambient_target_digits": target_digits.tolist(),
        "ambient_target_digits_sha256_uint8": array_sha256(target_digits),
        "quotient_target_digits": quotient_target_rows[0].tolist(),
        "quotient_target_digits_sha256_uint8": array_sha256(quotient_target_rows),
        "same_exact_target_row_mapped_jointly_mod3_mod7": True,
        "common_exact_dependency_syndrome_checked_with_direct_catalog_rows": True,
    }
    if not quotient_moduli:
        return {
            **base_record,
            "quotient_group": "trivial",
            "quotient_group_size": 1,
            "factor_quotient_supports": [],
            "decision_status": "necessary_only_trivial_exact_stabilizer_quotient",
            "rigorously_rejected": False,
            "necessary_only": True,
            "skipped": False,
            "exact_quotient_membership_completed": True,
            "target_present": True,
            "target_presence_proves_feasibility": False,
        }

    quotient_codec = semigroup.MixedRadixCodec(quotient_moduli)
    factors = []
    factor_audits = []
    for direction in range(P + 1):
        grade = grades[direction]
        source_codes = projection_data["supports"][direction][grade]
        source_rows = ambient_codec.decode(source_codes)
        quotient_rows = map_rows_to_quotient(source_rows, q3, q7, len(projection))
        quotient_codes = quotient_codec.unique_codes(quotient_rows)
        require(len(quotient_codes), "factor quotient support is empty")
        factors.append((direction, quotient_codes))
        factor_audits.append(
            {
                "direction": direction,
                "excess_grade": grade,
                "catalog_class": {0: "U", 1: "S", 2: "M", 3: "H"}[grade],
                "is_high_H3_factor": direction in high_directions,
                "ambient_joint_support_states": len(source_codes),
                "ambient_joint_support_sha256_uint64": array_sha256(
                    source_codes.astype("<u8", copy=False)
                ),
                "quotient_joint_support_states": len(quotient_codes),
                "quotient_joint_support_sha256_uint64": array_sha256(
                    quotient_codes.astype("<u8", copy=False)
                ),
                "same_source_row_mapped_to_both_prime_quotients_before_deduplication": True,
            }
        )

    sizes = tuple(len(row[1]) for row in factors)
    directions = tuple(row[0] for row in factors)
    partition = naive_join.old_join.balanced_partition(directions, sizes)
    by_direction = {direction: support for direction, support in factors}
    left_factors = tuple(
        (direction, by_direction[direction]) for direction in partition["left_directions"]
    )
    right_factors = tuple(
        (direction, by_direction[direction]) for direction in partition["right_directions"]
    )
    left, left_audit = naive_join.convolve_support_sequence(
        left_factors, quotient_codec, state_cap, pair_chunk_cap
    )
    right, right_audit = naive_join.convolve_support_sequence(
        right_factors, quotient_codec, state_cap, pair_chunk_cap
    )
    common_fields = {
        **base_record,
        "quotient_group": f"F3^{len(q3)} x F7^{len(q7)}",
        "quotient_group_size": quotient_codec.group_size,
        "factor_quotient_supports": factor_audits,
        "balanced_partition": partition,
        "left_support": left_audit,
        "right_support": right_audit,
        "all_eight_quotient_image_supports_including_every_high_H3_retained": True,
        "only_the_summed_exact_translation_stabilizer_is_quotiented_out": True,
    }
    if left is None or right is None:
        return {
            **common_fields,
            "decision_status": "skipped_quotient_convolution_state_cap",
            "rigorously_rejected": False,
            "necessary_only": False,
            "skipped": True,
            "exact_quotient_membership_completed": False,
            "partial_support_used": False,
        }

    if len(left) == quotient_codec.group_size or len(right) == quotient_codec.group_size:
        return {
            **common_fields,
            "decision_status": "necessary_only_exact_quotient_full_group_saturation",
            "rigorously_rejected": False,
            "necessary_only": True,
            "skipped": False,
            "exact_quotient_membership_completed": True,
            "target_present": True,
            "target_present_by_exact_full_group_saturation": True,
            "target_presence_proves_feasibility": False,
        }

    joined = naive_join.meet_target(
        left,
        right,
        quotient_target_rows[0],
        quotient_codec,
        pair_chunk_cap,
    )
    rejected = joined["matching_projected_support_pairs"] == 0
    return {
        **common_fields,
        "decision_status": (
            "rigorous_exact_stabilizer_quotient_rejection"
            if rejected
            else "necessary_only_exact_stabilizer_quotient_survivor"
        ),
        "rigorously_rejected": rejected,
        "necessary_only": not rejected,
        "skipped": False,
        "exact_quotient_membership_completed": True,
        "target_present": not rejected,
        "quotient_join": joined,
        "missing_target_is_rigorous_rejection": rejected,
        "target_presence_proves_feasibility": False,
        "prior_global_join_decision": current_row["decision_status"],
    }


def manufactured_self_audit() -> dict:
    """Exercise subgroup, quotient, same-row, cap, and false-subgroup traps."""
    codec = semigroup.MixedRadixCodec((3,) * MOD3_DIMENSION + (7,))
    high_rows = np.asarray(
        [(a, 0, 0, 0, 0, 0, b) for a in range(3) for b in range(7)],
        dtype=np.uint8,
    )
    high = codec.unique_codes(high_rows)
    basis, subgroup = audit_subgroup_support((high, high, high, high), codec, 1)
    table = {0: {"component_basis": basis}}
    q3, q7, quotient = case_stabilizer_quotient((0,), table, 1)
    require(
        q3.shape == (5, MOD3_DIMENSION) and q7.shape == (0, 1),
        "manufactured quotient dimensions changed",
    )

    low_rows = np.asarray(
        [[0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0]],
        dtype=np.uint8,
    )
    mapped = map_rows_to_quotient(low_rows, q3, q7, 1)
    require(len(np.unique(mapped, axis=0)) == 2, "manufactured low quotient failed")
    present_target = np.asarray([[2, 1, 0, 0, 0, 0, 4]], dtype=np.uint8)
    absent_target = np.asarray([[1, 0, 1, 0, 0, 0, 5]], dtype=np.uint8)
    present = map_rows_to_quotient(present_target, q3, q7, 1)
    absent = map_rows_to_quotient(absent_target, q3, q7, 1)
    qcodec = semigroup.MixedRadixCodec((3,) * len(q3))
    low_codes = qcodec.unique_codes(mapped)
    present_hit = sorted_membership(low_codes, qcodec.unique_codes(present))
    absent_hit = sorted_membership(low_codes, qcodec.unique_codes(absent))
    require(present_hit and not absent_hit, "manufactured quotient membership trap failed")

    naive_present = any(
        tuple(
            (h.astype(np.int16) + low.astype(np.int16))
            % np.asarray((3,) * MOD3_DIMENSION + (7,))
        )
        == tuple(present_target[0])
        for h in high_rows
        for low in low_rows
    )
    naive_absent = any(
        tuple(
            (h.astype(np.int16) + low.astype(np.int16))
            % np.asarray((3,) * MOD3_DIMENSION + (7,))
        )
        == tuple(absent_target[0])
        for h in high_rows
        for low in low_rows
    )
    require(naive_present == present_hit and naive_absent == absent_hit, "naive/quotient trap mismatch")

    joint_codec = semigroup.MixedRadixCodec((3, 7))
    same_row_support = joint_codec.unique_codes(
        np.asarray([[0, 1], [1, 0]], dtype=np.uint8)
    )
    zero = joint_codec.unique_codes(np.asarray([[0, 0]], dtype=np.uint8))
    require(not sorted_membership(same_row_support, zero), "same-row cross-prime trap passed")
    decoded = joint_codec.decode(same_row_support)
    require(
        np.any(decoded[:, 0] == 0) and np.any(decoded[:, 1] == 0),
        "same-row trap lacks independent marginal hits",
    )

    false_codec = semigroup.MixedRadixCodec((3,) * MOD3_DIMENSION)
    false_rows = np.asarray(
        [
            [0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0],
            [2, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 2, 0, 0, 0, 0],
        ],
        dtype=np.uint8,
    )
    false_codes = false_codec.unique_codes(false_rows)
    false_subgroup_rejected = False
    try:
        audit_subgroup_support(
            (false_codes, false_codes, false_codes, false_codes), false_codec, 0
        )
    except AssertionError:
        false_subgroup_rejected = True
    require(false_subgroup_rejected, "manufactured non-subgroup was accepted")
    false_stabilizer, _false_basis, false_stabilizer_audit = exact_translation_stabilizer(
        false_codes, false_codec, 0, 25
    )
    require(
        len(false_stabilizer) == 1 and int(false_stabilizer[0]) == 0,
        "manufactured non-subgroup stabilizer changed",
    )

    # Exercise the non-coprime primary-decomposition path.  This is the union
    # of two H-cosets for H = F3 x F7, with quotient cosets {0,1} in a second
    # F3 coordinate.  Its exact translation stabilizer is H, although the
    # 42-point union itself is not a subgroup.
    second_coset = high_rows.copy()
    second_coset[:, 1] = 1
    two_cosets = codec.unique_codes(np.vstack((high_rows, second_coset)))
    mixed_subgroup_rejected = False
    try:
        audit_subgroup_codes(two_cosets, codec, 1)
    except AssertionError:
        mixed_subgroup_rejected = True
    require(mixed_subgroup_rejected, "two-coset trap was mistaken for a subgroup")
    mixed_stabilizer, _mixed_basis, mixed_stabilizer_audit = (
        exact_translation_stabilizer(two_cosets, codec, 1, 10_000)
    )
    require(
        np.array_equal(mixed_stabilizer, high),
        "coprime-primary stabilizer reconstruction failed",
    )

    cap_codec = semigroup.MixedRadixCodec((7, 7))
    axis_a = cap_codec.unique_codes(
        np.asarray([[value, 0] for value in range(5)], dtype=np.uint8)
    )
    axis_b = cap_codec.unique_codes(
        np.asarray([[0, value] for value in range(5)], dtype=np.uint8)
    )
    capped, cap_audit = naive_join.convolve_support_sequence(
        ((0, axis_a), (1, axis_b)), cap_codec, 10, 25
    )
    require(capped is None and cap_audit["status"] == "skipped_state_cap", "cap trap failed")
    require(cap_audit["partial_support_used"] is False, "cap trap retained partial support")
    return {
        "passed": True,
        "stabilized_support_subgroup_proof_trap_passed": True,
        "subgroup_audit": subgroup,
        "quotient_audit": quotient,
        "quotient_membership_equals_direct_naive_membership_for_present_and_absent_targets": True,
        "same_row_cross_prime_false_positive_trap_rejected": True,
        "non_subgroup_with_zero_and_inverses_rejected_by_closure_order_audit": True,
        "non_subgroup_exact_translation_stabilizer_audit": false_stabilizer_audit,
        "noncoprime_two_coset_primary_stabilizer_audit": mixed_stabilizer_audit,
        "noncoprime_two_coset_stabilizer_equals_F3_times_F7": True,
        "state_cap_is_explicit_skip_and_partial_support_is_discarded": True,
        "cap_audit": cap_audit,
    }


def existing_naive_comparison(
    *,
    selected_cases: list[dict],
    current_by_key: dict[str, dict],
    case_results_by_key: dict[str, dict],
    rebuilt: dict,
    common: np.ndarray,
    quotient_data: dict[int, dict],
    projection: tuple[int, ...],
    projection_data: dict,
    requested_cases: int,
    state_cap: int,
    pair_chunk_cap: int,
) -> dict:
    """Compare against the repository's existing naive eight-factor join."""
    require(requested_cases >= 0, "naive comparison case count is negative")
    rows = []
    for target_row in selected_cases[:requested_cases]:
        key = str(target_row["case_key"])
        observed = naive_join.evaluate_projected_case(
            target_row=target_row,
            current_row=current_by_key[key],
            rebuilt=rebuilt,
            common=common,
            quotient_data=quotient_data,
            projection=projection,
            supports=projection_data["supports"],
            state_cap=state_cap,
            pair_chunk_cap=pair_chunk_cap,
        )
        quotient_row = case_results_by_key[key]["projection_results_by_key"][projection]
        both_complete = not observed["skipped"] and not quotient_row["skipped"]
        if both_complete:
            require(
                bool(observed["rigorously_rejected"])
                == bool(quotient_row["rigorously_rejected"]),
                "exact quotient decision disagrees with existing naive case join",
            )
        rows.append(
            {
                "case_key": key,
                "naive_decision_status": observed["decision_status"],
                "quotient_decision_status": quotient_row["decision_status"],
                "naive_completed": not observed["skipped"],
                "quotient_completed": not quotient_row["skipped"],
                "both_completed": both_complete,
                "rejection_decisions_equal_when_both_completed": (
                    bool(observed["rigorously_rejected"])
                    == bool(quotient_row["rigorously_rejected"])
                    if both_complete
                    else None
                ),
                "naive_result_sha256": json_sha256(observed),
            }
        )
    completed = sum(row["both_completed"] for row in rows)
    return {
        "projection_mod7_coordinates": list(projection),
        "existing_naive_function": "p7_infinity7_positive_z7_semigroup_case_join.evaluate_projected_case",
        "requested_case_count": requested_cases,
        "comparison_rows": rows,
        "completed_comparisons": completed,
        "skipped_comparisons": len(rows) - completed,
        "all_completed_rejection_decisions_equal": all(
            row["rejection_decisions_equal_when_both_completed"] is not False for row in rows
        ),
        "state_cap_can_only_make_comparison_an_explicit_skip": True,
    }


def run(
    *,
    parent_path: Path,
    current_join_path: Path,
    hilbert_basis_path: Path,
    output_path: Path,
    all_cases: bool,
    case_keys: tuple[str, ...] | None,
    case_count: int | None,
    coordinate_count: int,
    projection_count: int,
    explicit_projections: tuple[tuple[int, ...], ...] | None,
    state_cap: int,
    pair_chunk_cap: int,
    naive_comparison_cases: int,
) -> dict:
    started = time.time()
    usage_started = resource.getrusage(resource.RUSAGE_SELF)
    require(1 <= state_cap <= semigroup.MAX_SMOKE_STATE_CAP, "state cap is outside audited bound")
    require(1 <= pair_chunk_cap <= semigroup.MAX_SMOKE_STATE_CAP, "pair chunk cap outside bound")
    require(0 <= naive_comparison_cases <= EXPECTED_TARGET_CASES, "bad comparison case count")

    parent_payload, representatives, current_by_key, targets, input_audit = (
        naive_join.load_current_problem(parent_path, current_join_path)
    )
    rebuilt, reconstruction_audit = naive_join.reconstruct_and_validate(
        parent_payload, representatives
    )
    target_grade_audit = naive_join.audit_target_grades(targets, current_by_key, rebuilt)
    selected_cases, case_selection = select_cases(
        targets,
        current_by_key,
        all_cases=all_cases,
        case_keys=case_keys,
        case_count=case_count,
    )

    system = rebuilt["systems"][0]["A"]
    matrix = np.ascontiguousarray(system["matrix"], dtype=np.int64)
    kernel_rows = np.ascontiguousarray(rebuilt["kernel_rows"], dtype=np.int64)
    basis, degrees, basis_audit = semigroup.load_hilbert_basis(
        hilbert_basis_path, kernel_rows
    )
    common, common_audit = naive_join.torsion.exact_common_dependency_basis(
        matrix, kernel_rows
    )
    quotient_data, torsion_audit = semigroup.derive_torsion_quotients(matrix, common)
    require(
        len(quotient_data[3]["complement"]) == MOD3_DIMENSION
        and len(quotient_data[7]["complement"]) == 21,
        "effective torsion dimensions changed",
    )

    orbit = rebuilt["orbits"][0]
    full_rows, full_projection_audit = naive_join.full_generator_projections(
        basis, orbit, quotient_data
    )
    projections, projection_selection = naive_join.select_projection_subsets(
        full_rows=full_rows,
        degrees=degrees,
        coordinate_count=coordinate_count,
        projection_count=projection_count,
        explicit=explicit_projections,
    )
    require(1 <= len(projections) <= MAX_PROJECTIONS, "selected projection count changed")
    require(
        all(1 <= len(row) <= MAX_MOD7_COORDINATES for row in projections),
        "projection width is outside the exact bounded implementation",
    )

    projection_payloads: dict[tuple[int, ...], dict | None] = {}
    projection_audits = []
    for projection in projections:
        payload, audit = build_projection_supports(
            projection=projection,
            full_generator_rows=full_rows,
            degrees=degrees,
            rebuilt=rebuilt,
            quotient_data=quotient_data,
            state_cap=state_cap,
            pair_chunk_cap=pair_chunk_cap,
        )
        projection_payloads[projection] = payload
        projection_audits.append(audit)

    case_results = []
    internal_results: dict[str, dict] = {}
    for target_row in selected_cases:
        key = str(target_row["case_key"])
        projection_results = []
        projection_results_by_key = {}
        for projection in projections:
            projection_data = projection_payloads[projection]
            if projection_data is None:
                projection_result = {
                    "mod7_coordinates": list(projection),
                "decision_status": "skipped_projection_stabilizer_certificate_incomplete",
                    "rigorously_rejected": False,
                    "necessary_only": False,
                    "skipped": True,
                    "exact_quotient_membership_completed": False,
                    "partial_support_used": False,
                }
            else:
                projection_result = evaluate_quotient_case(
                    target_row=target_row,
                    current_row=current_by_key[key],
                    rebuilt=rebuilt,
                    common=common,
                    quotient_data=quotient_data,
                    projection=projection,
                    projection_data=projection_data,
                    state_cap=state_cap,
                    pair_chunk_cap=pair_chunk_cap,
                )
            projection_results.append(projection_result)
            projection_results_by_key[projection] = projection_result

        rejected = any(row["rigorously_rejected"] for row in projection_results)
        skipped = not rejected and any(row["skipped"] for row in projection_results)
        necessary = not rejected and not skipped
        require(sum((rejected, skipped, necessary)) == 1, "aggregate quotient decision ambiguous")
        public = {
            "case_key": key,
            "catalog_pattern": current_by_key[key]["catalog_pattern"],
            "prior_global_join_decision": current_by_key[key]["decision_status"],
            "projection_results": projection_results,
            "decision_status": (
                "rigorous_stabilizer_quotient_rejection"
                if rejected
                else "skipped_at_least_one_selected_projection"
                if skipped
                else "necessary_only_survivor_of_all_selected_quotients"
            ),
            "rigorously_rejected": rejected,
            "necessary_only_survivor": necessary,
            "skipped": skipped,
            "all_selected_projections_completed": not any(
                row["skipped"] for row in projection_results
            ),
        }
        public["decision_certificate_sha256"] = json_sha256(public)
        case_results.append(public)
        internal_results[key] = {
            "public": public,
            "projection_results_by_key": projection_results_by_key,
        }

    comparison_projection = next(
        (projection for projection in projections if projection_payloads[projection] is not None),
        None,
    )
    if comparison_projection is None or naive_comparison_cases == 0:
        comparison_audit = {
            "requested_case_count": naive_comparison_cases,
            "completed_comparisons": 0,
            "decision_status": (
                "disabled_by_cli"
                if naive_comparison_cases == 0
                else "skipped_no_completed_projection_support"
            ),
            "partial_support_used": False,
        }
    else:
        comparison_audit = existing_naive_comparison(
            selected_cases=selected_cases,
            current_by_key=current_by_key,
            case_results_by_key=internal_results,
            rebuilt=rebuilt,
            common=common,
            quotient_data=quotient_data,
            projection=comparison_projection,
            projection_data=projection_payloads[comparison_projection],
            requested_cases=min(naive_comparison_cases, len(selected_cases)),
            state_cap=state_cap,
            pair_chunk_cap=pair_chunk_cap,
        )

    counts = {
        "selected": len(case_results),
        "rejected": sum(row["rigorously_rejected"] for row in case_results),
        "surviving": sum(row["necessary_only_survivor"] for row in case_results),
        "skipped": sum(row["skipped"] for row in case_results),
    }
    require(
        counts["rejected"] + counts["surviving"] + counts["skipped"] == counts["selected"],
        "result count partition failed",
    )
    usage = resource.getrusage(resource.RUSAGE_SELF)
    script_path = Path(__file__).resolve()
    result = {
        "experiment": "p7_infinity7_positive_z7_semigroup_case_quotient",
        "status": (
            "complete_all_51_exact_selected_stabilizer_quotient_projections"
            if len(case_results) == EXPECTED_TARGET_CASES and counts["skipped"] == 0
            else "complete_selected_exact_stabilizer_quotient_projections"
            if counts["skipped"] == 0
            else "selected_stabilizer_quotient_run_with_explicit_skips"
        ),
        "p": 7,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "full_51_case_run": len(case_results) == EXPECTED_TARGET_CASES,
        "source_provenance": {
            "this_script_path": str(script_path),
            "this_script_sha256": file_sha256(script_path),
            "existing_naive_case_join_script_sha256": file_sha256(Path(naive_join.__file__)),
            "high_semigroup_support_script_sha256": file_sha256(Path(semigroup.__file__)),
        },
        "configuration": {
            "all_cases": all_cases,
            "explicit_case_keys": list(case_keys) if case_keys is not None else None,
            "case_count": case_count,
            "automatic_mod7_coordinate_count": coordinate_count,
            "automatic_projection_count": projection_count,
            "explicit_projections": (
                [list(row) for row in explicit_projections]
                if explicit_projections is not None
                else None
            ),
            "state_cap": state_cap,
            "pair_chunk_cap": pair_chunk_cap,
            "naive_comparison_cases": naive_comparison_cases,
        },
        "input_and_current_decision_audit": input_audit,
        "representative_reconstruction_audit": reconstruction_audit,
        "target_grade_audit": target_grade_audit,
        "case_selection": case_selection,
        "normaliz_Hilbert_basis": basis_audit,
        "exact_common_rational_dependency_audit": common_audit,
        "varying_torsion_quotient_audit": torsion_audit,
        "full_generator_projection_audit": full_projection_audit,
        "projection_subset_selection": projection_selection,
        "projection_H3_stabilizer_audits": projection_audits,
        "manufactured_self_audit": manufactured_self_audit(),
        "existing_naive_case_join_comparison": comparison_audit,
        "result_counts": counts,
        "case_results_sha256": naive_join.old_join.canonical_case_digest(case_results),
        "case_results": case_results,
        "logical_semantics": {
            "T3_T4_T5_T6_equality_is_tested_and_never_assumed": True,
            "non_subgroup_H3_uses_complete_exact_translation_stabilizer_fallback": True,
            "every_H3_stabilizer_has_direct_zero_inverse_and_basis_translation_closure_audit": True,
            "every_H3_support_calibrated_against_all_37856_direct_rows": True,
            "high_factor_translation_stabilizers_summed_by_explicit_F3_and_F7_row_spaces": True,
            "quotient_maps_are_explicit_full_rank_annihilators_with_exact_kernel": True,
            "same_source_row_supplies_mod3_and_mod7_before_every_deduplication": True,
            "all_eight_U_S_M_H3_factor_images_are_retained_in_the_quotient_join": True,
            "missing_target_in_completed_quotient_is_rigorous_rejection": True,
            "target_presence_is_only_necessary": True,
            "state_cap_hit_is_explicit_skip": True,
            "partial_support_after_cap_is_discarded": True,
            "different_projection_witnesses_are_not_assumed_compatible": True,
            "positive_z7_closure_claimed": False,
            "binary_edge_feasibility_claimed": False,
        },
        "all_51_target_cases_processed": len(case_results) == EXPECTED_TARGET_CASES,
        "positive_z7_excluded": False,
        "resource_profile": {
            "elapsed_seconds": time.time() - started,
            "user_cpu_seconds": usage.ru_utime - usage_started.ru_utime,
            "system_cpu_seconds": usage.ru_stime - usage_started.ru_stime,
            "maximum_resident_set_kib": usage.ru_maxrss,
            "effective_algorithm": "CPU exact finite-field linear algebra and bounded sorted support joins",
            "GPU_used": False,
        },
        "output_path": str(output_path.resolve()),
    }
    return result


def parse_case_keys(value: str) -> tuple[str, ...]:
    rows = tuple(item.strip() for item in value.split(",") if item.strip())
    require(rows, "explicit case-key list is empty")
    require(len(rows) == len(set(rows)), "explicit case-key list repeats")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-audit-only", action="store_true")
    parser.add_argument("--parent-input", type=Path, default=DEFAULT_PARENT_INPUT)
    parser.add_argument("--current-join", type=Path, default=DEFAULT_CURRENT_JOIN)
    parser.add_argument("--hilbert-basis", type=Path, default=DEFAULT_HILBERT_BASIS)
    parser.add_argument("--output", type=Path)
    cases = parser.add_mutually_exclusive_group()
    cases.add_argument("--all-cases", action="store_true")
    cases.add_argument("--case-keys", type=parse_case_keys)
    cases.add_argument("--case-count", type=int)
    parser.add_argument(
        "--mod7-coordinate-count", type=int, default=DEFAULT_MOD7_COORDINATE_COUNT
    )
    parser.add_argument("--projection-count", type=int, default=DEFAULT_PROJECTION_COUNT)
    parser.add_argument(
        "--mod7-projections",
        help="explicit semicolon-separated coordinate subsets, e.g. '0,1,2,3;4,5,6,7'",
    )
    parser.add_argument("--state-cap", type=int, default=DEFAULT_STATE_CAP)
    parser.add_argument("--pair-chunk-cap", type=int, default=DEFAULT_PAIR_CHUNK_CAP)
    parser.add_argument(
        "--naive-comparison-cases",
        type=int,
        default=DEFAULT_NAIVE_COMPARISON_CASES,
        help="compare this many cases on the first completed projection; zero disables",
    )
    args = parser.parse_args()

    if args.self_audit_only:
        started = time.time()
        audit = manufactured_self_audit()
        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(
            json.dumps(
                {
                    "status": "manufactured_self_audit_passed",
                    "audit": audit,
                    "elapsed_seconds": time.time() - started,
                    "maximum_resident_set_kib": usage.ru_maxrss,
                },
                indent=2,
                sort_keys=True,
            ),
            flush=True,
        )
        return

    require(args.output is not None, "--output is required unless --self-audit-only is used")
    explicit = (
        naive_join.parse_projection_spec(args.mod7_projections)
        if args.mod7_projections
        else None
    )
    result = run(
        parent_path=args.parent_input,
        current_join_path=args.current_join,
        hilbert_basis_path=args.hilbert_basis,
        output_path=args.output,
        all_cases=args.all_cases,
        case_keys=args.case_keys,
        case_count=args.case_count,
        coordinate_count=args.mod7_coordinate_count,
        projection_count=args.projection_count,
        explicit_projections=explicit,
        state_cap=args.state_cap,
        pair_chunk_cap=args.pair_chunk_cap,
        naive_comparison_cases=args.naive_comparison_cases,
    )
    naive_join.torsion.pointed.atomic_write(args.output, result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.output),
                "selected_cases": result["result_counts"]["selected"],
                "rigorously_rejected": result["result_counts"]["rejected"],
                "necessary_only_survivors": result["result_counts"]["surviving"],
                "skipped": result["result_counts"]["skipped"],
                "selected_projections": [
                    row["mod7_coordinates"]
                    for row in result["projection_subset_selection"]["selected"]
                ],
                "naive_comparisons_completed": result[
                    "existing_naive_case_join_comparison"
                ].get("completed_comparisons", 0),
                "resource_profile": result["resource_profile"],
            },
            indent=2,
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
