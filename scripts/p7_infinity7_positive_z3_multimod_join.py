#!/usr/bin/env python3
"""Extract and test the exact positive p=7 infinity+7 z=3 mod-7 tuples.

The audited mod-seven run leaves exactly two labeled orbit/mean cases.  Each
case has three complete 1,764-row Johnson-slice catalogs and exactly four
catalog tuples satisfying all 135 dependencies modulo seven.  This script
reconstructs those leaves from :mod:`p7_infinity7_positive_z3_mod7_join`,
extracts every matching tuple with its three catalog-row indices, rebuilds
the corresponding *integer* 281-entry right side, and tests that same right
side against complete left-nullspace bases modulo 3, 5, 7, and 11.

The pair side of each 1+2 extraction is streamed in chunks of at most 100,000
syndrome rows.  A failure at any additional prime is a rigorous exclusion of
that mod-seven tuple.  The z=3 branch is declared closed if and only if all
eight extracted mod-seven tuples fail modulo at least one of 3, 5, and 11.
Passing every tested prime remains only a necessary condition for an edge
lift.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z3_mod7_join as audited  # noqa: E402


DEFAULT_INPUT = Path("/tmp/p7_inf7_z3_full.json")
MODULI = (3, 5, 7, 11)
ADDITIONAL_MODULI = (3, 5, 11)
EXPECTED_RANKS = {3: 161, 5: 167, 7: 146, 11: 167}
EXPECTED_DEPENDENCY_DIMENSIONS = {3: 120, 5: 114, 7: 135, 11: 114}
EXPECTED_SURVIVOR_KEYS = ((0, 30), (8, 1))
EXPECTED_MATCHES_PER_CASE = 4
EXPECTED_TOTAL_MOD7_TUPLES = 8
EXPECTED_INPUT_DECISION_SHA256 = (
    "c39fb7f530a6380c09d0bf300d6d249df2304370867066e6ce74de62906f275f"
)
MAX_PAIR_CHUNK_ROWS = 100_000
NEGATIVE_MOD7 = np.asarray((0, 6, 5, 4, 3, 2, 1), dtype=np.uint8)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def uint8_sha256(values: np.ndarray) -> str:
    canonical = np.ascontiguousarray(values, dtype=np.uint8)
    return sha256_bytes(canonical.tobytes())


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256_bytes(encoded)


def load_and_audit_input(path: Path) -> tuple[dict, list[dict], dict]:
    encoded = path.read_bytes()
    payload = json.loads(encoded)
    expected_scalars = {
        "experiment": "p7_infinity7_positive_z3_mod7_join",
        "status": "complete_rigorous_mod7_necessary_sieve_with_survivors",
        "p": 7,
        "c_H": 1,
        "infinity_in_boundary": True,
        "finite_boundary_points": 7,
        "z": 3,
        "phase": 0,
        "full_run": True,
        "processed_orbits": 10,
        "full_orbit_count": 10,
        "processed_exact_mean_leaves": 400,
        "processed_weighted_boundary_allocation_cases": 225_792,
        "rejected_weighted_boundary_allocation_cases": 225_008,
        "surviving_cases": 2,
        "surviving_weighted_boundary_allocation_cases": 784,
        "weighted_exact_mod7_catalog_tuples": 3_136,
        "z3_branch_excluded": False,
    }
    for key, expected in expected_scalars.items():
        require(payload.get(key) == expected, f"input field {key!r} changed")
    require(
        payload.get("all_case_decisions_sha256")
        == EXPECTED_INPUT_DECISION_SHA256,
        "input decision certificate changed",
    )
    require(
        payload.get("modular_passing_is_edge_feasibility") is False,
        "input overstates modular feasibility",
    )
    require(
        payload["mean_leaf_coverage"].get("all_400_corrected_leaves_covered")
        is True,
        "input does not certify all corrected mean leaves",
    )
    require(
        payload["mean_leaf_coverage"].get("exact_mean_leaves") == 400,
        "input corrected mean-leaf census changed",
    )
    require(
        sum(int(value) for value in payload["rejected_kind_histogram"].values())
        == 398,
        "input does not reject exactly the other 398 mean leaves",
    )

    survivors = sorted(
        payload.get("survivor_cases", []),
        key=lambda row: (int(row["orbit_index"]), int(row["orbit_leaf_index"])),
    )
    keys = tuple(
        (int(row["orbit_index"]), int(row["orbit_leaf_index"]))
        for row in survivors
    )
    require(keys == EXPECTED_SURVIVOR_KEYS, f"unexpected input survivors: {keys}")
    for row in survivors:
        require(row.get("leaf_kind") == "residue00_three_catalog",
                "input survivor is not an exact three-catalog leaf")
        require(row.get("catalog_pattern") == "M1764^3",
                "input survivor catalog pattern changed")
        require(row.get("high_direction") is None,
                "input survivor unexpectedly contains a high direction")
        require(row.get("high_direction_block_relaxed") is False,
                "input survivor unexpectedly used a high-block relaxation")
        require(
            row.get("mod7_test")
            == "exact_full_135_three_catalog_1plus2_join",
            "input survivor did not use the exact full mod-seven join",
        )
        require(row.get("matching_remaining_catalog_tuples") == 4,
                "input survivor does not report four matching tuples")
        require(row.get("exact_mod7_catalog_tuples") == 4,
                "input exact mod-seven tuple count changed")
        require(row["join"].get("catalog_sizes") == [1764, 1764, 1764],
                "input survivor is not a 1764^3 join")
        require(row["join"].get("cartesian_catalog_tuples") == 1764**3,
                "input Cartesian catalog count changed")

    return payload, survivors, {
        "path": str(path),
        "sha256": sha256_bytes(encoded),
        "decision_sha256": payload["all_case_decisions_sha256"],
        "full_run": True,
        "processed_exact_mean_leaves": payload["processed_exact_mean_leaves"],
        "processed_weighted_boundary_allocation_cases": payload[
            "processed_weighted_boundary_allocation_cases"
        ],
        "survivor_keys": [list(key) for key in keys],
        "reported_mod7_tuples_per_case": EXPECTED_MATCHES_PER_CASE,
        "reported_total_mod7_tuples": EXPECTED_TOTAL_MOD7_TUPLES,
    }


def modular_systems(
    matrix: np.ndarray,
    input_payload: dict,
) -> tuple[dict[int, np.ndarray], dict]:
    require(matrix.shape == (281, 1225), "integer system shape changed")
    full = audited.equation_matrix()
    require(full.shape == (282, 1225), "source equation matrix shape changed")
    rebuilt = np.concatenate((full[:1], full[2:]), axis=0)
    require(np.array_equal(matrix, rebuilt),
            "281-row system is not exactly source rows 0,2,...,281")
    require(np.all(matrix[0] == 1), "edge-count row changed")

    matrix_hash = audited.matrix_sha256(matrix)
    require(
        input_payload["linear_system"].get("matrix_sha256") == matrix_hash,
        "input and reconstructed integer matrices differ",
    )
    calibration = (17 * np.arange(matrix.shape[1], dtype=np.int64) + 3) % 2
    manufactured_rhs = matrix.astype(np.int64) @ calibration

    bases: dict[int, np.ndarray] = {}
    rows = []
    for modulus in MODULI:
        rank, dependencies = audited.left_dependencies(matrix, modulus)
        dependencies = np.ascontiguousarray(dependencies, dtype=np.int64)
        expected_dimension = matrix.shape[0] - rank
        require(rank == EXPECTED_RANKS[modulus],
                f"rank changed modulo {modulus}: {rank}")
        require(
            expected_dimension == EXPECTED_DEPENDENCY_DIMENSIONS[modulus],
            f"dependency dimension changed modulo {modulus}",
        )
        require(
            dependencies.shape == (expected_dimension, matrix.shape[0]),
            f"dependency basis shape changed modulo {modulus}",
        )
        require(audited.modular_rank(matrix, modulus) == rank,
                f"independent matrix rank audit failed modulo {modulus}")
        require(
            audited.modular_rank(dependencies, modulus) == expected_dimension,
            f"dependency basis is not independent modulo {modulus}",
        )
        require(
            not np.any(
                dependencies @ (matrix.astype(np.int64) % modulus) % modulus
            ),
            f"left-nullspace audit failed modulo {modulus}",
        )
        require(
            not np.any(dependencies @ (manufactured_rhs % modulus) % modulus),
            f"manufactured consistent RHS was rejected modulo {modulus}",
        )
        basis_hash = uint8_sha256(dependencies)
        if modulus == 7:
            require(
                input_payload["linear_system"].get("rank") == rank
                and input_payload["linear_system"].get(
                    "left_dependency_dimension"
                )
                == expected_dimension,
                "input mod-seven rank metadata changed",
            )
            require(
                input_payload["linear_system"].get("dependency_sha256")
                == basis_hash,
                "input and reconstructed mod-seven dependency bases differ",
            )
        bases[modulus] = dependencies
        rows.append(
            {
                "modulus": modulus,
                "rank": rank,
                "dependency_dimension": expected_dimension,
                "rank_plus_dependency_dimension": rank + expected_dimension,
                "basis_rank": audited.modular_rank(dependencies, modulus),
                "left_null_audit": True,
                "manufactured_rhs_audit": True,
                "dependency_basis_sha256_uint8": basis_hash,
            }
        )

    return bases, {
        "construction": (
            "the 282-row integer source matrix with only distinguished-edge "
            "row 1 deleted"
        ),
        "source_shape": list(full.shape),
        "shape": list(matrix.shape),
        "edge_variables": matrix.shape[1],
        "edge_count_rhs": audited.EDGE_COUNT,
        "direction_block_offset": 1,
        "direction_block_width": 35,
        "integer_matrix_dtype": str(matrix.dtype),
        "integer_matrix_sha256": matrix_hash,
        "source_integer_matrix_sha256": audited.matrix_sha256(full),
        "same_integer_matrix_used_for_every_modulus": True,
        "complete_left_dependency_bases": rows,
    }


def reconstruct_survivor(
    input_row: dict,
    orbits: list[dict],
    leaves_by_orbit: list[list[dict]],
) -> tuple[int, int, dict, dict]:
    orbit_index = int(input_row["orbit_index"])
    leaf_index = int(input_row["orbit_leaf_index"])
    require(0 <= orbit_index < len(orbits), "survivor orbit index is outside evidence")
    require(0 <= leaf_index < len(leaves_by_orbit[orbit_index]),
            "survivor leaf index is outside reconstructed leaves")
    orbit = orbits[orbit_index]
    leaf = leaves_by_orbit[orbit_index][leaf_index]

    comparisons = {
        "source_orbit_index": int(orbit["source_orbit_index"]),
        "orbit_size": int(orbit["size"]),
        "representative_finite_field": list(orbit["representative"]),
        "direction_types": list(audited.DIRECTION_TYPES),
        "b_values": list(orbit["b_values"]),
        "undetermined_directions": list(orbit["undetermined"]),
        "undetermined_type_counts": {
            str(key): value
            for key, value in sorted(orbit["undetermined_type_counts"].items())
        },
        "residue": leaf["residue"],
        "leaf_kind": leaf["kind"],
        "catalog_pattern": leaf["catalog_pattern"],
        "q_values": None
        if leaf["q_values"] is None
        else list(leaf["q_values"]),
        "residue_increments": None
        if leaf["residue_increments"] is None
        else list(leaf["residue_increments"]),
        "scaled_means": list(leaf["means"]),
        "catalog_roles": list(leaf["roles"]),
        "variable_directions": list(leaf["variable_directions"]),
        "high_direction": leaf["high_direction"],
    }
    for key, rebuilt in comparisons.items():
        require(input_row.get(key) == rebuilt,
                f"input survivor {orbit_index}/{leaf_index} field {key} changed")
    require(leaf["kind"] == "residue00_three_catalog",
            "reconstructed survivor is not a three-catalog leaf")
    require(leaf["high_direction"] is None,
            "reconstructed survivor has a relaxed high direction")
    require(len(leaf["variable_directions"]) == 3,
            "reconstructed survivor does not have three variable directions")
    return orbit_index, leaf_index, orbit, leaf


def build_catalogs(orbit: dict, leaf: dict) -> tuple[dict[int, np.ndarray], list[dict]]:
    catalogs: dict[int, np.ndarray] = {}
    rows = []
    for direction, (mask, b, mean, role) in enumerate(
        zip(orbit["masks"], orbit["b_values"], leaf["means"], leaf["roles"])
    ):
        require(role in ("fixed", "M"),
                "survivor contains a non-exact catalog role")
        catalog = np.ascontiguousarray(
            audited.mapped_catalog(int(mask), int(mean)), dtype=np.int16
        )
        expected_size = 1 if role == "fixed" else audited.catalog_size("M", int(b))
        require(len(catalog) == expected_size,
                f"direction {direction} catalog size changed")
        if role == "M":
            require(expected_size == 1764,
                    f"direction {direction} is not a 1,764-row catalog")
        require(catalog.shape == (expected_size, 35),
                f"direction {direction} catalog shape changed")
        require(np.all((0 <= catalog) & (catalog <= 13)),
                f"direction {direction} catalog left the integer score range")
        catalogs[direction] = catalog
        rows.append(
            {
                "direction": direction,
                "odd_fibre_mask": int(mask),
                "b": int(b),
                "scaled_mean": int(mean),
                "role": role,
                "rows": len(catalog),
                "catalog_sha256_uint8": uint8_sha256(catalog),
            }
        )
    return catalogs, rows


def mod7_join_data(
    dependencies: np.ndarray,
    orbit: dict,
    leaf: dict,
    catalogs: dict[int, np.ndarray],
) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
    base = dependencies[:, 0].astype(np.int64) * audited.EDGE_COUNT
    contributions = []
    variable_directions = tuple(
        int(value) for value in leaf["variable_directions"]
    )
    variable = set(variable_directions)
    require(
        tuple(direction for direction in range(audited.P + 1) if direction in variable)
        == variable_directions,
        "variable-direction order is not canonical",
    )
    factory = audited.ContributionFactory(dependencies)
    for direction in range(audited.P + 1):
        catalog = catalogs[direction]
        values = catalog.astype(np.int64)
        block = dependencies[
            :, 1 + 35 * direction : 1 + 35 * (direction + 1)
        ]
        direct = np.ascontiguousarray(block @ (13 - values).T % 7, dtype=np.uint8)
        imported = factory.get(
            direction, int(orbit["masks"][direction]), int(leaf["means"][direction])
        )
        require(np.array_equal(direct, imported),
                f"direction {direction} contribution reconstruction failed")
        if direction in variable:
            contributions.append(direct)
        else:
            require(direct.shape[1] == 1,
                    f"fixed direction {direction} is not unique")
            base += direct[:, 0].astype(np.int64)
    require(len(contributions) == 3, "expected exactly three variable contributions")
    require(all(matrix.shape == (135, 1764) for matrix in contributions),
            "mod-seven contribution dimensions changed")
    return (
        np.ascontiguousarray(base % 7, dtype=np.uint8),
        tuple(contributions),
    )


def extract_three_catalog_matches(
    base: np.ndarray,
    contributions: tuple[np.ndarray, ...],
    max_pair_chunk_rows: int = MAX_PAIR_CHUNK_ROWS,
) -> tuple[list[tuple[int, int, int]], dict]:
    """Return all row-index triples solving base+c0+c1+c2=0 modulo seven."""
    require(len(contributions) == 3, "extractor needs exactly three catalogs")
    dimension = len(base)
    require(base.shape == (dimension,) and dimension > 0, "join base has bad shape")
    require(max_pair_chunk_rows > 0, "pair chunk bound must be positive")
    require(all(matrix.shape[0] == dimension for matrix in contributions),
            "catalog syndrome dimensions disagree")
    sizes = tuple(int(matrix.shape[1]) for matrix in contributions)
    require(all(size > 0 for size in sizes), "empty catalog passed to extractor")
    require(sizes[2] <= max_pair_chunk_rows,
            "third catalog alone exceeds the requested pair chunk bound")

    single_rows = np.ascontiguousarray(contributions[0].T, dtype=np.uint8)
    single_keys = audited.row_keys(single_rows)
    order = np.argsort(single_keys, kind="stable")
    sorted_keys = single_keys[order]
    require(len(sorted_keys) == sizes[0], "single-side signature count changed")

    third_rows = np.ascontiguousarray(contributions[2].T, dtype=np.uint8)
    second_block_width = max(1, max_pair_chunk_rows // sizes[2])
    matches: list[tuple[int, int, int]] = []
    streamed_pairs = 0
    peak_rows = 0
    chunks = 0
    for start in range(0, sizes[1], second_block_width):
        stop = min(start + second_block_width, sizes[1])
        second_rows = np.ascontiguousarray(
            contributions[1][:, start:stop].T, dtype=np.uint8
        )
        needed = second_rows[:, None, :] + third_rows[None, :, :]
        np.remainder(needed, 7, out=needed)
        needed = np.ascontiguousarray(needed.reshape(-1, dimension), dtype=np.uint8)
        require(len(needed) <= max_pair_chunk_rows,
                "pair stream exceeded its row bound")
        np.add(needed, base[None, :], out=needed)
        np.remainder(needed, 7, out=needed)
        needed[:] = NEGATIVE_MOD7[needed]

        needed_keys = audited.row_keys(needed)
        lower = np.searchsorted(sorted_keys, needed_keys, side="left")
        upper = np.searchsorted(sorted_keys, needed_keys, side="right")
        hit_rows = np.flatnonzero(lower < upper)
        for flat_index in hit_rows.tolist():
            second_index = start + flat_index // sizes[2]
            third_index = flat_index % sizes[2]
            for sorted_position in range(
                int(lower[flat_index]), int(upper[flat_index])
            ):
                first_index = int(order[sorted_position])
                match = (first_index, second_index, third_index)
                syndrome = base.astype(np.int16).copy()
                for matrix, row_index in zip(contributions, match):
                    syndrome += matrix[:, row_index].astype(np.int16)
                require(not np.any(syndrome % 7),
                        "extracted row-index tuple is not a mod-seven match")
                matches.append(match)
        streamed_pairs += len(needed)
        peak_rows = max(peak_rows, len(needed))
        chunks += 1

    matches.sort()
    require(streamed_pairs == sizes[1] * sizes[2],
            "pair extractor missed streamed catalog pairs")
    require(len(set(matches)) == len(matches),
            "pair extractor emitted duplicate row-index tuples")
    return matches, {
        "partition": [[0], [1, 2]],
        "catalog_sizes": list(sizes),
        "materialized_single_rows": sizes[0],
        "streamed_pair_rows": streamed_pairs,
        "streamed_pair_chunks": chunks,
        "streamed_pair_peak_rows": peak_rows,
        "configured_pair_chunk_bound_rows": max_pair_chunk_rows,
        "raw_base_seven_signature_bytes": dimension,
        "duplicate_single_signatures_preserve_all_row_indices": True,
        "matching_row_index_tuples": len(matches),
        "row_index_tuple_sha256": json_sha256(matches),
    }


def extractor_self_audit() -> dict:
    base = np.asarray([0, 0], dtype=np.uint8)
    contributions = (
        np.asarray([[0, 0, 1], [0, 0, 0]], dtype=np.uint8),
        np.asarray([[0, 1, 6], [0, 0, 0]], dtype=np.uint8),
        np.asarray([[0, 1, 6], [0, 0, 0]], dtype=np.uint8),
    )
    extracted, metadata = extract_three_catalog_matches(
        base, contributions, max_pair_chunk_rows=4
    )
    brute = []
    for indices in itertools.product(*(range(matrix.shape[1]) for matrix in contributions)):
        syndrome = base.astype(np.int16).copy()
        for matrix, row_index in zip(contributions, indices):
            syndrome += matrix[:, row_index].astype(np.int16)
        if not np.any(syndrome % 7):
            brute.append(tuple(int(value) for value in indices))
    require(extracted == brute, "row-index extractor failed direct brute-force audit")
    require(len(extracted) > len(set(index[0] for index in extracted)),
            "extractor self-audit did not exercise duplicate signatures")
    require(metadata["streamed_pair_chunks"] > 1,
            "extractor self-audit did not exercise pair chunking")
    return {
        "passed": True,
        "direct_brute_force_match": True,
        "duplicate_signature_multiplicity_exercised": True,
        "multi_chunk_streaming_exercised": True,
        "matching_tuples": len(extracted),
        "matching_tuple_sha256": json_sha256(extracted),
    }


def integer_rhs_for_tuple(
    catalogs: dict[int, np.ndarray],
    variable_directions: tuple[int, ...],
    row_indices: tuple[int, int, int],
) -> np.ndarray:
    selected = dict(zip(variable_directions, row_indices))
    rhs = np.empty(281, dtype=np.int16)
    rhs[0] = audited.EDGE_COUNT
    for direction in range(audited.P + 1):
        catalog = catalogs[direction]
        row_index = selected.get(direction, 0)
        if direction not in selected:
            require(len(catalog) == 1,
                    f"non-variable direction {direction} is not unique")
        require(0 <= row_index < len(catalog),
                f"catalog row index outside direction {direction}")
        block = slice(1 + 35 * direction, 1 + 35 * (direction + 1))
        rhs[block] = 13 - catalog[row_index]
    require(np.all((0 <= rhs[1:]) & (rhs[1:] <= 13)),
            "integer bad-count right side left the range 0..13")
    return rhs


def test_integer_tuple(
    rhs: np.ndarray,
    bases: dict[int, np.ndarray],
) -> tuple[dict[str, dict], list[int]]:
    tests: dict[str, dict] = {}
    eliminating = []
    for modulus in MODULI:
        dependencies = bases[modulus]
        syndrome = np.ascontiguousarray(
            dependencies @ (rhs.astype(np.int64) % modulus) % modulus,
            dtype=np.uint8,
        )
        nonzero = np.flatnonzero(syndrome)
        passed = not len(nonzero)
        if modulus == 7:
            require(passed, "extracted mod-seven tuple failed the direct integer audit")
        elif not passed:
            eliminating.append(modulus)
        tests[str(modulus)] = {
            "dependency_dimension": len(dependencies),
            "passes_all_dependencies": passed,
            "nonzero_dependency_coordinates": int(len(nonzero)),
            "first_nonzero_dependency_coordinate": (
                None if not len(nonzero) else int(nonzero[0])
            ),
            "syndrome_sha256_uint8": uint8_sha256(syndrome),
            "syndrome": syndrome.tolist(),
        }
    return tests, eliminating


def run(input_path: Path) -> dict:
    started = time.time()
    input_payload, input_survivors, input_audit = load_and_audit_input(input_path)
    matrix, imported_mod7_dependencies, imported_linear = (
        audited.translation_equivariant_system()
    )
    bases, linear_audit = modular_systems(matrix, input_payload)
    require(np.array_equal(imported_mod7_dependencies, bases[7]),
            "audited module returned a different mod-seven basis")
    require(imported_linear["matrix_sha256"] == linear_audit["integer_matrix_sha256"],
            "audited module matrix metadata changed")

    orbits, orbit_audit = audited.evidence_z3_orbits()
    leaves_by_orbit, leaf_audit = audited.exact_mean_leaves(orbits)
    require(len(orbits) == 10 and sum(map(len, leaves_by_orbit)) == 400,
            "audited z=3 reconstruction coverage changed")
    extraction_audit = extractor_self_audit()

    case_rows = []
    tuple_rows = []
    for input_row in input_survivors:
        orbit_index, leaf_index, orbit, leaf = reconstruct_survivor(
            input_row, orbits, leaves_by_orbit
        )
        catalogs, catalog_audit = build_catalogs(orbit, leaf)
        base, contributions = mod7_join_data(bases[7], orbit, leaf, catalogs)
        matches, join_audit = extract_three_catalog_matches(base, contributions)
        require(len(matches) == EXPECTED_MATCHES_PER_CASE,
                f"case {orbit_index}/{leaf_index} did not reconstruct four tuples")
        require(len(matches) == int(input_row["exact_mod7_catalog_tuples"]),
                f"case {orbit_index}/{leaf_index} disagrees with input tuple count")

        variable_directions = tuple(int(value) for value in leaf["variable_directions"])
        case_tuple_rows = []
        for case_tuple_index, row_indices in enumerate(matches):
            rhs = integer_rhs_for_tuple(catalogs, variable_directions, row_indices)
            tests, eliminating = test_integer_tuple(rhs, bases)
            selected_hashes = {
                str(direction): uint8_sha256(catalogs[direction][row_index])
                for direction, row_index in zip(variable_directions, row_indices)
            }
            record = {
                "orbit_index": orbit_index,
                "orbit_leaf_index": leaf_index,
                "case_tuple_index": case_tuple_index,
                "variable_directions": list(variable_directions),
                "catalog_row_indices_variable_order": list(row_indices),
                "catalog_row_indices_by_direction": {
                    str(direction): row_index
                    for direction, row_index in zip(variable_directions, row_indices)
                },
                "selected_catalog_row_sha256_uint8_by_direction": selected_hashes,
                "integer_rhs_encoding": "281 uint8 values: edge count, then 8x35 bad counts",
                "integer_rhs_sha256_uint8": uint8_sha256(rhs),
                "integer_rhs": rhs.tolist(),
                "modular_dependency_tests": tests,
                "eliminating_additional_primes": eliminating,
                "fails_at_least_one_additional_prime": bool(eliminating),
                "survives_all_additional_primes": not bool(eliminating),
            }
            case_tuple_rows.append(record)
            tuple_rows.append(record)

        case_rows.append(
            {
                "orbit_index": orbit_index,
                "source_orbit_index": int(orbit["source_orbit_index"]),
                "orbit_leaf_index": leaf_index,
                "orbit_size": int(orbit["size"]),
                "representative_finite_field": list(orbit["representative"]),
                "scaled_means": list(leaf["means"]),
                "variable_directions": list(variable_directions),
                "catalog_audit": catalog_audit,
                "mod7_extraction": join_audit,
                "input_reported_mod7_tuples": int(
                    input_row["exact_mod7_catalog_tuples"]
                ),
                "extracted_mod7_tuples": len(matches),
                "tuple_decisions": case_tuple_rows,
            }
        )

    require(len(tuple_rows) == EXPECTED_TOTAL_MOD7_TUPLES,
            "total extracted mod-seven tuple count changed")
    row_index_certificate = [
        [
            int(row["orbit_index"]),
            int(row["orbit_leaf_index"]),
            *[int(value) for value in row["catalog_row_indices_variable_order"]],
        ]
        for row in tuple_rows
    ]
    decision_certificate = [
        {
            "orbit": row["orbit_index"],
            "leaf": row["orbit_leaf_index"],
            "indices": row["catalog_row_indices_variable_order"],
            "passes": {
                modulus: row["modular_dependency_tests"][str(modulus)][
                    "passes_all_dependencies"
                ]
                for modulus in MODULI
            },
        }
        for row in tuple_rows
    ]
    additional_survivors = [
        row for row in tuple_rows if row["survives_all_additional_primes"]
    ]
    every_mod7_tuple_eliminated = all(
        row["fails_at_least_one_additional_prime"] for row in tuple_rows
    )
    z3_excluded = every_mod7_tuple_eliminated
    require(
        z3_excluded
        == (
            len(tuple_rows) == EXPECTED_TOTAL_MOD7_TUPLES
            and all(
                not all(
                    row["modular_dependency_tests"][str(modulus)][
                        "passes_all_dependencies"
                    ]
                    for modulus in ADDITIONAL_MODULI
                )
                for row in tuple_rows
            )
        ),
        "z=3 exclusion is not exactly the required all-tuples condition",
    )

    audited_source = Path(audited.__file__).resolve()
    this_source = Path(__file__).resolve()
    return {
        "experiment": "p7_infinity7_positive_z3_multimod_join",
        "status": (
            "complete_rigorous_multimod_exclusion"
            if z3_excluded
            else "complete_rigorous_multimod_necessary_sieve_with_survivors"
        ),
        "p": 7,
        "c_H": 1,
        "infinity_in_boundary": True,
        "finite_boundary_points": 7,
        "z": 3,
        "phase": 0,
        "input_audit": input_audit,
        "audited_z3_module": {
            "path": str(audited_source.relative_to(ROOT)),
            "sha256": sha256_bytes(audited_source.read_bytes()),
        },
        "this_script": {
            "path": str(this_source.relative_to(ROOT)),
            "sha256": sha256_bytes(this_source.read_bytes()),
        },
        "orbit_reconstruction": orbit_audit,
        "mean_leaf_reconstruction": leaf_audit,
        "linear_system": linear_audit,
        "extractor_self_audit": extraction_audit,
        "tested_moduli": list(MODULI),
        "additional_elimination_moduli": list(ADDITIONAL_MODULI),
        "same_integer_rhs_tested_at_every_modulus": True,
        "input_surviving_mean_cases": len(input_survivors),
        "extracted_mod7_catalog_tuples": len(tuple_rows),
        "expected_mod7_catalog_tuples": EXPECTED_TOTAL_MOD7_TUPLES,
        "mod7_row_index_certificate_sha256": json_sha256(row_index_certificate),
        "multimod_decision_certificate_sha256": json_sha256(decision_certificate),
        "additional_prime_surviving_tuples": len(additional_survivors),
        "all_mod7_tuples_fail_at_least_one_additional_prime": (
            every_mod7_tuple_eliminated
        ),
        "z3_branch_excluded": z3_excluded,
        "exclusion_iff_all_mod7_tuples_fail_an_additional_prime": True,
        "modular_passing_is_edge_feasibility": False,
        "survivor_cases": additional_survivors,
        "cases": case_rows,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    out = run(args.input)
    audited.atomic_write(args.output, out)
    if not args.quiet:
        summary = {
            key: value
            for key, value in out.items()
            if key not in {"cases", "survivor_cases"}
        }
        print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
