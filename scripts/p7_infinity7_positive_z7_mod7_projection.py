#!/usr/bin/env python3
"""Rigorous necessary mod-seven projection sieve for positive p=7, z=7.

The two boundary representatives are loaded from the committed exact orbit
census.  Every one of the 2,160 exact mean leaves is generated from the two
allowed type residues (zero and four), rather than from a hand-written case
list.  The generated leaves are checked against the corrected residue and
catalog-pattern censuses before any sieve decision is made.

For each leaf, every non-enumerated high direction block is annihilated in
the 135-dimensional left-dependency space of the translation-equivariant
281-by-1225 system.  If more than two complete small/medium catalogs remain,
all retained pairs are tested while the other direction blocks are also
annihilated.  A zero-hit projection is a rigorous exclusion because this is
a relaxation of the omitted catalogs.  A passing leaf is only a necessary
mod-seven survivor and is emitted honestly; it is not an edge lift.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_infinity7_positive_z2_mod7_join import (  # noqa: E402
    ContributionFactory,
    DIRECTION_TYPES,
    LABELS,
    canonical_catalog,
    matrix_sha256,
    modular_rank,
    modular_right_nullspace,
    translation_equivariant_system,
)


P = 7
MODULUS = 7
EDGE_COUNT = 29
ORBIT_EVIDENCE = ROOT / "evidence" / "p7_infinity7_positive_zge2_orbits.json"
EXPECTED_ORBIT_CERTIFICATE = "485c006b2320854c91ccff615d2aed013bef2c0683959ce125bf1abcf5e391ba"
EXPECTED_RESIDUES = {"00": 1400, "04": 340, "40": 340, "44": 80}
EXPECTED_LEAVES = 2160
EXPECTED_WEIGHTED_CASES = 60_480
EXPECTED_PATTERN_CENSUS = {
    "00": Counter(
        {
            (0, 0, 7): 8,
            (1, 0, 5): 120,
            (1, 0, 4): 104,
            (1, 0, 3): 32,
            (2, 0, 3): 336,
            (2, 0, 2): 384,
            (2, 0, 1): 192,
            (2, 0, 0): 32,
            (3, 0, 1): 144,
            (3, 0, 0): 48,
        }
    ),
    "single4": Counter(
        {
            (1, 3, 4): 8,
            (1, 3, 3): 32,
            (2, 2, 3): 48,
            (2, 3, 2): 96,
            (2, 3, 1): 192,
            (2, 3, 0): 64,
            (3, 2, 1): 144,
            (3, 2, 0): 48,
            (3, 3, 0): 48,
        }
    ),
    "44": Counter({(2, 6, 0): 32, (3, 5, 0): 48}),
}
ALLOWED_CONDITIONED_DIMENSIONS = {
    0: {135},
    1: {112},
    2: {91, 92},
    3: {73, 74},
    4: {57, 58},
    5: {42},
    6: {28},
    7: {14},
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def atomic_write(path: Path, payload: dict) -> None:
    """Atomically replace path with fsynced, canonical JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


@functools.lru_cache(maxsize=None)
def weak_compositions(total: int, parts: int) -> tuple[tuple[int, ...], ...]:
    if parts == 1:
        return ((total,),)
    return tuple(
        (first, *tail)
        for first in range(total + 1)
        for tail in weak_compositions(total - first, parts - 1)
    )


def boundary_masks(boundary: tuple[int, ...]) -> tuple[int, ...]:
    rows = []
    for labels in LABELS:
        mask = 0
        for point in boundary:
            mask ^= 1 << labels[point]
        rows.append(mask)
    return tuple(rows)


def load_z7_orbits() -> tuple[list[dict], dict]:
    raw_bytes = ORBIT_EVIDENCE.read_bytes()
    source = json.loads(raw_bytes)
    require(source["experiment"] == "p7_infinity7_positive_zge2_orbits", "wrong orbit evidence")
    require(source["status"] == "complete_exact_pair_transversal_orbit_census", "orbit census incomplete")
    require(source["all_required_audits_passed"] is True, "orbit evidence audits did not pass")
    require(source["p"] == P and source["c_H"] == 1, "orbit evidence has wrong scope")
    require(source["census"]["boundary_count_by_z"]["7"] == 56, "z=7 boundary count changed")
    require(source["census"]["orbit_count_by_z"]["7"] == 2, "z=7 orbit count changed")
    require(
        source["orbit_audit"]["orbit_certificate_sha256"] == EXPECTED_ORBIT_CERTIFICATE,
        "orbit certificate changed",
    )
    require(
        tuple(int(row["quadratic_type"]) for row in source["directions"])
        == tuple(DIRECTION_TYPES),
        "direction types disagree with local geometry",
    )

    orbits = []
    for row in source["orbits"]:
        if int(row["z"]) != 7:
            continue
        representative = tuple(int(value) for value in row["representative_finite_field"])
        masks = boundary_masks(representative)
        require(masks == tuple(int(value) for value in row["direction_masks"]), "stored masks changed")
        b_values = tuple(mask.bit_count() for mask in masks)
        undetermined = tuple(index for index, b in enumerate(b_values) if b == P)
        require(len(undetermined) == 7 and sorted(b_values) == [1] + [7] * 7, "not a z=7 line")
        require(tuple(row["undetermined_directions"]) == undetermined, "undetermined directions changed")
        require(int(row["size"]) == 28 and int(row["pair_transversal_multiplicity"]) == 21, "z=7 orbit weight changed")
        orbits.append(
            {
                "source_orbit_index": int(row["orbit_index"]),
                "branch_orbit_index": int(row["branch_orbit_index"]),
                "representative": representative,
                "size": int(row["size"]),
                "masks": masks,
                "b_values": b_values,
                "undetermined": undetermined,
            }
        )
    orbits.sort(key=lambda row: row["branch_orbit_index"])
    require([row["source_orbit_index"] for row in orbits] == [102, 103], "wrong z=7 source orbits")
    require([row["branch_orbit_index"] for row in orbits] == [0, 1], "wrong z=7 branch indices")
    require(sum(row["size"] for row in orbits) == 56, "z=7 orbit weights do not cover branch")
    return orbits, {
        "path": str(ORBIT_EVIDENCE.relative_to(ROOT)),
        "file_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "orbit_certificate_sha256": EXPECTED_ORBIT_CERTIFICATE,
        "source_orbit_indices": [102, 103],
        "orbit_sizes": [28, 28],
        "boundary_count": 56,
        "locally_recomputed_direction_masks": True,
    }


def type_options(orbit: dict, quadratic_type: int, residue: int) -> tuple[dict, ...]:
    directions = tuple(index for index, value in enumerate(DIRECTION_TYPES) if value == quadratic_type)
    require(len(directions) == 4, "direction types are not four by four")
    k = sum(orbit["b_values"][direction] == P for direction in directions)
    total = k if residue == 0 else k - 2
    if total < 0:
        return ()
    rows = []
    for composition in weak_compositions(total, 4):
        q_values = {direction: q for direction, q in zip(directions, composition)}
        means = {}
        for direction in directions:
            b = orbit["b_values"][direction]
            floor = (0 if b == P else 8) if residue == 0 else (4 if b == P else 12)
            means[direction] = floor + 8 * q_values[direction]
        require(sum(means.values()) == 32, "quadratic-type mean budget changed")
        rows.append({"residue": residue, "q_values": q_values, "means": means})
    expected = math.comb(k + 3, 3) if residue == 0 else math.comb(k + 1, 3)
    require(len(rows) == expected, "type composition count changed")
    return tuple(rows)


def exact_mean_leaves(orbits: list[dict]) -> tuple[list[list[dict]], dict]:
    by_orbit: list[list[dict]] = []
    residue_histogram: Counter[str] = Counter()
    pattern_histograms = {key: Counter() for key in EXPECTED_PATTERN_CENSUS}
    per_orbit_residues = []
    for orbit in orbits:
        leaves = []
        local_residues: Counter[str] = Counter()
        for residue_minus, residue_plus in itertools.product((0, 4), repeat=2):
            for minus, plus in itertools.product(
                type_options(orbit, -1, residue_minus),
                type_options(orbit, 1, residue_plus),
            ):
                q_values = [0] * 8
                means = [0] * 8
                residues = [0] * 8
                for option in (minus, plus):
                    for direction, q in option["q_values"].items():
                        q_values[direction] = int(q)
                        means[direction] = int(option["means"][direction])
                        residues[direction] = int(option["residue"])

                classes = []
                levels = []
                for direction in range(8):
                    # This unified level makes every non-enumerated block high
                    # exactly at level >= 2.  Residue-four q=0 is level one.
                    level = q_values[direction] + residues[direction] // 4
                    levels.append(level)
                    if level == 0:
                        catalog_class = "U"
                    elif residues[direction] == 0 and level == 1:
                        catalog_class = "M"
                    elif residues[direction] == 4 and level == 1:
                        catalog_class = "S"
                    else:
                        catalog_class = "H"
                    classes.append(catalog_class)
                require(all((value == "H") == (level >= 2) for value, level in zip(classes, levels)), "high-level classification changed")
                require(Counter(classes).total() == 8, "catalog classes do not cover directions")

                residue_pair = f"{residue_minus}{residue_plus}"
                pattern = (classes.count("H"), classes.count("S"), classes.count("M"))
                family = "00" if residue_pair == "00" else "44" if residue_pair == "44" else "single4"
                residue_histogram[residue_pair] += 1
                local_residues[residue_pair] += 1
                pattern_histograms[family][pattern] += 1
                leaves.append(
                    {
                        "residue_pair": residue_pair,
                        "q_values": tuple(q_values),
                        "scaled_means": tuple(means),
                        "catalog_levels": tuple(levels),
                        "catalog_classes": tuple(classes),
                        "high_directions": tuple(index for index, value in enumerate(classes) if value == "H"),
                        "enumerated_directions": tuple(index for index, value in enumerate(classes) if value in ("S", "M")),
                        "pattern": pattern,
                    }
                )
        require(len(leaves) == 1080, "a z=7 orbit does not have 1,080 leaves")
        by_orbit.append(leaves)
        per_orbit_residues.append(dict(sorted(local_residues.items())))

    require(dict(sorted(residue_histogram.items())) == EXPECTED_RESIDUES, "corrected residue census changed")
    require(sum(len(rows) for rows in by_orbit) == EXPECTED_LEAVES, "2,160-leaf coverage changed")
    for family, expected in EXPECTED_PATTERN_CENSUS.items():
        require(pattern_histograms[family] == expected, f"{family} catalog-pattern census changed")
    require(
        per_orbit_residues
        == [
            {"00": 700, "04": 140, "40": 200, "44": 40},
            {"00": 700, "04": 200, "40": 140, "44": 40},
        ],
        "per-orbit residue census changed",
    )
    return by_orbit, {
        "exact_mean_leaves": EXPECTED_LEAVES,
        "leaves_per_orbit": [1080, 1080],
        "residue_pair_histogram": dict(sorted(residue_histogram.items())),
        "per_orbit_residue_pair_histograms": per_orbit_residues,
        "pattern_histograms": {
            family: {f"H{h}_S{s}_M{m}": count for (h, s, m), count in sorted(rows.items())}
            for family, rows in pattern_histograms.items()
        },
        "weighted_boundary_allocation_cases": EXPECTED_WEIGHTED_CASES,
        "high_definition": "catalog_level=q+residue/4 >= 2; residue-zero q>=2 and residue-four q>=1",
        "residues_two_and_six_absent": True,
    }


def catalog_audit() -> dict:
    expected = {(7, 0): 1, (7, 4): 56, (7, 8): 1764, (1, 8): 1, (1, 12): 56, (1, 16): 1764}
    observed = {(b, mean): len(canonical_catalog(b, mean)) for b, mean in expected}
    require(observed == expected, f"z=7 catalog sizes changed: {observed}")
    return {f"b{b}_mean{mean}": count for (b, mean), count in sorted(observed.items())}


def raw_keys(rows: np.ndarray) -> np.ndarray:
    rows = np.ascontiguousarray(rows, dtype=np.uint8)
    require(rows.ndim == 2 and rows.shape[1] > 0, "signature matrix has bad shape")
    return rows.view(np.dtype((np.void, rows.shape[1]))).reshape(-1)


def count_matches(available: np.ndarray, needed: np.ndarray) -> int:
    keys, counts = np.unique(raw_keys(available), return_counts=True)
    wanted = raw_keys(needed)
    positions = np.searchsorted(keys, wanted)
    valid_indices = np.flatnonzero(positions < len(keys))
    if not len(valid_indices):
        return 0
    equal = keys[positions[valid_indices]] == wanted[valid_indices]
    hits = valid_indices[equal]
    return int(np.sum(counts[positions[hits]], dtype=np.int64))


def exact_projection_join(base: np.ndarray, contributions: tuple[np.ndarray, ...]) -> tuple[int, dict]:
    """Exact zero-syndrome count for zero, one, or two retained catalogs."""
    require(base.ndim == 1 and 0 < len(base) <= 135, "projection base has bad dimension")
    require(len(contributions) <= 2, "projection retained more than two catalogs")
    for contribution in contributions:
        require(contribution.ndim == 2 and contribution.shape[0] == len(base), "catalog projection shape changed")
    sizes = tuple(int(row.shape[1]) for row in contributions)
    if not contributions:
        count = int(not np.any(base % MODULUS))
    elif len(contributions) == 1:
        needed = np.ascontiguousarray(
            (-base[None, :].astype(np.int16)) % MODULUS,
            dtype=np.uint8,
        )
        count = count_matches(contributions[0].T, needed)
    else:
        first, second = contributions
        if first.shape[1] > second.shape[1]:
            first, second = second, first
        needed = np.ascontiguousarray(
            (-base[None, :].astype(np.int16) - second.T.astype(np.int16)) % MODULUS,
            dtype=np.uint8,
        )
        count = count_matches(first.T, needed)
    return count, {
        "retained_catalog_count": len(contributions),
        "catalog_sizes": list(sizes),
        "cartesian_catalog_tuples": math.prod(sizes) if sizes else 1,
        "matching_projected_catalog_tuples": count,
        "signature_encoding": f"{len(base)} raw base-seven bytes",
    }


def join_self_audit() -> dict:
    cases = (
        (np.array([0, 0], dtype=np.uint8), ()),
        (np.array([1, 2], dtype=np.uint8), (np.array([[6, 6, 1], [5, 5, 0]], dtype=np.uint8),)),
        (
            np.array([0, 0], dtype=np.uint8),
            (
                np.array([[0, 0, 1], [0, 0, 2]], dtype=np.uint8),
                np.array([[0, 0, 6], [0, 0, 5]], dtype=np.uint8),
            ),
        ),
    )
    rows = []
    for base, contributions in cases:
        observed, metadata = exact_projection_join(base, contributions)
        brute = 0
        choices = itertools.product(*(range(row.shape[1]) for row in contributions)) if contributions else [()]
        for indices in choices:
            syndrome = base.astype(np.int16).copy()
            for row, index in zip(contributions, indices):
                syndrome += row[:, index]
            brute += int(not np.any(syndrome % MODULUS))
        require(observed == brute and observed > 0, "projection join failed positive brute-force audit")
        rows.append({"brute_force_count": brute, **metadata})
    require(rows[1]["matching_projected_catalog_tuples"] == 2, "one-catalog duplicate audit changed")
    require(rows[2]["matching_projected_catalog_tuples"] > 1, "two-catalog duplicate audit changed")
    return {"passed": True, "positive_hits_and_duplicate_multiplicity_exercised": True, "cases": rows}


class ProjectionFactory:
    def __init__(self, dependencies: np.ndarray):
        self.dependencies = dependencies.astype(np.int64)
        self.contributions = ContributionFactory(self.dependencies)
        self.conditioners: dict[tuple[int, ...], tuple[np.ndarray, dict]] = {}
        self.projected_catalogs: dict[tuple[tuple[int, ...], int, int, int], np.ndarray] = {}

    def conditioner(self, omitted: tuple[int, ...]) -> tuple[np.ndarray, dict]:
        omitted = tuple(sorted(omitted))
        require(len(set(omitted)) == len(omitted) and all(0 <= d < 8 for d in omitted), "bad omitted set")
        if omitted not in self.conditioners:
            if omitted:
                columns = np.concatenate(
                    [np.arange(1 + 35 * direction, 1 + 35 * (direction + 1)) for direction in omitted]
                )
                block = self.dependencies[:, columns]
                coefficients, block_rank = modular_right_nullspace(block.T, MODULUS)
            else:
                columns = np.empty(0, dtype=np.int64)
                coefficients = np.eye(135, dtype=np.int64)
                block_rank = 0
            coefficients = np.ascontiguousarray(coefficients, dtype=np.int64)
            conditioned = coefficients @ self.dependencies % MODULUS
            dimension = len(coefficients)
            require(dimension == 135 - block_rank, "annihilator rank-nullity changed")
            require(dimension in ALLOWED_CONDITIONED_DIMENSIONS[len(omitted)], "conditioned dimension changed")
            require(modular_rank(coefficients, MODULUS) == dimension, "conditioner coefficient rank changed")
            require(modular_rank(conditioned, MODULUS) == dimension, "conditioned dependencies lost rank")
            if omitted:
                require(not np.any(conditioned[:, columns]), "conditioner did not annihilate every omitted block")
            metadata = {
                "omitted_directions": list(omitted),
                "omitted_direction_count": len(omitted),
                "omitted_block_rank": int(block_rank),
                "conditioned_dependency_dimension": dimension,
                "all_omitted_blocks_annihilated": True,
                "coefficient_sha256": matrix_sha256(coefficients.astype(np.uint8)),
            }
            self.conditioners[omitted] = (coefficients, metadata)
        return self.conditioners[omitted]

    def catalog(self, omitted: tuple[int, ...], direction: int, mask: int, mean: int) -> np.ndarray:
        key = (tuple(sorted(omitted)), direction, mask, mean)
        if key not in self.projected_catalogs:
            coefficients, _metadata = self.conditioner(key[0])
            full = self.contributions.get(direction, mask, mean).astype(np.int64)
            self.projected_catalogs[key] = np.ascontiguousarray(coefficients @ full % MODULUS, dtype=np.uint8)
        return self.projected_catalogs[key]

    def audit(self) -> dict:
        histogram: Counter[tuple[int, int]] = Counter()
        for _omitted, (_coefficients, row) in self.conditioners.items():
            histogram[(row["omitted_direction_count"], row["conditioned_dependency_dimension"])] += 1
        return {
            "allowed_dimensions_by_omitted_direction_count": {
                str(count): sorted(values) for count, values in ALLOWED_CONDITIONED_DIMENSIONS.items()
            },
            "observed_conditioner_histogram": {
                f"omit{count}_dim{dimension}": value
                for (count, dimension), value in sorted(histogram.items())
            },
            "conditioner_count": len(self.conditioners),
            "projected_catalog_cache_entries": len(self.projected_catalogs),
            "all_conditioners_rank_and_annihilation_audited": True,
        }


def retained_sets(enumerated: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    if len(enumerated) <= 2:
        return (enumerated,)
    return tuple(itertools.combinations(enumerated, 2))


def relaxed_state(orbit: dict, leaf: dict, dependencies: np.ndarray, factory: ProjectionFactory) -> dict:
    classes = leaf["catalog_classes"]
    high = tuple(leaf["high_directions"])
    enumerated = tuple(leaf["enumerated_directions"])
    require(all(classes[d] == "H" for d in high), "a high block escaped annihilation")
    require(all(classes[d] in ("S", "M") for d in enumerated), "bad enumerated catalog class")

    fixed = (dependencies[:, 0].astype(np.int64) * EDGE_COUNT) % MODULUS
    for direction, catalog_class in enumerate(classes):
        if catalog_class != "U":
            continue
        contribution = factory.contributions.get(
            direction, orbit["masks"][direction], leaf["scaled_means"][direction]
        )
        require(contribution.shape[1] == 1, "fixed direction catalog is not unique")
        fixed = (fixed + contribution[:, 0]) % MODULUS

    plans = []
    for retained in retained_sets(enumerated):
        omitted = tuple(sorted(set(high) | (set(enumerated) - set(retained))))
        require(set(high) <= set(omitted), "a high direction was retained")
        require(set(retained).isdisjoint(omitted), "retained and omitted directions overlap")
        require(set(high) | set(enumerated) == set(omitted) | set(retained), "variable block coverage gap")
        coefficients, conditioning = factory.conditioner(omitted)
        plans.append((
            -int(conditioning["conditioned_dependency_dimension"]),
            retained,
            omitted,
            coefficients,
            conditioning,
        ))
    plans.sort(key=lambda row: (row[0], row[1]))

    projection_rows = []
    for _negative_dimension, retained, omitted, coefficients, conditioning in plans:
        base = np.ascontiguousarray(coefficients @ fixed % MODULUS, dtype=np.uint8)
        contributions = tuple(
            factory.catalog(
                omitted,
                direction,
                orbit["masks"][direction],
                leaf["scaled_means"][direction],
            )
            for direction in retained
        )
        expected_sizes = tuple(56 if classes[direction] == "S" else 1764 for direction in retained)
        require(tuple(row.shape[1] for row in contributions) == expected_sizes, "retained catalog size changed")
        count, join = exact_projection_join(base, contributions)
        row = {
            "retained_directions": list(retained),
            "retained_catalog_classes": [classes[direction] for direction in retained],
            "omitted_directions": list(omitted),
            "conditioned_dependency_dimension": conditioning["conditioned_dependency_dimension"],
            **join,
        }
        projection_rows.append(row)
        if count == 0:
            return {
                "passes_all_projections": False,
                "planned_projection_count": len(plans),
                "tested_projection_count": len(projection_rows),
                "projections": projection_rows,
                "first_failed_projection": row,
            }
    require(len(projection_rows) == len(plans), "survivor did not pass every planned retained pair")
    return {
        "passes_all_projections": True,
        "planned_projection_count": len(plans),
        "tested_projection_count": len(projection_rows),
        "projections": projection_rows,
        "first_failed_projection": None,
    }


def smoke_selection(leaves_by_orbit: list[list[dict]]) -> list[tuple[int, int]]:
    selected = []
    seen = set()
    for orbit_index, leaves in enumerate(leaves_by_orbit):
        for leaf_index, leaf in enumerate(leaves):
            key = (orbit_index, leaf["residue_pair"], leaf["pattern"])
            if key not in seen:
                seen.add(key)
                selected.append((orbit_index, leaf_index))
    arities = {
        min(2, len(leaves_by_orbit[orbit_index][leaf_index]["enumerated_directions"]))
        for orbit_index, leaf_index in selected
    }
    require(arities == {0, 1, 2}, "smoke selection does not exercise all join arities")
    require({orbit for orbit, _leaf in selected} == {0, 1}, "smoke selection misses an orbit")
    return selected


def render_leaf(orbit_index: int, leaf_index: int, orbit: dict, leaf: dict, decision: dict) -> dict:
    return {
        "branch_orbit_index": orbit_index,
        "source_orbit_index": orbit["source_orbit_index"],
        "orbit_leaf_index": leaf_index,
        "orbit_size": orbit["size"],
        "representative_finite_field": list(orbit["representative"]),
        "residue_pair_minus_plus": leaf["residue_pair"],
        "q_values": list(leaf["q_values"]),
        "scaled_means": list(leaf["scaled_means"]),
        "catalog_levels": list(leaf["catalog_levels"]),
        "catalog_classes": list(leaf["catalog_classes"]),
        "catalog_pattern_H_S_M": list(leaf["pattern"]),
        "high_directions_relaxed_by_full_block_annihilation": list(leaf["high_directions"]),
        "enumerated_directions": list(leaf["enumerated_directions"]),
        "necessary_mod7_survivor_only": bool(decision["passes_all_projections"]),
        "planned_projection_count": decision["planned_projection_count"],
        "tested_projection_count": decision["tested_projection_count"],
        "projections": decision["projections"],
    }


def run(smoke_test: bool = False) -> dict:
    started = time.time()
    require(Counter(DIRECTION_TYPES) == Counter({-1: 4, 1: 4}), "direction type census changed")
    matrix, dependencies, linear_audit = translation_equivariant_system()
    require(matrix.shape == (281, 1225), "translation-equivariant matrix shape changed")
    require(linear_audit["rank"] == 146 and dependencies.shape == (135, 281), "rank/dependency census changed")
    require(linear_audit["direction_block_offset"] == 1 and linear_audit["edge_count_rhs"] == 29, "offset/base changed")
    orbits, orbit_source = load_z7_orbits()
    leaves_by_orbit, leaf_audit = exact_mean_leaves(orbits)
    catalogs = catalog_audit()
    join_audit = join_self_audit()
    factory = ProjectionFactory(dependencies)

    selected = smoke_selection(leaves_by_orbit) if smoke_test else [
        (orbit_index, leaf_index)
        for orbit_index, leaves in enumerate(leaves_by_orbit)
        for leaf_index in range(len(leaves))
    ]
    state_cache: dict[tuple, dict] = {}
    survivors = []
    rejection_samples = []
    decision_digest = hashlib.sha256()
    processed_residues: Counter[str] = Counter()
    processed_patterns: Counter[tuple[int, int, int]] = Counter()
    rejected_residues: Counter[str] = Counter()
    logical_projection_histogram: Counter[tuple[int, int, int]] = Counter()
    per_orbit = [Counter() for _ in orbits]

    for orbit_index, leaf_index in selected:
        orbit = orbits[orbit_index]
        leaf = leaves_by_orbit[orbit_index][leaf_index]
        processed_residues[leaf["residue_pair"]] += 1
        processed_patterns[leaf["pattern"]] += 1
        per_orbit[orbit_index]["processed"] += 1
        state_key = (
            orbit_index,
            tuple(leaf["catalog_classes"]),
            tuple(None if value == "H" else mean for value, mean in zip(leaf["catalog_classes"], leaf["scaled_means"])),
        )
        if state_key not in state_cache:
            state_cache[state_key] = relaxed_state(orbit, leaf, dependencies, factory)
        decision = state_cache[state_key]
        for projection in decision["projections"]:
            logical_projection_histogram[
                (
                    int(projection["retained_catalog_count"]),
                    len(projection["omitted_directions"]),
                    int(projection["conditioned_dependency_dimension"]),
                )
            ] += 1
        digest_row = {
            "orbit": orbit_index,
            "leaf": leaf_index,
            "residue": leaf["residue_pair"],
            "q": leaf["q_values"],
            "classes": leaf["catalog_classes"],
            "passing": decision["passes_all_projections"],
            "tested": decision["tested_projection_count"],
            "failure": decision["first_failed_projection"],
        }
        decision_digest.update(json.dumps(digest_row, sort_keys=True, separators=(",", ":")).encode())
        decision_digest.update(b"\n")
        if decision["passes_all_projections"]:
            per_orbit[orbit_index]["surviving"] += 1
            survivors.append(render_leaf(orbit_index, leaf_index, orbit, leaf, decision))
        else:
            per_orbit[orbit_index]["rejected"] += 1
            rejected_residues[leaf["residue_pair"]] += 1
            if len(rejection_samples) < 32:
                sample = render_leaf(orbit_index, leaf_index, orbit, leaf, decision)
                sample["projections"] = [decision["first_failed_projection"]]
                rejection_samples.append(sample)

    full_run = not smoke_test
    processed = len(selected)
    rejected = processed - len(survivors)
    require(processed == rejected + len(survivors), "leaf decision census mismatch")
    if full_run:
        require(processed == EXPECTED_LEAVES, "full run did not process all 2,160 leaves")
        require(dict(sorted(processed_residues.items())) == EXPECTED_RESIDUES, "full processed residue census changed")
        require(sum(row["processed"] for row in per_orbit) == EXPECTED_LEAVES, "per-orbit processing census changed")
        require(all(row["processed"] == 1080 for row in per_orbit), "full run missed an orbit leaf")
    excluded = full_run and not survivors
    weighted_processed = sum(orbits[orbit_index]["size"] for orbit_index, _leaf in selected)
    weighted_survivors = sum(int(row["orbit_size"]) for row in survivors)

    return {
        "experiment": "p7_infinity7_positive_z7_mod7_projection",
        "status": (
            "complete_rigorous_mod7_projection_exclusion"
            if excluded
            else "complete_rigorous_mod7_necessary_sieve_with_survivors"
            if full_run
            else "bounded_smoke_test_only"
        ),
        "p": P,
        "c_H": 1,
        "infinity_in_boundary": True,
        "finite_boundary_points": 7,
        "z": 7,
        "phase": 0,
        "linear_system": linear_audit,
        "orbit_source": orbit_source,
        "mean_leaf_coverage": leaf_audit,
        "catalog_row_counts": catalogs,
        "catalog_source": "complete exact Johnson-slice catalogs imported from the audited z=2 implementation",
        "join_self_audit": join_audit,
        "conditioning_audit": factory.audit(),
        "sieve_method": {
            "high_blocks": "annihilate the full 35-column direction block (a superset relaxation of the exact high catalog)",
            "more_than_two_enumerated_catalogs": "test every retained pair and annihilate every other variable direction block",
            "at_most_two_enumerated_catalogs": "retain all and exact-join complete catalogs",
            "rejection_semantics": "one zero-hit necessary projection rigorously excludes the leaf",
            "survivor_semantics": "passed every planned projection; necessary mod-seven survivor only, not an edge lift",
        },
        "smoke_test": smoke_test,
        "smoke_limitations": (
            "Processes one representative of each orbit/residue/catalog-pattern class; it validates plumbing and assertions but cannot exclude z=7."
            if smoke_test
            else None
        ),
        "full_run": full_run,
        "full_exact_mean_leaves": EXPECTED_LEAVES,
        "processed_exact_mean_leaves": processed,
        "rejected_exact_mean_leaves": rejected,
        "surviving_exact_mean_leaves": len(survivors),
        "processed_residue_pair_histogram": dict(sorted(processed_residues.items())),
        "rejected_residue_pair_histogram": dict(sorted(rejected_residues.items())),
        "processed_pattern_histogram": {
            f"H{h}_S{s}_M{m}": count for (h, s, m), count in sorted(processed_patterns.items())
        },
        "processed_weighted_boundary_allocation_cases": weighted_processed,
        "surviving_weighted_boundary_allocation_cases": weighted_survivors,
        "full_weighted_boundary_allocation_cases": EXPECTED_WEIGHTED_CASES,
        "computed_relaxed_states": len(state_cache),
        "logical_projection_test_histogram": {
            f"retain{retain}_omit{omit}_dim{dimension}": count
            for (retain, omit, dimension), count in sorted(logical_projection_histogram.items())
        },
        "all_case_decisions_sha256": decision_digest.hexdigest(),
        "modular_passing_is_edge_feasibility": False,
        "exact_edge_lift_claimed": False,
        "z7_branch_excluded": excluded,
        "per_orbit_summary": [
            {
                "branch_orbit_index": index,
                "source_orbit_index": orbits[index]["source_orbit_index"],
                "orbit_size": orbits[index]["size"],
                "processed_leaves": row["processed"],
                "rejected_leaves": row["rejected"],
                "surviving_leaves": row["surviving"],
            }
            for index, row in enumerate(per_orbit)
        ],
        "rejection_samples": rejection_samples,
        "survivor_cases": survivors,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--smoke-test",
        action="store_true",
        help="run a bounded pattern-representative subset; never claims z=7 closure",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    result = run(smoke_test=args.smoke_test)
    atomic_write(args.output, result)
    if not args.quiet:
        summary = {
            key: value
            for key, value in result.items()
            if key not in {"per_orbit_summary", "rejection_samples", "survivor_cases"}
        }
        print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
