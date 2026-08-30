#!/usr/bin/env python3
"""Exact orbit census for positive p=7 infinity-plus-seven with z >= 2.

For a seven-point affine boundary, ``b_d=7`` means that every fibre in
direction ``d`` contains exactly one boundary point.  Thus a boundary with
at least two such directions is a permutation graph for each pair of its
undetermined directions.  This script enumerates the exact incidence set

    C(8, 2) * 7! = 141,120

by direction-pair transversals; it never scans ``C(49, 7)`` subsets.  A
boundary with ``z`` undetermined directions occurs exactly ``C(z, 2)``
times.  The resulting distinct boundaries are decomposed under the full
2,352-element square affine-semilinear group of GF(49).

All census, multiplicity, group-action, mask, and orbit-size identities are
hard assertions.  The JSON output contains one canonical representative per
orbit and enough direction data to reconstruct every reported branch.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)


P = 7
Q = P * P
FINITE_BOUNDARY_SIZE = 7
EXPECTED_BOUNDARY_COUNTS = {2: 123_480, 3: 5_488, 7: 56}
EXPECTED_ORBIT_SIZE_HISTOGRAMS = {
    2: {588: 18, 1_176: 52, 2_352: 22},
    3: {392: 6, 784: 4},
    7: {28: 2},
}
EXPECTED_ORBIT_COUNTS = {2: 92, 3: 10, 7: 2}
EXPECTED_Z2_ORBIT_TYPE_SPLIT = {"same_type": 48, "split_type": 44}
EXPECTED_GROUP_SIZE = 2_352
EXPECTED_PAIR_INCIDENCES = math.comb(P + 1, 2) * math.factorial(P)


def atomic_write(path: Path, payload: dict[str, object]) -> None:
    """Atomically replace ``path`` with canonical, indented JSON."""
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


def boundary_digest(boundaries: set[tuple[int, ...]]) -> str:
    """Digest sorted seven-byte GF(49) representatives unambiguously."""
    digest = hashlib.sha256()
    for boundary in sorted(boundaries):
        if len(boundary) != FINITE_BOUNDARY_SIZE or any(not 0 <= u < Q for u in boundary):
            raise AssertionError("invalid boundary passed to digest")
        digest.update(bytes(boundary))
        digest.update(b"\xff")
    return digest.hexdigest()


def rows_digest(rows: list[dict[str, object]]) -> str:
    certificate = [
        [int(row["z"]), row["representative_finite_field"], int(row["size"])]
        for row in rows
    ]
    return hashlib.sha256(
        json.dumps(certificate, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def direction_metadata() -> tuple[
    tuple[tuple[int, int], ...], tuple[int, ...], tuple[tuple[int, ...], ...]
]:
    directions = projective_directions(P)
    data = tuple(field_direction_data(P, direction) for direction in directions)
    types = tuple(int(eps) for eps, _labels in data)
    labels = tuple(tuple(int(value) for value in row) for _eps, row in data)
    if len(directions) != P + 1 or Counter(types) != {-1: 4, 1: 4}:
        raise AssertionError("unexpected p=7 projective direction/type census")
    if any(len(row) != Q or Counter(row) != {value: P for value in range(P)} for row in labels):
        raise AssertionError("a projective direction does not have seven fibres of size seven")
    return directions, types, labels


def pair_intersection_table(
    first: tuple[int, int],
    second: tuple[int, int],
    first_labels: tuple[int, ...],
    second_labels: tuple[int, ...],
) -> tuple[tuple[int, ...], ...]:
    """Point with prescribed labels for two distinct linear forms."""
    r0, s0 = first
    r1, s1 = second
    determinant = (r0 * s1 - r1 * s0) % P
    if determinant == 0:
        raise AssertionError("distinct projective directions were dependent")
    inverse = pow(determinant, -1, P)
    table = []
    for first_label in range(P):
        row = []
        for second_label in range(P):
            x = (first_label * s1 - s0 * second_label) * inverse % P
            y = (r0 * second_label - first_label * r1) * inverse % P
            point = x + P * y
            if (
                first_labels[point] != first_label
                or second_labels[point] != second_label
            ):
                raise AssertionError("linear-system intersection mislabeled")
            row.append(point)
        table.append(tuple(row))
    out = tuple(table)
    if len({point for row in out for point in row}) != Q:
        raise AssertionError("direction-pair intersections do not cover GF(49)")
    return out


def odd_fibre_masks(
    boundary: tuple[int, ...], labels: tuple[tuple[int, ...], ...]
) -> tuple[int, ...]:
    masks = []
    for direction_labels in labels:
        mask = 0
        for point in boundary:
            mask ^= 1 << direction_labels[point]
        masks.append(mask)
    return tuple(masks)


def enumerate_pair_transversals(
    directions: tuple[tuple[int, int], ...],
    types: tuple[int, ...],
    labels: tuple[tuple[int, ...], ...],
) -> tuple[
    dict[tuple[int, ...], int],
    dict[int, set[tuple[int, ...]]],
    list[dict[str, object]],
    dict[str, object],
]:
    """Generate all pair-transversal incidences and audit their quotient."""
    direction_pairs = tuple(itertools.combinations(range(P + 1), 2))
    pair_index = {pair: index for index, pair in enumerate(direction_pairs)}
    origin_bits: dict[tuple[int, ...], int] = {}
    pair_rows: list[dict[str, object]] = []
    incidence_total = 0

    for index, (first, second) in enumerate(direction_pairs):
        table = pair_intersection_table(
            directions[first], directions[second], labels[first], labels[second]
        )
        local_boundaries: set[tuple[int, ...]] = set()
        for permutation in itertools.permutations(range(P)):
            boundary = tuple(sorted(table[value][permutation[value]] for value in range(P)))
            if len(boundary) != FINITE_BOUNDARY_SIZE or len(set(boundary)) != FINITE_BOUNDARY_SIZE:
                raise AssertionError("a permutation graph was not a seven-point set")
            if boundary in local_boundaries:
                raise AssertionError("one direction pair generated a boundary twice")
            local_boundaries.add(boundary)
            old_bits = origin_bits.get(boundary, 0)
            bit = 1 << index
            if old_bits & bit:
                raise AssertionError("duplicate pair-transversal incidence")
            origin_bits[boundary] = old_bits | bit
        if len(local_boundaries) != math.factorial(P):
            raise AssertionError("a direction pair did not generate exactly 7! graphs")
        incidence_total += len(local_boundaries)
        pair_rows.append(
            {
                "pair_index": index,
                "directions": [first, second],
                "direction_types": [types[first], types[second]],
                "permutation_graphs": len(local_boundaries),
            }
        )

    if len(direction_pairs) != 28 or incidence_total != EXPECTED_PAIR_INCIDENCES:
        raise AssertionError("pair-transversal incidence census changed")

    branches = {z: set() for z in EXPECTED_BOUNDARY_COUNTS}
    z_histogram: Counter[int] = Counter()
    type_split_histograms: dict[int, Counter[tuple[int, int]]] = {
        z: Counter() for z in EXPECTED_BOUNDARY_COUNTS
    }
    multiplicity_histograms: dict[int, Counter[int]] = {
        z: Counter() for z in EXPECTED_BOUNDARY_COUNTS
    }
    pair_z_histograms = [Counter() for _pair in direction_pairs]

    for boundary, observed_bits in origin_bits.items():
        masks = odd_fibre_masks(boundary, labels)
        b_values = tuple(mask.bit_count() for mask in masks)
        if any(b not in (1, 3, 5, 7) for b in b_values):
            raise AssertionError("seven-point boundary had an impossible odd-fibre count")
        undetermined = tuple(index for index, b in enumerate(b_values) if b == P)
        z = len(undetermined)
        if z not in branches:
            raise AssertionError(f"unexpected z={z} in pair-transversal census")

        expected_bits = 0
        for pair in itertools.combinations(undetermined, 2):
            expected_bits |= 1 << pair_index[pair]
        if observed_bits != expected_bits:
            raise AssertionError("pair origins differ from all undetermined-direction pairs")
        if observed_bits.bit_count() != math.comb(z, 2):
            raise AssertionError("pair-transversal multiplicity is not C(z,2)")

        branches[z].add(boundary)
        z_histogram[z] += 1
        split = (
            sum(types[direction] == -1 for direction in undetermined),
            sum(types[direction] == 1 for direction in undetermined),
        )
        type_split_histograms[z][split] += 1
        multiplicity_histograms[z][observed_bits.bit_count()] += 1
        bits = observed_bits
        while bits:
            low_bit = bits & -bits
            pair_z_histograms[low_bit.bit_length() - 1][z] += 1
            bits ^= low_bit

    if dict(z_histogram) != EXPECTED_BOUNDARY_COUNTS:
        raise AssertionError(f"z census changed: {dict(z_histogram)}")
    if any(branches[z] & branches[w] for z, w in itertools.combinations(branches, 2)):
        raise AssertionError("z branches overlap")
    if sum(len(branch) for branch in branches.values()) != len(origin_bits):
        raise AssertionError("z branches do not partition generated boundaries")
    weighted_total = sum(len(branches[z]) * math.comb(z, 2) for z in branches)
    if weighted_total != incidence_total:
        raise AssertionError("weighted z census does not recover pair incidences")

    for row, histogram in zip(pair_rows, pair_z_histograms):
        if sum(histogram.values()) != math.factorial(P):
            raise AssertionError("per-pair z incidence census is incomplete")
        row["z_incidence_histogram"] = {
            str(z): count for z, count in sorted(histogram.items())
        }

    audit = {
        "generation_method": "28 direction pairs times all 7! permutation graphs",
        "uses_choose_49_7_scan": False,
        "unordered_direction_pairs": len(direction_pairs),
        "permutations_per_pair": math.factorial(P),
        "pair_transversal_incidences": incidence_total,
        "distinct_boundaries": len(origin_bits),
        "boundary_count_by_z": {
            str(z): z_histogram[z] for z in sorted(z_histogram)
        },
        "incidence_multiplicity_by_z": {
            str(z): {
                str(multiplicity): count
                for multiplicity, count in sorted(multiplicity_histograms[z].items())
            }
            for z in sorted(multiplicity_histograms)
        },
        "undetermined_type_split_histogram_by_z": {
            str(z): {
                f"minus_{minus}_plus_{plus}": count
                for (minus, plus), count in sorted(type_split_histograms[z].items())
            }
            for z in sorted(type_split_histograms)
        },
        "weighted_incidence_reconstruction": weighted_total,
        "every_pair_has_7_factorial_distinct_graphs": True,
        "pair_origin_masks_equal_all_undetermined_pairs": True,
        "multiplicity_equals_binomial_z_2": True,
        "z_branches_are_disjoint_and_complete": True,
        "boundary_sha256_by_z": {
            str(z): boundary_digest(branches[z]) for z in sorted(branches)
        },
    }
    return origin_bits, branches, pair_rows, audit


def multiplicative_order(value: int, mul) -> int:
    product = 1
    for order in range(1, Q):
        product = mul(product, value)
        if product == 1:
            return order
    raise AssertionError("nonzero GF(49) element had no multiplicative order")


def square_affine_semilinear_group(
    directions: tuple[tuple[int, int], ...],
    types: tuple[int, ...],
    labels: tuple[tuple[int, ...], ...],
) -> tuple[tuple[tuple[int, ...], ...], dict[str, object]]:
    """Construct and exactly audit x -> beta + alpha*x^(7^e), chi(alpha)=1."""
    q, mul, add, chi, frob, _norm, irreducible_a, irreducible_b = field_ctx(P)
    if q != Q:
        raise AssertionError("field context has the wrong order")
    squares = tuple(alpha for alpha in range(1, q) if chi(alpha) == 1)
    if len(squares) != (Q - 1) // 2:
        raise AssertionError("GF(49) square subgroup has the wrong size")

    parameter_count = 0
    by_permutation: dict[tuple[int, ...], tuple[int, int, int]] = {}
    for frobenius_power in (0, 1):
        conjugates = tuple(frob(u) if frobenius_power else u for u in range(q))
        for alpha in squares:
            linear = tuple(mul(alpha, conjugates[u]) for u in range(q))
            for beta in range(q):
                parameter_count += 1
                permutation = tuple(add(beta, linear[u]) for u in range(q))
                if len(set(permutation)) != q:
                    raise AssertionError("affine-semilinear parameter was not bijective")
                if permutation in by_permutation:
                    raise AssertionError("distinct affine-semilinear parameters coincided")
                by_permutation[permutation] = (beta, alpha, frobenius_power)

    permutations = tuple(sorted(by_permutation))
    identity = tuple(range(q))
    if (
        parameter_count != EXPECTED_GROUP_SIZE
        or len(permutations) != EXPECTED_GROUP_SIZE
        or identity not in by_permutation
    ):
        raise AssertionError("square affine-semilinear group size changed")

    square_generator = next(
        alpha for alpha in squares if multiplicative_order(alpha, mul) == len(squares)
    )
    generators = (
        tuple(add(1, u) for u in range(q)),
        tuple(add(P, u) for u in range(q)),
        tuple(mul(square_generator, u) for u in range(q)),
        tuple(frob(u) for u in range(q)),
    )
    permutation_set = set(permutations)
    for permutation in permutations:
        for generator in generators:
            composition = tuple(permutation[generator[u]] for u in range(q))
            if composition not in permutation_set:
                raise AssertionError("affine-semilinear set is not generator-closed")

    def subtract(u: int, v: int) -> int:
        return (u % P - v % P) % P + P * ((u // P - v // P) % P)

    for permutation in permutations:
        base = permutation[0]
        if any(
            chi(subtract(permutation[u], base)) != chi(u) for u in range(1, q)
        ):
            raise AssertionError("group element failed exact Paley-sign preservation")

    fibres = tuple(
        tuple(
            tuple(u for u in range(q) if labels[direction][u] == fibre)
            for fibre in range(P)
        )
        for direction in range(P + 1)
    )

    def partition_signature(rows: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
        return tuple(sorted(tuple(sorted(row)) for row in rows))

    partition_lookup = {
        partition_signature(fibres[direction]): direction
        for direction in range(P + 1)
    }
    if len(partition_lookup) != P + 1:
        raise AssertionError("projective direction partitions are not distinct")

    induced_actions: set[tuple[int, ...]] = set()
    for permutation in permutations:
        action = []
        for source in range(P + 1):
            transformed = partition_signature(
                tuple(tuple(permutation[u] for u in fibre) for fibre in fibres[source])
            )
            target = partition_lookup.get(transformed)
            if target is None:
                raise AssertionError("group element did not permute affine directions")
            if types[target] != types[source]:
                raise AssertionError("group element changed a direction's quadratic type")
            action.append(target)
        if sorted(action) != list(range(P + 1)):
            raise AssertionError("induced direction action is not a permutation")
        induced_actions.add(tuple(action))

    permutation_digest = hashlib.sha256()
    for permutation in permutations:
        permutation_digest.update(bytes(permutation))
        permutation_digest.update(b"\xff")
    action_digest = hashlib.sha256()
    for action in sorted(induced_actions):
        action_digest.update(bytes(action))
        action_digest.update(b"\xff")

    audit = {
        "field": "GF(7)[w] with w^2 = a*w + b",
        "field_encoding": "c0 + 7*c1",
        "irreducible_polynomial_parameters": {
            "a": irreducible_a,
            "b": irreducible_b,
        },
        "action": "x -> beta + alpha*x^(7^e)",
        "translations": Q,
        "square_multipliers": len(squares),
        "frobenius_powers": 2,
        "parameter_count": parameter_count,
        "unique_permutation_count": len(permutations),
        "group_size": len(permutations),
        "square_subgroup_generator": square_generator,
        "square_subgroup_generator_order": multiplicative_order(square_generator, mul),
        "generator_count": len(generators),
        "closed_under_right_composition_by_generators": True,
        "all_parameterizations_unique": True,
        "all_maps_bijective": True,
        "identity_present": True,
        "orbit_of_zero_size": len({permutation[0] for permutation in permutations}),
        "all_preserve_paley_difference_signs": True,
        "all_permute_projective_directions": True,
        "all_preserve_direction_types": True,
        "induced_direction_action_count": len(induced_actions),
        "permutation_sha256": permutation_digest.hexdigest(),
        "induced_direction_action_sha256": action_digest.hexdigest(),
    }
    if audit["orbit_of_zero_size"] != Q:
        raise AssertionError("affine group is not transitive on finite points")
    return permutations, audit


def representative_record(
    representative: tuple[int, ...],
    orbit_size: int,
    z: int,
    branch_orbit_index: int,
    origin_bits: dict[tuple[int, ...], int],
    direction_pairs: tuple[tuple[int, int], ...],
    directions: tuple[tuple[int, int], ...],
    types: tuple[int, ...],
    labels: tuple[tuple[int, ...], ...],
) -> dict[str, object]:
    masks = odd_fibre_masks(representative, labels)
    undetermined = tuple(
        direction for direction, mask in enumerate(masks) if mask.bit_count() == P
    )
    if len(undetermined) != z:
        raise AssertionError("representative z does not match its branch")
    source_pairs = [
        index
        for index in range(len(direction_pairs))
        if origin_bits[representative] & (1 << index)
    ]
    if [direction_pairs[index] for index in source_pairs] != list(
        itertools.combinations(undetermined, 2)
    ):
        raise AssertionError("representative pair origins are incomplete")

    direction_rows = []
    for index, (direction, eps, mask) in enumerate(zip(directions, types, masks)):
        direction_rows.append(
            {
                "direction_index": index,
                "linear_form": list(direction),
                "quadratic_type": eps,
                "odd_fibre_mask": mask,
                "odd_fibres": [fibre for fibre in range(P) if mask & (1 << fibre)],
                "b": mask.bit_count(),
                "undetermined": mask == (1 << P) - 1,
            }
        )
    return {
        "branch_orbit_index": branch_orbit_index,
        "z": z,
        "representative_finite_field": list(representative),
        "representative_vertices": [0, *(point + 1 for point in representative)],
        "size": orbit_size,
        "stabilizer_size": EXPECTED_GROUP_SIZE // orbit_size,
        "direction_masks": list(masks),
        "direction_types": list(types),
        "undetermined_directions": list(undetermined),
        "undetermined_direction_types": [types[index] for index in undetermined],
        "undetermined_type_counts": {
            "-1": sum(types[index] == -1 for index in undetermined),
            "1": sum(types[index] == 1 for index in undetermined),
        },
        "pair_transversal_multiplicity": origin_bits[representative].bit_count(),
        "pair_transversal_source_pair_indices": source_pairs,
        "pair_transversal_source_pairs": [
            list(direction_pairs[index]) for index in source_pairs
        ],
        "direction_rows": direction_rows,
    }


def decompose_orbits(
    branches: dict[int, set[tuple[int, ...]]],
    origin_bits: dict[tuple[int, ...], int],
    permutations: tuple[tuple[int, ...], ...],
    directions: tuple[tuple[int, int], ...],
    types: tuple[int, ...],
    labels: tuple[tuple[int, ...], ...],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    direction_pairs = tuple(itertools.combinations(range(P + 1), 2))
    rows: list[dict[str, object]] = []
    size_histograms: dict[int, Counter[int]] = {}
    covered_by_z: dict[int, set[tuple[int, ...]]] = {}
    z2_orbit_type_split: Counter[str] = Counter()

    for z in sorted(branches):
        branch = branches[z]
        remaining = set(branch)
        covered: set[tuple[int, ...]] = set()
        branch_rows = []
        while remaining:
            representative = min(remaining)
            orbit = {
                tuple(sorted(permutation[point] for point in representative))
                for permutation in permutations
            }
            if min(orbit) != representative:
                raise AssertionError("selected orbit representative was not canonical")
            if not orbit <= branch:
                raise AssertionError("group orbit left its exact z branch")
            if orbit & covered:
                raise AssertionError("computed group orbits overlap")
            if EXPECTED_GROUP_SIZE % len(orbit):
                raise AssertionError("orbit size does not divide group order")
            branch_rows.append(
                representative_record(
                    representative,
                    len(orbit),
                    z,
                    len(branch_rows),
                    origin_bits,
                    direction_pairs,
                    directions,
                    types,
                    labels,
                )
            )
            covered.update(orbit)
            remaining.difference_update(orbit)

        if covered != branch:
            raise AssertionError("orbit partition does not cover its z branch")
        covered_by_z[z] = covered
        size_histogram = Counter(int(row["size"]) for row in branch_rows)
        if len(branch_rows) != EXPECTED_ORBIT_COUNTS[z]:
            raise AssertionError(f"z={z} orbit count changed: {len(branch_rows)}")
        if dict(size_histogram) != EXPECTED_ORBIT_SIZE_HISTOGRAMS[z]:
            raise AssertionError(f"z={z} orbit-size histogram changed: {dict(size_histogram)}")
        if sum(size * count for size, count in size_histogram.items()) != len(branch):
            raise AssertionError("orbit-size sum does not recover branch census")
        if z == 2:
            for row in branch_rows:
                undetermined_types = tuple(int(value) for value in row["undetermined_direction_types"])
                classification = (
                    "same_type" if len(set(undetermined_types)) == 1 else "split_type"
                )
                z2_orbit_type_split[classification] += 1
            if dict(z2_orbit_type_split) != EXPECTED_Z2_ORBIT_TYPE_SPLIT:
                raise AssertionError(
                    "z=2 orbit direction-type split changed: "
                    f"{dict(z2_orbit_type_split)}"
                )
        size_histograms[z] = size_histogram
        rows.extend(branch_rows)

    if len(rows) != sum(EXPECTED_ORBIT_COUNTS.values()):
        raise AssertionError("total orbit count changed")
    for global_index, row in enumerate(rows):
        row["orbit_index"] = global_index

    audit = {
        "orbit_count": len(rows),
        "orbit_count_by_z": {
            str(z): sum(int(row["z"]) == z for row in rows) for z in sorted(branches)
        },
        "orbit_size_histogram_by_z": {
            str(z): {
                str(size): count for size, count in sorted(size_histograms[z].items())
            }
            for z in sorted(size_histograms)
        },
        "orbit_size_sum_by_z": {
            str(z): sum(size * count for size, count in size_histograms[z].items())
            for z in sorted(size_histograms)
        },
        "z2_orbit_direction_type_split": {
            key: z2_orbit_type_split[key] for key in ("same_type", "split_type")
        },
        "z2_orbit_direction_type_split_matches_48_plus_44": True,
        "all_orbit_sizes_divide_group_order": True,
        "orbits_are_pairwise_disjoint": True,
        "orbits_cover_every_generated_boundary": all(
            covered_by_z[z] == branches[z] for z in branches
        ),
        "representatives_are_lexicographically_canonical": True,
        "orbit_certificate_sha256": rows_digest(rows),
    }
    return rows, audit


def generate() -> dict[str, object]:
    started = time.time()
    directions, types, labels = direction_metadata()
    origin_bits, branches, pair_rows, pair_audit = enumerate_pair_transversals(
        directions, types, labels
    )
    permutations, group_audit = square_affine_semilinear_group(
        directions, types, labels
    )
    orbits, orbit_audit = decompose_orbits(
        branches, origin_bits, permutations, directions, types, labels
    )

    required_totals = {
        str(z): len(branches[z]) for z in sorted(EXPECTED_BOUNDARY_COUNTS)
    }
    required_orbits = {
        str(z): sum(int(row["z"]) == z for row in orbits)
        for z in sorted(EXPECTED_ORBIT_COUNTS)
    }
    if required_totals != {str(z): count for z, count in EXPECTED_BOUNDARY_COUNTS.items()}:
        raise AssertionError("required boundary totals failed")
    if required_orbits != {str(z): count for z, count in EXPECTED_ORBIT_COUNTS.items()}:
        raise AssertionError("required orbit totals failed")

    return {
        "experiment": "p7_infinity7_positive_zge2_orbits",
        "status": "complete_exact_pair_transversal_orbit_census",
        "p": P,
        "c_H": 1,
        "boundary": "infinity plus seven finite GF(49) points",
        "finite_boundary_size": FINITE_BOUNDARY_SIZE,
        "z_definition": "number of projective directions with odd-fibre mask 127 (b=7)",
        "scope_z": sorted(EXPECTED_BOUNDARY_COUNTS),
        "directions": [
            {
                "direction_index": index,
                "linear_form": list(direction),
                "quadratic_type": types[index],
            }
            for index, direction in enumerate(directions)
        ],
        "census": {
            "pair_transversal_incidences": EXPECTED_PAIR_INCIDENCES,
            "distinct_boundaries": sum(EXPECTED_BOUNDARY_COUNTS.values()),
            "boundary_count_by_z": required_totals,
            "orbit_count": sum(EXPECTED_ORBIT_COUNTS.values()),
            "orbit_count_by_z": required_orbits,
            "orbit_size_histogram_by_z": orbit_audit["orbit_size_histogram_by_z"],
        },
        "pair_transversal_audit": pair_audit,
        "direction_pair_audits": pair_rows,
        "group_audit": group_audit,
        "orbit_audit": orbit_audit,
        "all_required_audits_passed": True,
        "orbits": orbits,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the exact positive p=7 infinity-plus-seven z=2,3,7 "
            "boundary orbits from pair transversals."
        )
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = generate()
    atomic_write(args.output, result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.output),
                "census": result["census"],
                "group_size": result["group_audit"]["group_size"],
                "all_required_audits_passed": result["all_required_audits_passed"],
                "elapsed_seconds": result["elapsed_seconds"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
