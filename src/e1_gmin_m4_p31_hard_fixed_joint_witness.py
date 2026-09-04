"""Exact joint-gate witness and boundary rejection at ``p=31``.

The displayed sixteen localized-Mobius halves simultaneously realize the
hard-fixed parallel profile, the required two-bit selector correction, and a
physical direction-zero cancellation.  Thus those three necessary gates are
jointly feasible; treating them as mutually incompatible cannot close E1.

The witness nevertheless does not lift to the required common graph.  Its
prescribed collision shape occurs in exactly one pair of halves, with thirty
scalar-equivalent centre lifts.  For every lift, the full vertex-boundary
system for the other fourteen centres and the fixed antipodal edge has
coefficient rank 225 and augmented rank 226 over GF(2).  This excludes only
this sixteen-half witness.  It does not close the hard-fixed family, E1, or
residual (ii).
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from itertools import combinations

from e1_gmin_m4_inversion_antisymmetric_radon import (
    _negative_edge,
    edge_radon_image,
)
from e1_gmin_m4_mobius_half_symmetric import paley_edge_sign
from e1_gmin_m4_p31_mobius_collision_key import (
    annihilator_direction_index,
)
from e1_gmin_m4_p31_top_mobius_boundary_parity import (
    half_kernel_sigma_formula,
)
from scripts import residual_branch_c_auxiliary_transverse_gpu as auxiliary
from scripts.residual_branch_c_center_boundary_gf2 import (
    BoundaryColumnCache,
    boundary_system,
    boundary_system_columns,
)


P = 31
FIXED_DIRECTION_INDEX = 0
HALF_CHOICES = (
    (0, 15, 3),
    (1, 9, 9),
    (2, 24, 21),
    (3, 26, 7),
    (7, 4, 22),
    (9, 22, 23),
    (10, 16, 27),
    (15, 21, 6),
    (16, 2, 1),
    (21, 12, 14),
    (22, 28, 20),
    (24, 3, 30),
    (28, 29, 26),
    (29, 31, 27),
    (30, 10, 21),
    (31, 7, 27),
)
EXPECTED_RAW_PROFILE = (
    15, 15, 14, 14, 15, 16, 16, 14,
    16, 14, 14, 16, 15, 16, 16, 14,
    14, 16, 16, 16, 16, 14, 14, 16,
    14, 16, 15, 16, 14, 14, 15, 14,
)
EXPECTED_AGGREGATE_SIGNATURE = 0x01000401
EXPECTED_CORRECTION_SUPPORT = (10, 24)
COLLISION_HALF_INDICES = (1, 5)
COLLISION_ENDPOINT_SUPPORT = (10, 24)
CANONICAL_COLLISION_CENTERS = (1, 22)
CANONICAL_COLLISION_ORBIT = ((3, 9), (3, 27))
CANONICAL_COLLISION_COEFFICIENTS = (1, -1)
CLEAN_CENTERS = (7, 1, 12, 6, 12, 22, 28, 26, 12, 3, 21, 3, 12, 25, 2, 16)
CLEAN_FIXED_EDGE = ((0, 1), (0, 30))
EXPECTED_CLEAN_GRAPH_SHA256 = (
    "6e924c3acb493799f7951a6ab75e22a2628d452d4654643944c1f3871a75a6a4"
)
EXPECTED_BOUNDARY_MATRIX_SHA256 = (
    "c32f2cecfe8c378f290fc3cfae70e2fa2e85046bae7981572c6537c1147f4207"
)


def _design() -> tuple[auxiliary.HalfChoice, ...]:
    return tuple(auxiliary.HalfChoice(*row) for row in HALF_CHOICES)


def _endpoint_support(orbit) -> tuple[int, int]:
    support = tuple(
        sorted(annihilator_direction_index(P, point) for point in orbit)
    )
    if len(set(support)) != 2:
        raise ArithmeticError("a collision orbit lost its two endpoint rays")
    return support


def _aggregate_signature(design: tuple[auxiliary.HalfChoice, ...]) -> int:
    signature = 0
    for choice in design:
        for kernel_index, kernel in enumerate(auxiliary.DIRECTIONS):
            if half_kernel_sigma_formula(choice.target, choice.auxiliary, kernel) == -1:
                signature ^= 1 << kernel_index
    return signature


def _prescribed_collision_seeds(design, orbit_cache) -> tuple[dict[str, object], ...]:
    seeds = []
    matching_pairs = set()
    for first, second in combinations(range(len(design)), 2):
        for first_center in range(1, P):
            first_map = orbit_cache[first][first_center - 1]
            for second_center in range(1, P):
                second_map = orbit_cache[second][second_center - 1]
                for orbit in set(first_map) & set(second_map):
                    if (
                        first_map[orbit] == -second_map[orbit]
                        and auxiliary.frozen._spatial_direction_index(orbit)
                        == FIXED_DIRECTION_INDEX
                        and _endpoint_support(orbit) == COLLISION_ENDPOINT_SUPPORT
                    ):
                        matching_pairs.add((first, second))
                        seeds.append(
                            {
                                "half_indices": (first, second),
                                "centers": (first_center, second_center),
                                "orbit": orbit,
                                "coefficients": (
                                    first_map[orbit],
                                    second_map[orbit],
                                ),
                            }
                        )
    seeds.sort(
        key=lambda row: (row["half_indices"], row["centers"], row["orbit"])
    )
    if matching_pairs != {COLLISION_HALF_INDICES} or len(seeds) != P - 1:
        raise ArithmeticError("the unique prescribed collision pair changed")
    return tuple(seeds)


def _clean_graph_replay(design, orbit_cache) -> dict[str, object]:
    total = Counter()
    for half_index, center in enumerate(CLEAN_CENTERS):
        total.update(orbit_cache[half_index][center - 1])
    cancelled = tuple(sorted(orbit for orbit, value in total.items() if value == 0))
    surviving = {orbit: value for orbit, value in total.items() if value}
    if (
        cancelled != (CANONICAL_COLLISION_ORBIT,)
        or len(surviving) != 478
        or any(abs(value) != 1 for value in surviving.values())
    ):
        raise ArithmeticError("the clean physical collision replay changed")
    graph = [
        orbit if coefficient == 1 else _negative_edge(P, orbit)
        for orbit, coefficient in surviving.items()
    ]
    if CLEAN_FIXED_EDGE in graph:
        raise ArithmeticError("the fixed edge collided with a surviving half edge")
    graph = tuple(sorted((*graph, CLEAN_FIXED_EDGE)))
    if len(graph) != 479 or len(set(graph)) != 479:
        raise ArithmeticError("the physical witness is not a simple 479-edge graph")

    image = edge_radon_image(
        P, {edge: paley_edge_sign(P, edge) for edge in graph}
    )
    parallel = tuple(
        auxiliary.DIRECTION_SIGNS[index] * image.get(("P", index), 0)
        for index in range(P + 1)
    )
    expected = tuple(
        value - int(index == FIXED_DIRECTION_INDEX)
        for index, value in enumerate(EXPECTED_RAW_PROFILE)
    )
    digest = hashlib.sha256(
        json.dumps(graph, separators=(",", ":")).encode()
    ).hexdigest()
    if parallel != expected or digest != EXPECTED_CLEAN_GRAPH_SHA256:
        raise ArithmeticError("the physical graph profile or digest changed")
    return {
        "centers": CLEAN_CENTERS,
        "fixed_edge": CLEAN_FIXED_EDGE,
        "cancelled_orbit": cancelled[0],
        "graph_edge_count": len(graph),
        "graph_sha256": digest,
        "final_parallel_profile": parallel,
        "exact_simple_graph_replay": True,
    }


def joint_witness_boundary_certificate() -> dict[str, object]:
    """Replay joint feasibility and the exact GF(2) rejection of this witness."""
    design = _design()
    raw_profile = auxiliary.raw_parallel_profile(design)
    aggregate_signature = _aggregate_signature(design)
    correction_signature = aggregate_signature ^ (1 << FIXED_DIRECTION_INDEX)
    correction_support = tuple(
        index for index in range(P + 1) if correction_signature >> index & 1
    )
    if (
        tuple(choice.target_index for choice in design) != auxiliary.HARD
        or len({choice.auxiliary_index for choice in design}) != len(design)
        or raw_profile != EXPECTED_RAW_PROFILE
        or aggregate_signature != EXPECTED_AGGREGATE_SIGNATURE
        or correction_support != EXPECTED_CORRECTION_SUPPORT
    ):
        raise ArithmeticError("the hard-fixed profile/signature witness changed")

    orbit_cache = auxiliary._orbit_cache(design)
    seeds = _prescribed_collision_seeds(design, orbit_cache)
    canonical_seed = seeds[0]
    if (
        canonical_seed["centers"] != CANONICAL_COLLISION_CENTERS
        or canonical_seed["orbit"] != CANONICAL_COLLISION_ORBIT
        or canonical_seed["coefficients"] != CANONICAL_COLLISION_COEFFICIENTS
    ):
        raise ArithmeticError("the canonical physical collision changed")
    graph = _clean_graph_replay(design, orbit_cache)

    halves = tuple((choice.target, choice.auxiliary) for choice in design)
    boundary_cache = BoundaryColumnCache()
    systems = tuple(
        boundary_system_columns(
            halves, seed, FIXED_DIRECTION_INDEX, boundary_cache
        )
        for seed in seeds
    )
    rank_histogram = Counter(
        (row["coefficient_rank"], row["augmented_rank"], row["consistent"])
        for row in systems
    )
    matrix_hashes = {
        row["matrix_augmented_columns_sha256"] for row in systems
    }
    if (
        rank_histogram != {(225, 226, False): P - 1}
        or systems[0]["matrix_augmented_columns_sha256"]
        != EXPECTED_BOUNDARY_MATRIX_SHA256
    ):
        raise ArithmeticError("the full vertex-boundary rejection changed")

    independent_row_replay = boundary_system(
        halves, canonical_seed, FIXED_DIRECTION_INDEX
    )
    if (
        independent_row_replay["coefficient_rank"] != 225
        or independent_row_replay["augmented_rank"] != 226
        or independent_row_replay["consistent"] is not False
        or independent_row_replay["contradiction_equation_count"] == 0
    ):
        raise ArithmeticError("row/column boundary eliminations disagree")

    return {
        "schema": "e1_p31_hard_fixed_joint_witness_boundary_v1",
        "classification": (
            "exact E1/common-graph witness for three joint gates, then exact "
            "GF(2) boundary rejection of that witness"
        ),
        "p": P,
        "half_choices_target_auxiliary_scale": HALF_CHOICES,
        "raw_parallel_profile": raw_profile,
        "aggregate_signature_hex": f"{aggregate_signature:08x}",
        "correction_signature_support": correction_support,
        "prescribed_collision_pair_half_indices": COLLISION_HALF_INDICES,
        "prescribed_collision_lift_count": len(seeds),
        "joint_profile_signature_physical_collision_feasible": True,
        "clean_physical_graph": graph,
        "boundary_system": {
            "collision_lifts_checked": len(systems),
            "variable_count": systems[0]["variable_count"],
            "equation_count": systems[0]["equation_count"],
            "coefficient_rank": systems[0]["coefficient_rank"],
            "augmented_rank": systems[0]["augmented_rank"],
            "all_collision_lifts_inconsistent": True,
            "canonical_matrix_augmented_columns_sha256": systems[0][
                "matrix_augmented_columns_sha256"
            ],
            "distinct_scalar_matrix_hash_count": len(matrix_hashes),
            "independent_row_contradiction_equation_count": independent_row_replay[
                "contradiction_equation_count"
            ],
            "independent_row_matrix_sha256": independent_row_replay[
                "matrix_augmented_rows_sha256"
            ],
        },
        "displayed_sixteen_half_witness_excluded": True,
        "scope": (
            "excludes only the displayed hard-fixed p=31 sixteen-half witness; "
            "other labelled half families and arbitrary lifts remain open"
        ),
        "e1_closed": False,
        "residual_ii_closed": False,
    }
