#!/usr/bin/env python3
"""Prop. 15.715 -- close positive p7 infinity+7 with z=1.

There are exactly four mean allocations for each z=1 boundary.  A complete
CUDA scan projects them onto 23 exact mod-seven dependencies, and the host
checks every projected boundary against the complete catalogs on all 135
dependencies.  None of the 6,324,528 z=1 boundaries survives.
"""
from __future__ import annotations

import functools
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15714 import (
    _array_sha256,
    p7_positive_infinity_plus_seven_z0_exclusion,
)


ROOT = Path(__file__).resolve().parents[1]
SCAN = ROOT / "evidence" / "p7_infinity7_positive_z1_v100.json"
AUDIT = ROOT / "evidence" / "p7_infinity7_positive_z1_v100_audit.json"
PROJECTED_RANKS = (
    ROOT / "evidence" / "p7_infinity7_positive_z1_projected_ranks.json"
)

SCAN_SHA256 = "a965cb6242118daba652dad2053dec06d4dccb59729d7ab83fb65a81fb487c12"
AUDIT_SHA256 = "edcfb01b515b2f0141f6ef652bf25624067bcefb7dbc49984ddaf5d2d97e4176"
PROJECTED_RANKS_SHA256 = (
    "e44edb4c1689af9ca859f4410ae33f24a1f0613217580abe0781858aa3dc128d"
)
EXPECTED_PROJECTED_RANK_SHA256 = (
    "23de1f85d34f641d06279e8cdbc17fe6615fcd98198885e18f695d1812982b4c"
)
PRESERVED_PRIMARY_FULL_OUTPUT_SHA256 = (
    "75b37dcea1c7862677738f9f5e67d2b47f67c73b1e2b207064500dc891942464"
)

EXPECTED_INPUT_SHA256 = {
    "equation_matrix": "32b378e8bd6c55deb9b6b546c73ee869b69a5e0d7037f0ed474be5ae882fbc1a",
    "left_dependencies": "0405fad25d2295ed722bd8ee15ebd6592907d8abdf9d97ddfa866500816dbad2",
    "base": "1bd5ca1015bbb735222fe100b1e8b41bb54b13df2aa8667480b11108ad4cea05",
    "floor_tables": "2b66d3b0184b8bb3cc85452ebffa91a819aa26dcf83eac59cae2091c08277656",
    "finite_labels": "96c23cbbee81215029c045ddc326ef7950db23fdc61036d4714f798d9db8e895",
    "direction_types": "a71e5bf6a9f0ea29944e99ec7a721144096416ced212c5f534c2f4a633d48fd5",
    "zero_tables": "72642567c6d5dfe031174a2690df4f0be1501921e353d923c934ca7f10e07d79",
    "selected_rows": "a3f3d4d6546d519e09bd0ff6f7b957b8b59327e90d9a2c7b8fd2d645d0325d09",
    "catalog_keys": "6350a652f1d7fa725883171874813fd0a370ddcb894436ac7edf25371b2a88bf",
    "catalog_tails": "c9cad7108ff2344cb681121314066b9fc5823cfa15aceb9ef81af915856f0b3c",
    "catalog_counts": "8a453446b3f6c7a85f19e4782f0b9add3dff1e242a255ea256f59abf65a719a6",
}

EXPECTED_SOURCE_SHA256 = {
    "scripts/p7_infinity7_positive_z1_mod7_gpu.py": (
        "2e03ddbc8f27a0e00fef2ee603571c6fdd50a4a42e9785806c61f2237e082019"
    ),
    "scripts/p7_infinity7_positive_z0_mod7_gpu.py": (
        "f748005d9be3286094cf7693941544d762f8bef3ae30a27c19b32fb781e7d951"
    ),
    "scripts/p7_size6_positive_infinity_mod7_gpu.py": (
        "d0e54d2749a1fcd2841674134301fef241acbed2abdc5777aedfc5c36e87330a"
    ),
    "scripts/p7_unsaturated_mod7_batch.py": (
        "598c5012a080677df0bb6cde00d3beb256ad5cf9a7cdf36cb5dc46591e6a6f27"
    ),
    "scripts/p7_unsaturated_modular_catalog_filter.py": (
        "f9b2781984ab3e2336977d43b657fe337bb09b37baea600b6fdb1f94483d135a"
    ),
    "src/e1_gmin_m4_prop15632.py": (
        "eda17b867dfd9654eb69d5b9f8b6dbb1b9791dd4ef79acc3060ddbf980822e0c"
    ),
}


def _read_pinned_json(path: Path, expected_sha256: str) -> dict[str, object]:
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise ArithmeticError(f"p7 positive z1 evidence bytes changed: {path.name}")
    return json.loads(raw)


@functools.cache
def _rebuild_z1_inputs() -> dict[str, object]:
    from e1_gmin_m4_prop15632 import field_direction_data, projective_directions
    from p7_infinity7_positive_z1_mod7_gpu import (
        catalog_tables,
        independent_dependency_rows,
        zero_contributions,
    )
    from p7_size6_positive_infinity_mod7_gpu import dependency_tables, finite_labels
    from p7_unsaturated_modular_catalog_filter import equation_matrix, left_dependencies

    matrix = equation_matrix()
    rank, dependencies = left_dependencies(matrix, 7)
    base, floor_tables, linear = dependency_tables()
    labels = finite_labels()
    direction_types = np.asarray(
        [field_direction_data(7, direction)[0] for direction in projective_directions(7)],
        dtype=np.int8,
    )
    zero_tables = zero_contributions(dependencies)
    selected_rows = np.stack(
        [independent_dependency_rows(dependencies, direction) for direction in range(8)]
    )
    catalog_keys, catalog_tails, catalog_counts, catalog_histogram = catalog_tables(
        dependencies, selected_rows
    )
    return {
        "matrix": matrix,
        "rank": rank,
        "dependencies": dependencies,
        "base": base,
        "floor_tables": floor_tables,
        "linear": linear,
        "labels": labels,
        "direction_types": direction_types,
        "zero_tables": zero_tables,
        "selected_rows": selected_rows,
        "catalog_keys": catalog_keys,
        "catalog_tails": catalog_tails,
        "catalog_counts": catalog_counts,
        "catalog_histogram": catalog_histogram,
    }


@functools.cache
def _current_z1_input_integrity() -> dict[str, object]:
    """Rebuild and pin the projection and full-dependency mathematical inputs."""
    inputs = _rebuild_z1_inputs()
    actual = {
        "equation_matrix": _array_sha256(inputs["matrix"]),
        "left_dependencies": _array_sha256(inputs["dependencies"]),
        "base": _array_sha256(inputs["base"]),
        "floor_tables": _array_sha256(inputs["floor_tables"]),
        "finite_labels": _array_sha256(inputs["labels"]),
        "direction_types": _array_sha256(inputs["direction_types"]),
        "zero_tables": _array_sha256(inputs["zero_tables"]),
        "selected_rows": _array_sha256(inputs["selected_rows"]),
        "catalog_keys": _array_sha256(inputs["catalog_keys"]),
        "catalog_tails": _array_sha256(inputs["catalog_tails"]),
        "catalog_counts": _array_sha256(inputs["catalog_counts"]),
    }
    if actual != EXPECTED_INPUT_SHA256:
        raise ArithmeticError("p7 positive z1 mathematical input arrays changed")
    if (
        int(inputs["rank"]) != 147
        or inputs["catalog_histogram"] != {"1764": 232, "2233": 280}
    ):
        raise ArithmeticError("p7 positive z1 input reconstruction changed")
    source_sha256 = {
        relative: hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        for relative in EXPECTED_SOURCE_SHA256
    }
    if source_sha256 != EXPECTED_SOURCE_SHA256:
        raise ArithmeticError("p7 positive z1 scan source changed")
    return {
        "array_sha256": actual,
        "source_sha256": source_sha256,
        "recomputed_from_current_checkout": True,
    }


def _unrank_boundary(rank: int) -> tuple[int, ...]:
    """Small CPU combinadic decoder kept separate from the CUDA kernel."""
    out = []
    next_value = 0
    for position in range(7):
        remaining = 6 - position
        for candidate in range(next_value, 49 - remaining):
            ways = math.comb(49 - candidate - 1, remaining)
            if rank < ways:
                out.append(candidate)
                next_value = candidate + 1
                break
            rank -= ways
    if len(out) != 7:
        raise ArithmeticError("p7 z1 projected rank failed CPU combinadic decoding")
    return tuple(out)


def _packed_target_present(
    inputs: dict[str, object], direction: int, mask: int, target: np.ndarray
) -> bool:
    selected = inputs["selected_rows"][direction]
    values = np.asarray(target[selected], dtype=np.int64)
    key = 0
    for value in values[:22]:
        key = 7 * key + int(value)
    tail = int(values[22])
    count = int(inputs["catalog_counts"][direction, mask])
    keys = inputs["catalog_keys"][direction, mask, :count]
    tails = inputs["catalog_tails"][direction, mask, :count]
    return bool(np.any((keys == key) & (tails == tail)))


def _cpu_projection_passes(rank: int) -> bool:
    """Replay the 23-coordinate projection for one stored rank on the CPU."""
    inputs = _rebuild_z1_inputs()
    boundary = _unrank_boundary(rank)
    labels = inputs["labels"]
    masks = []
    for direction in range(8):
        mask = 0
        for point in boundary:
            mask ^= 1 << int(labels[direction, point])
        masks.append(mask)
    undetermined = [
        direction for direction, mask in enumerate(masks) if mask.bit_count() == 7
    ]
    if len(undetermined) != 1:
        return False
    u = undetermined[0]

    fixed = inputs["base"].astype(np.int64).copy()
    for direction, mask in enumerate(masks):
        if direction != u:
            fixed += inputs["floor_tables"][direction, mask]
    if _packed_target_present(inputs, u, masks[u], (-fixed) % 7):
        return True

    for elevated in range(8):
        if (
            elevated == u
            or inputs["direction_types"][elevated]
            != inputs["direction_types"][u]
        ):
            continue
        fixed = inputs["base"].astype(np.int64).copy() + inputs["zero_tables"][u]
        for direction, mask in enumerate(masks):
            if direction != u and direction != elevated:
                fixed += inputs["floor_tables"][direction, mask]
        if _packed_target_present(inputs, elevated, masks[elevated], (-fixed) % 7):
            return True
    return False


def _validate_projected_rank_certificate(
    scan: dict[str, object], audit: dict[str, object]
) -> tuple[list[int], dict[str, object]]:
    certificate = _read_pinned_json(PROJECTED_RANKS, PROJECTED_RANKS_SHA256)
    ranks = [int(rank) for rank in certificate["projected_survivor_ranks"]]
    digest = hashlib.sha256(
        json.dumps(ranks, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if (
        certificate["experiment"]
        != "p7_infinity7_positive_z1_projected_survivor_rank_certificate"
        or certificate["status"]
        != "full_rank_list_recovered_from_preserved_primary_run_output"
        or int(certificate["p"]) != 7
        or int(certificate["rank_count"]) != 1_326
        or certificate["rank_encoding"]
        != "sorted zero-based combinadic ranks in C(49,7)"
        or certificate["rank_digest_encoding"]
        != "SHA-256 of compact ASCII JSON list with separators comma and colon"
        or certificate["primary_summary_evidence"] != str(SCAN.relative_to(ROOT))
        or certificate["primary_summary_sha256"] != SCAN_SHA256
        or certificate["different_grid_rerun_summary_evidence"]
        != str(AUDIT.relative_to(ROOT))
        or certificate["different_grid_rerun_summary_sha256"] != AUDIT_SHA256
        or certificate["preserved_primary_full_output_sha256"]
        != PRESERVED_PRIMARY_FULL_OUTPUT_SHA256
        or certificate["summaries_report_same_digest"] is not True
        or certificate["independent_implementation_validated"] is not False
        or len(ranks) != 1_326
        or len(set(ranks)) != 1_326
        or ranks != sorted(ranks)
        or not all(0 <= rank < math.comb(49, 7) for rank in ranks)
        or digest != EXPECTED_PROJECTED_RANK_SHA256
        or certificate["projected_survivor_rank_sha256"] != digest
        or scan["projected_survivor_rank_sha256"] != digest
        or audit["projected_survivor_rank_sha256"] != digest
        or scan["first_projected_survivor_ranks"] != ranks[:64]
        or audit["first_projected_survivor_ranks"] != ranks[:64]
    ):
        raise ArithmeticError("p7 positive z1 full projected-rank certificate changed")
    if not all(_cpu_projection_passes(rank) for rank in ranks):
        raise ArithmeticError("stored p7 z1 projected rank fails current CPU projection replay")
    return ranks, {
        "rank_count": len(ranks),
        "rank_sha256": digest,
        "all_ranks_replayed_against_current_cpu_projection": True,
        "recovery_status": certificate["status"],
        "preserved_primary_full_output_sha256": (
            PRESERVED_PRIMARY_FULL_OUTPUT_SHA256
        ),
        "independent_implementation_validated": False,
    }


def p7_positive_infinity_plus_seven_z1_exclusion() -> dict[str, object]:
    previous = p7_positive_infinity_plus_seven_z0_exclusion()
    scan = _read_pinned_json(SCAN, SCAN_SHA256)
    audit = _read_pinned_json(AUDIT, AUDIT_SHA256)
    input_integrity = _current_z1_input_integrity()
    projected_ranks, projected_rank_validation = _validate_projected_rank_certificate(
        scan, audit
    )
    linear = scan["linear_system"]
    if (
        scan["experiment"] != "p7_infinity7_positive_z1_mod7_gpu"
        or scan["status"] != "complete_projected_then_exact_mod_seven_z1_exhaustion"
        or int(scan["all_boundaries"]) != 85_900_584
        or int(scan["checked_boundaries"]) != 85_900_584
        or int(scan["z1_boundaries"]) != 6_324_528
        or int(scan["mean_allocation_count_per_boundary"]) != 4
        or int(scan["projected_dependency_count"]) != 23
        or int(scan["projected_survivors"]) != 1_326
        or scan["projected_survivor_rank_sha256"] != EXPECTED_PROJECTED_RANK_SHA256
        or int(scan["all_dependency_survivors"]) != 0
        or scan["all_dependency_survivor_ranks"] != []
        or scan["z1_branch_excluded"] is not True
        or int(linear["equations"]) != 282
        or int(linear["edge_variables"]) != 1_225
        or int(linear["rank_mod_7"]) != 147
        or int(linear["left_dependency_dimension"]) != 135
        or linear["left_null_audit"] is not True
        or scan["catalog_row_histogram_by_direction_mask"]
        != {"1764": 232, "2233": 280}
        or int(scan["blocks"]) != 65_535
    ):
        raise ArithmeticError("p7 positive z1 complete scan changed")

    matching_keys = (
        "all_boundaries",
        "checked_boundaries",
        "z1_boundaries",
        "mean_allocation_count_per_boundary",
        "projected_dependency_count",
        "projected_survivors",
        "projected_survivor_rank_sha256",
        "first_projected_survivor_ranks",
        "all_dependency_survivors",
        "all_dependency_survivor_ranks",
        "catalog_row_histogram_by_direction_mask",
        "linear_system",
    )
    if (
        any(audit[key] != scan[key] for key in matching_keys)
        or int(audit["blocks"]) != 32_768
    ):
        raise ArithmeticError("different-grid p7 z1 rerun changed")

    actual_before = int(previous["actual_boundary_count_after_z0_exclusion"])
    actual_after = actual_before - int(scan["z1_boundaries"])
    projected_before = previous["remaining_projected_undetermined_direction_histogram"]
    if projected_before != {1: 300, 2: 280, 3: 210, 7: 2}:
        raise ArithmeticError("pre-15.715 projected profile envelope changed")
    projected_after = {key: value for key, value in projected_before.items() if key != 1}
    projected_count_after = sum(projected_after.values())
    if actual_after != 129_024 or projected_count_after != 492:
        raise ArithmeticError("post-15.715 positive branch count changed")

    return {
        "proposition": "15.715",
        "p": 7,
        "boundary": "infinity plus seven finite points",
        "product_sign": "positive",
        "undetermined_direction_count_excluded": 1,
        "actual_boundary_count_before": actual_before,
        "z1_boundaries_excluded": int(scan["z1_boundaries"]),
        "mean_allocations_per_boundary": int(scan["mean_allocation_count_per_boundary"]),
        "actual_boundary_count_after_z1_exclusion": actual_after,
        "projected_b_profile_count_before": int(previous["projected_b_profile_count_after"]),
        "projected_b_profiles_excluded_here": int(projected_before[1]),
        "projected_b_profile_count_after": projected_count_after,
        "remaining_projected_undetermined_direction_histogram": projected_after,
        "projected_mod7_boundary_candidates": int(scan["projected_survivors"]),
        "projected_mod7_boundary_candidate_ranks_stored": len(projected_ranks),
        "full_mod7_survivors": int(scan["all_dependency_survivors"]),
        "modulus": 7,
        "projected_dependency_count": int(scan["projected_dependency_count"]),
        "left_dependency_count": int(linear["left_dependency_dimension"]),
        "scan_evidence": str(SCAN.relative_to(ROOT)),
        "different_grid_rerun_evidence": str(AUDIT.relative_to(ROOT)),
        "projected_rank_certificate_evidence": str(PROJECTED_RANKS.relative_to(ROOT)),
        "evidence_sha256": {
            "primary_scan": SCAN_SHA256,
            "different_grid_rerun": AUDIT_SHA256,
            "projected_rank_certificate": PROJECTED_RANKS_SHA256,
        },
        "scan_script": "scripts/p7_infinity7_positive_z1_mod7_gpu.py",
        "input_integrity": input_integrity,
        "projected_rank_validation": projected_rank_validation,
        "different_grid_rerun_validated": True,
        "independent_implementation_validated": False,
        "validation_scope": (
            "complete projected-then-host-exact CUDA scan plus a same-implementation "
            "different-grid rerun; all 1326 stored projected ranks replay the current "
            "CPU projection, but no independent implementation has validated absence "
            "of GPU projection false negatives"
        ),
        "positive_z1_branch_closed": True,
        "positive_p7_infinity_plus_seven_closed": False,
        "negative_p7_infinity_plus_seven_changed": False,
        "top_level_gates_changed": False,
        "proved_by_complete_exact_finite_scan": True,
    }


def main() -> None:
    theorem = p7_positive_infinity_plus_seven_z1_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15715.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.715: positive p7 infinity+7 z=1 boundaries "
        f"{theorem['z1_boundaries_excluded']} -> 0"
    )


if __name__ == "__main__":
    main()
