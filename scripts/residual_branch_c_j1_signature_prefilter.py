#!/usr/bin/env python3
"""Replay the bounded p31 j=1 selector-signature prefilter.

This script checks the ``f=1,d=1`` correction pattern on the frozen
10,000-vertex auxiliary-design BFS checkpoint.  It deliberately does not
enumerate new designs.  The checkpoint is a bounded prefix of a component,
so excluding every stored vertex is neither an exhaustive component result
nor a global closure of residual (ii).

For an exact half design let ``G`` be the XOR of its 32 adaptive-kernel
signature words.  The j=1 ledger in this branch has a fixed antipodal edge in
direction ``F``, one unused doubled origin orbit in direction ``N``, and two
nonorigin inversion-orbit cancellations.  Each cancellation word has weight
at most two.  Hence a necessary condition is

    G + e_F + e_N = v_1 + v_2,

and therefore ``min_N wt(G + e_F + e_N) <= 4``.  The scan below evaluates
this exact zero-cost necessary condition from the closed half-signature
formula; it does not construct edges or search centres.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from typing import Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_inversion_antisymmetric_radon import (  # noqa: E402
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (  # noqa: E402
    _parallel_formula,
    paley_direction_sign,
)
from e1_gmin_m4_p31_top_mobius_boundary_parity import (  # noqa: E402
    half_kernel_sigma_formula,
)
from io_atomic import write_json_atomic  # noqa: E402


P = 31
FIXED_DIRECTION_INDEX = 5
DIRECTIONS = tuple(projective_functionals(P))
DIRECTION_SIGNS = tuple(
    paley_direction_sign(P, direction) for direction in DIRECTIONS
)
HARD_TARGET_ORDER = tuple(
    index for index, sign in enumerate(DIRECTION_SIGNS) if sign == 1
)
AUXILIARY_DIRECTION_SET = (
    1,
    2,
    7,
    9,
    10,
    13,
    14,
    15,
    16,
    21,
    22,
    24,
    28,
    29,
    30,
    31,
)
RAW_PARALLEL_PROFILE = (
    15,
    14,
    14,
    15,
    16,
    16,
    16,
    14,
    16,
    14,
    14,
    16,
    16,
    15,
    15,
    14,
    14,
    16,
    16,
    16,
    16,
    14,
    14,
    16,
    14,
    16,
    16,
    16,
    14,
    14,
    14,
    14,
)
P2_EXACT_COVER_COUNTS = tuple(
    (
        RAW_PARALLEL_PROFILE[index]
        - int(index in HARD_TARGET_ORDER)
        - int(index in AUXILIARY_DIRECTION_SET)
    )
    // 2
    for index in range(P + 1)
)

COMPONENT_SCHEMA = "residual_branch_c_auxiliary_component_v1"
OPTION_SCHEMA = "resii_p31_j1_closed_option_catalog_v1"
OUTPUT_SCHEMA = "resii_p31_j1_signature_prefilter_v1"
EXPECTED_COMPONENT_FILE_SHA256 = (
    "196588c21a37c7788565b64c5b2a7dbfcafaedbd864dadf7e51b8b278895ae5b"
)
EXPECTED_OPTION_FILE_SHA256 = (
    "6c2c9dd2ca12d007f865c4499dbe038ef53215688eb43bb3759aef7d39daa599"
)
EXPECTED_OPTION_RECORD_SHA256 = (
    "aa33d94449b90b8770ec819044552f8a575465203fc0ed29d9e69dcef37764ff"
)
EXPECTED_DESIGN_COUNT = 10_000
EXPECTED_FRONTIER_COUNT = 9_018
EXPECTED_TRANSFORMED_DISTANCE_HISTOGRAM = {
    6: 18,
    8: 99,
    10: 624,
    12: 1_549,
    14: 2_821,
    16: 2_482,
    18: 1_669,
    20: 616,
    22: 116,
    24: 6,
}
EXPECTED_WEIGHT_SIX_COMPONENT_IDS = (5497, 8971, 9077)

Option = tuple[int, int, int, int]
DesignKey = tuple[tuple[int, int], ...]
OptionLookup = Mapping[tuple[int, int, int], tuple[int, int]]


def _scaled(direction: tuple[int, int], scalar: int) -> tuple[int, int]:
    return scalar * direction[0] % P, scalar * direction[1] % P


@lru_cache(maxsize=1)
def closed_option_catalog() -> tuple[tuple[tuple[Option, ...], ...], str]:
    """Build all 7,260 allowed options from the exact closed formula.

    An option record is ``(auxiliary_index, scale, P2_mask, g_mask)``.
    The returned hash is over those binary records in hard-target order, not
    over a particular JSON serialization.
    """
    targets: list[tuple[Option, ...]] = []
    digest = hashlib.sha256()
    for target_index in HARD_TARGET_ORDER:
        target = DIRECTIONS[target_index]
        rows = []
        for auxiliary_index in AUXILIARY_DIRECTION_SET:
            if auxiliary_index == target_index:
                continue
            for scale in range(1, P):
                auxiliary = _scaled(DIRECTIONS[auxiliary_index], scale)
                profile = tuple(
                    _parallel_formula(P, target, auxiliary, direction)
                    for direction in DIRECTIONS
                )
                p2_mask = sum(
                    int(value == 2) << index
                    for index, value in enumerate(profile)
                )
                if p2_mask.bit_count() != 14:
                    raise ArithmeticError(
                        "a half option lost its fourteen P=2 directions"
                    )
                g_mask = sum(
                    int(half_kernel_sigma_formula(target, auxiliary, kernel) == -1)
                    << kernel_index
                    for kernel_index, kernel in enumerate(DIRECTIONS)
                )
                row = auxiliary_index, scale, p2_mask, g_mask
                rows.append(row)
                digest.update(
                    bytes((target_index, auxiliary_index, scale))
                    + p2_mask.to_bytes(4, "little")
                    + g_mask.to_bytes(4, "little")
                )
        targets.append(tuple(rows))
    result = tuple(targets)
    expected_counts = (480, 450, 450, 480) + (450,) * 12
    if tuple(map(len, result)) != expected_counts:
        raise ArithmeticError("the closed option counts changed")
    return result, digest.hexdigest()


def _read_pinned_json(
    path: Path, expected_sha256: str, description: str
) -> tuple[dict[str, object], str]:
    raw = path.read_bytes()
    actual_sha256 = hashlib.sha256(raw).hexdigest()
    if actual_sha256 != expected_sha256:
        raise ValueError(
            f"{description} SHA256 mismatch: {actual_sha256} != {expected_sha256}"
        )
    payload = json.loads(raw)
    if not isinstance(payload, dict):
        raise ValueError(f"{description} must be a JSON object")
    return payload, actual_sha256


def validated_option_lookup(
    payload: Mapping[str, object],
) -> tuple[dict[tuple[int, int, int], tuple[int, int]], tuple[int, ...]]:
    """Validate the serialized options against the independent formula."""
    if payload.get("schema") != OPTION_SCHEMA or payload.get("p") != P:
        raise ValueError("the option catalog has the wrong schema or prime")
    if tuple(payload.get("hard", ())) != HARD_TARGET_ORDER:
        raise ValueError("the hard-target order changed")
    if tuple(payload.get("auxiliary_set", ())) != AUXILIARY_DIRECTION_SET:
        raise ValueError("the forced auxiliary set changed")
    if tuple(payload.get("raw", ())) != RAW_PARALLEL_PROFILE:
        raise ValueError("the raw parallel profile changed")
    if tuple(payload.get("n2", ())) != P2_EXACT_COVER_COUNTS:
        raise ValueError("the P=2 exact-cover counts changed")

    closed, record_sha256 = closed_option_catalog()
    if record_sha256 != EXPECTED_OPTION_RECORD_SHA256:
        raise ArithmeticError("the independently rebuilt option-record hash changed")
    if payload.get("option_catalog_sha256") != record_sha256:
        raise ValueError("the option catalog's internal record hash changed")
    serialized = payload.get("options")
    if not isinstance(serialized, list) or len(serialized) != len(closed):
        raise ValueError("the option catalog has the wrong target count")
    parsed = tuple(
        tuple(tuple(int(value) for value in row) for row in rows)
        for rows in serialized
    )
    if parsed != closed:
        raise ValueError("the serialized option table disagrees with the formula")

    lookup: dict[tuple[int, int, int], tuple[int, int]] = {}
    for target_index, rows in zip(HARD_TARGET_ORDER, closed, strict=True):
        for auxiliary_index, scale, p2_mask, g_mask in rows:
            lookup[target_index, auxiliary_index, scale] = p2_mask, g_mask
    return lookup, tuple(map(len, closed))


def parse_design_key(value: object) -> DesignKey:
    if not isinstance(value, list) or len(value) != len(HARD_TARGET_ORDER):
        raise ValueError("a p31 design key needs sixteen entries")
    key = tuple(tuple(int(entry) for entry in row) for row in value)
    if any(len(row) != 2 for row in key):
        raise ValueError("each design-key row must be (auxiliary, scale)")
    if {row[0] for row in key} != set(AUXILIARY_DIRECTION_SET):
        raise ValueError("a design key changed the auxiliary SDR")
    if any(not 1 <= row[1] < P for row in key):
        raise ValueError("a relative scale left F_31^*")
    return key


def aggregate_signature_and_p2_counts(
    key: DesignKey, lookup: OptionLookup
) -> tuple[int, tuple[int, ...]]:
    aggregate = 0
    p2_counts = [0] * (P + 1)
    for target_index, (auxiliary_index, scale) in zip(
        HARD_TARGET_ORDER, key, strict=True
    ):
        try:
            p2_mask, g_mask = lookup[target_index, auxiliary_index, scale]
        except KeyError as error:
            raise ValueError(f"design key contains an unknown option: {error}") from error
        aggregate ^= g_mask
        for direction in range(P + 1):
            p2_counts[direction] += (p2_mask >> direction) & 1
    return aggregate, tuple(p2_counts)


def transformed_distances(
    aggregate_signature: int, fixed_direction_index: int = FIXED_DIRECTION_INDEX
) -> tuple[int, ...]:
    """Return ``wt(G+e_F+e_N)`` for every possible origin direction N."""
    if not 0 <= fixed_direction_index <= P:
        raise ValueError("the fixed direction index is outside P^1(F_31)")
    return tuple(
        (
            aggregate_signature
            ^ (1 << fixed_direction_index)
            ^ (1 << origin_direction)
        ).bit_count()
        for origin_direction in range(P + 1)
    )


def analyze_catalogs(
    component: Mapping[str, object],
    options: Mapping[str, object],
    *,
    component_file_sha256: str,
    option_file_sha256: str,
) -> dict[str, object]:
    """Analyze exactly the pinned checkpoint already loaded in memory."""
    if component.get("schema") != COMPONENT_SCHEMA or component.get("p") != P:
        raise ValueError("the component checkpoint has the wrong schema or prime")
    if tuple(component.get("hard_target_order", ())) != HARD_TARGET_ORDER:
        raise ValueError("the component hard-target order changed")
    if tuple(component.get("auxiliary_direction_set", ())) != (
        AUXILIARY_DIRECTION_SET
    ):
        raise ValueError("the component auxiliary set changed")
    if tuple(component.get("base_raw_parallel_profile", ())) != (
        RAW_PARALLEL_PROFILE
    ):
        raise ValueError("the component raw profile changed")
    if component.get("component_exhausted") is not False:
        raise ValueError("this replay is pinned to the nonexhausted checkpoint")
    if int(component.get("frontier_design_count", -1)) != EXPECTED_FRONTIER_COUNT:
        raise ValueError("the pinned BFS frontier count changed")

    lookup, option_counts = validated_option_lookup(options)
    serialized_keys = component.get("all_design_keys")
    if not isinstance(serialized_keys, list):
        raise ValueError("the component checkpoint has no design-key list")
    if len(serialized_keys) != EXPECTED_DESIGN_COUNT:
        raise ValueError("the pinned design count changed")

    histogram: Counter[int] = Counter()
    weight_six_records = []
    excluded = 0
    signature_stream = hashlib.sha256()
    for component_id, serialized_key in enumerate(serialized_keys):
        key = parse_design_key(serialized_key)
        aggregate, p2_counts = aggregate_signature_and_p2_counts(key, lookup)
        if p2_counts != P2_EXACT_COVER_COUNTS:
            raise ArithmeticError(
                f"component {component_id} failed the exact profile replay"
            )
        signature_stream.update(aggregate.to_bytes(4, "little"))
        distances = transformed_distances(aggregate)
        minimum = min(distances)
        histogram[minimum] += 1
        excluded += int(minimum > 4)
        if aggregate.bit_count() == 6:
            weight_six_records.append(
                {
                    "component_id": component_id,
                    "aggregate_signature_hex": f"{aggregate:08x}",
                    "aggregate_signature_support": [
                        index for index in range(P + 1) if aggregate >> index & 1
                    ],
                    "fixed_direction_in_signature": bool(
                        aggregate >> FIXED_DIRECTION_INDEX & 1
                    ),
                    "minimum_after_fixed_and_origin_double": minimum,
                    "minimizing_origin_directions": [
                        index
                        for index, value in enumerate(distances)
                        if value == minimum
                    ],
                    "design_key": key,
                }
            )

    actual_histogram = dict(sorted(histogram.items()))
    if actual_histogram != EXPECTED_TRANSFORMED_DISTANCE_HISTOGRAM:
        raise ArithmeticError("the transformed-distance histogram changed")
    if tuple(row["component_id"] for row in weight_six_records) != (
        EXPECTED_WEIGHT_SIX_COMPONENT_IDS
    ):
        raise ArithmeticError("the weight-six counterexample list changed")
    if excluded != EXPECTED_DESIGN_COUNT:
        raise ArithmeticError("a checkpoint design passed the necessary prefilter")

    return {
        "schema": OUTPUT_SCHEMA,
        "classification": (
            "exhaustive prefilter replay over this bounded 10,000-vertex BFS "
            "checkpoint only; the BFS component is not exhausted; not global"
        ),
        "p": P,
        "component_catalog_schema": COMPONENT_SCHEMA,
        "component_catalog_file_sha256": component_file_sha256,
        "component_discovered_design_count": EXPECTED_DESIGN_COUNT,
        "component_processed_design_count": int(
            component.get("processed_design_count", -1)
        ),
        "component_frontier_design_count": EXPECTED_FRONTIER_COUNT,
        "component_exhausted": False,
        "option_catalog_schema": OPTION_SCHEMA,
        "option_catalog_file_sha256": option_file_sha256,
        "option_catalog_record_sha256": EXPECTED_OPTION_RECORD_SHA256,
        "option_counts_by_target": option_counts,
        "option_count": sum(option_counts),
        "fixed_direction_index": FIXED_DIRECTION_INDEX,
        "forced_auxiliary_direction_set": AUXILIARY_DIRECTION_SET,
        "raw_parallel_profile": RAW_PARALLEL_PROFILE,
        "p2_exact_cover_counts": P2_EXACT_COVER_COUNTS,
        "necessary_equation": "G + e_F + e_N = v_collision_1 + v_collision_2",
        "one_nonorigin_collision_signature_max_weight": 2,
        "two_collision_xor_max_weight": 4,
        "transformed_distance": "min_N wt(G + e_F + e_N)",
        "minimum_transformed_distance_histogram": actual_histogram,
        "checkpoint_design_count": EXPECTED_DESIGN_COUNT,
        "excluded_checkpoint_design_count": excluded,
        "passing_checkpoint_design_count": EXPECTED_DESIGN_COUNT - excluded,
        "aggregate_signature_stream_sha256": signature_stream.hexdigest(),
        "discarded_conjecture": "every exact-profile design has wt(G) >= 8",
        "discarded_conjecture_counterexample_count": len(weight_six_records),
        "weight_six_aggregate_designs": weight_six_records,
        "bounded_checkpoint_f1_d1_excluded": True,
        "full_bfs_component_exhaustively_enumerated": False,
        "global_f1_d1_branch_closed": False,
        "j1_f3_d0_branch_addressed": False,
        "residual_ii_closed": False,
    }


def replay(component_path: Path, option_path: Path) -> dict[str, object]:
    component, component_sha256 = _read_pinned_json(
        component_path,
        EXPECTED_COMPONENT_FILE_SHA256,
        "component checkpoint",
    )
    options, option_sha256 = _read_pinned_json(
        option_path,
        EXPECTED_OPTION_FILE_SHA256,
        "option catalog",
    )
    return analyze_catalogs(
        component,
        options,
        component_file_sha256=component_sha256,
        option_file_sha256=option_sha256,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("component", type=Path)
    parser.add_argument("options", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = replay(args.component, args.options)
    write_json_atomic(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
