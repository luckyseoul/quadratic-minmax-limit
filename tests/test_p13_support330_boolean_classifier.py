import json
from copy import deepcopy

import pytest

from scripts.p13_support330_boolean_classifier import (
    ANCHORED_NOGOOD_COUNT,
    CATALOG_SIZE,
    DOMAIN_SIZE,
    IDENTITY_COUNT,
    SUPPORT_SIZE,
    all_shard_assignments,
    build_classifier_model,
    build_only_payload,
    catalog_arithmetic_certificate,
    merge_shard_payloads,
    positive_worker_count,
    shard_partition,
    support330_candidate_catalog,
)


def test_support330_candidate_catalog_arithmetic_and_anchor_filter():
    supports, forms = support330_candidate_catalog()
    row = catalog_arithmetic_certificate()

    assert len(supports) == len(set(supports)) == CATALOG_SIZE == 364
    assert {len(support) for support in supports} == {SUPPORT_SIZE}
    assert row["point_count"] == DOMAIN_SIZE == 1716
    assert row["pair_coordinate_count"] == 78
    assert row["third_difference_identity_count"] == IDENTITY_COUNT == 1638
    assert row["identity_nullspace_is_exact_degree_at_most_two_space"] is True
    assert row["support_density"] == "330/1716=5/26"
    assert row["family_counts"] == {
        "omitted_pair": 78,
        "all_equal_triple": 286,
    }
    assert row["anchored_candidate_count"] == ANCHORED_NOGOOD_COUNT == 70
    assert row["anchored_family_counts"] == {
        "omitted_pair": 15,
        "all_equal_triple": 55,
    }
    assert row["b2_combined_target_offsets"] == {
        "omitted_pair": 3,
        "all_equal_triple": 5,
    }
    assert all(row["lift_target_identities_verified"].values())
    assert row["every_candidate_satisfies_all_1638_identities"] is True
    assert row["candidate_catalog_verified"] is True
    assert row["catalog_exhaustiveness_claimed"] is False
    symmetry = row["S13_symmetry_certificate"]
    assert symmetry["generator_count"] == 12
    assert symmetry["catalog_images_checked"] == 12 * CATALOG_SIZE
    assert symmetry["catalog_closed_by_generator"] == [True] * 12
    assert symmetry["catalog_closed_under_generated_S13"] is True
    assert symmetry["anchor_orbit_size"] == DOMAIN_SIZE
    assert symmetry["anchor_orbit_is_all_J13_7"] is True
    assert symmetry["anchor_reduction_verified"] is True
    assert len(row["candidate_catalog_sha256"]) == 64
    assert len(row["anchored_candidate_catalog_sha256"]) == 64
    assert {form["family"] for form in forms} == {
        "omitted_pair",
        "all_equal_triple",
    }


def test_full_model_has_exact_constraints_and_no_theorem_before_solve():
    model, _values, metadata = build_classifier_model()
    assert model.Validate() == ""
    assert metadata["boolean_variable_count"] == DOMAIN_SIZE
    assert metadata["constraint_count"] == 1710
    assert metadata["third_difference_equality_count"] == IDENTITY_COUNT
    assert metadata["support_equation_count"] == 1
    assert metadata["anchor_equation_count"] == 1
    assert metadata["anchored_nogood_count"] == ANCHORED_NOGOOD_COUNT
    assert metadata["shard_fixing_count"] == 0
    assert metadata["partition"] == {
        "scheme": "binary_prefix_of_anchored_boolean_assignment",
        "shard_index": 0,
        "shard_count": 1,
        "prefix_width": 0,
        "fixed_point_indices": [],
        "fixed_values": [],
        "all_70_nogoods_retained_in_every_shard": True,
        "nogoods_are_not_sharded": True,
        "partition_is_disjoint_and_exhaustive": True,
    }
    assert len(metadata["model_textproto_sha256"]) == 64

    payload = build_only_payload(metadata)
    assert payload["result_status"] == "MODEL_BUILT_NOT_SOLVED"
    assert payload["classification"]["full_catalog_exhaustive"] is False
    assert payload["classification"]["incomplete"] is True
    json.dumps(payload)


def test_binary_prefix_shards_are_exact_and_never_split_nogoods():
    assert all_shard_assignments(4) == ((0, 0), (1, 0), (0, 1), (1, 1))
    model, _values, metadata = build_classifier_model(3, 4)
    partition = metadata["partition"]
    assert model.Validate() == ""
    assert metadata["constraint_count"] == 1712
    assert metadata["anchored_nogood_count"] == 70
    assert metadata["shard_fixing_count"] == 2
    assert partition["fixed_point_indices"] == [1, 2]
    assert partition["fixed_values"] == [1, 1]
    assert partition["all_70_nogoods_retained_in_every_shard"] is True
    assert partition["nogoods_are_not_sharded"] is True
    assert partition["partition_is_disjoint_and_exhaustive"] is True

    with pytest.raises(ValueError, match="power of two"):
        shard_partition(0, 3)
    with pytest.raises(ValueError, match="0 <= shard-index"):
        shard_partition(4, 4)

    assert positive_worker_count("1") == 1
    assert positive_worker_count(32) == 32
    with pytest.raises(Exception, match="workers must be positive"):
        positive_worker_count("0")


def test_shard_merger_requires_complete_exact_cover():
    rows = []
    for index in range(2):
        _model, _values, metadata = build_classifier_model(index, 2)
        payload = build_only_payload(metadata)
        payload["solver"] = {"status": "INFEASIBLE"}
        payload["witness"] = None
        rows.append(payload)

    merged = merge_shard_payloads(rows)
    assert merged["result_status"] == "COMPLETE_EXACT_SHARD_COVER_INFEASIBILITY"
    assert merged["partition"]["covered_shard_indices"] == [0, 1]
    assert merged["partition"]["all_70_nogoods_retained_in_every_shard"] is True
    assert merged["classification"]["full_catalog_exhaustive"] is True

    with pytest.raises(ValueError, match="complete cover"):
        merge_shard_payloads(rows[:1])
    duplicate = [rows[0], deepcopy(rows[0])]
    with pytest.raises(ValueError, match="duplicate shard"):
        merge_shard_payloads(duplicate)

    forged_partition = deepcopy(rows)
    forged_partition[1]["model"]["partition"]["fixed_values"] = [0]
    with pytest.raises(ValueError, match="exact binary-prefix shard"):
        merge_shard_payloads(forged_partition)
