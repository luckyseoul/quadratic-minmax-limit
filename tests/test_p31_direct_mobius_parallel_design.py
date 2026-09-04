from e1_gmin_m4_p31_direct_mobius_parallel_design import (
    centered_physical_graph,
    centered_physical_parallel_design_certificate,
    direct_parallel_design_certificate,
    transverse_compact_l1_diagnostic,
)


def test_direct_half_profiles_hit_the_exact_top_parallel_ledger() -> None:
    row = direct_parallel_design_certificate()
    assert row["proved"]
    assert row["profiles_replayed_from_physical_edges"]
    assert row["auxiliary_hard_count"] == 14
    assert row["auxiliary_opposite_count"] == 2
    assert len(set(row["auxiliary_direction_indices"])) == 16
    assert row["fixed_direction_index"] == 5
    assert row["cancellation_direction_index"] == 5
    assert row["final_hard_multiset"] == (14,) * 14 + (15,) * 2
    assert row["final_opposite_multiset"] == (15,) * 3 + (16,) * 13


def test_public_graph_api_replays_hash_and_canonical_centers() -> None:
    row = centered_physical_graph()
    assert row["proved"]
    assert row["edge_count"] == 479
    assert len(row["edges"]) == len(set(row["edges"])) == 479
    assert row["graph_sha256"] == (
        "c0b32bdf228401ba5ffe68be543b9e6fddb31f86594ff953e1d290a6faeeae0d"
    )
    assert len(row["hard_target_centers"]) == 16
    assert len(
        {record["target_direction_index"] for record in row["hard_target_centers"]}
    ) == 16


def test_frozen_graph_honestly_fails_the_transverse_atom_budget() -> None:
    row = transverse_compact_l1_diagnostic()
    assert row["proved"]
    assert row["violating_row_count"] == 32
    assert row["minimum_l1_excess"] == 122
    assert row["maximum_l1_excess"] == 194
    assert row["frozen_graph_fails_necessary_compact_l1_budget"]


def test_center_scaling_realizes_one_clean_nonorigin_cancellation() -> None:
    row = centered_physical_parallel_design_certificate()
    assert row["proved"]
    assert row["unique_pair_intersection_count"] == 1
    assert row["pair_intersections"][0]["halves"] == (2, 13)
    assert row["cancelled_orbit"] == ((2, 25), (29, 1))
    assert row["cancellation_is_nonorigin"]
    assert row["cancellation_direction_index"] == 5
    assert row["fixed_direction_index"] == 5
    assert row["surviving_nonfixed_orbit_count"] == 478
    assert row["graph_edge_count"] == 479
    assert row["origin_edge_count"] == 16
    assert row["same_sign_overlap_count"] == 0
    assert row["triple_overlap_count"] == 0
    assert row["exact_physical_edge_radon_replay"]
    assert row["final_hard_multiset"] == (14,) * 14 + (15,) * 2
    assert row["final_opposite_multiset"] == (15,) * 3 + (16,) * 13
