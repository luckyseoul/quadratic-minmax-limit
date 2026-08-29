from e1_gmin_m4_prop15696 import (
    AGGREGATE_ROWS,
    p19_b16_aggregate_degree_certificate,
    p19_b16_floor_and_kernel_certificate,
    p19_b16_solver_shard_certificate,
    p19_slack_twenty_b16_exclusion,
)


def test_b16_floor_equality_has_rank_169_and_two_integral_forms():
    row = p19_b16_floor_and_kernel_certificate()
    assert row["type_floor_sum"] == row["type_budget"] == 200
    assert row["forced_pointwise_layers"] == {7: 0, 8: 1, 10: 1}
    assert row["rank_witness_layer_histogram"] == {7: 16, 8: 33, 10: 120}
    assert row["rank_mod_two"] == row["therefore_rational_rank"] == 169
    assert row["explicit_rational_kernel_dimension"] == 2
    assert row["normal_forms"] == {
        "022": "A=1-z+u+v-2uv",
        "400": "A=1+3z-u-v-2zu-2zv+2uv",
    }
    assert row["canonical_pair_target_sums"] == {"b2": 1, "022": -19, "400": -19}


def test_aggregate_capacity_recomputes_exactly_ten_shards():
    row = p19_b16_aggregate_degree_certificate()
    assert tuple(row["admissible_rows"]) == AGGREGATE_ROWS
    assert row["admissible_infinity_degrees"] == [
        2, 8, 10, 12, 18, 20, 28, 30, 38, 48
    ]


def test_all_exact_edge_lift_shards_are_infeasible():
    row = p19_b16_solver_shard_certificate()
    assert row["normal_form_orbits"] == ["022", "400"]
    assert row["shards_per_orbit"] == 10
    assert row["shard_count"] == 20
    assert row["raw_shard_count"] == 22
    assert row["split_logical_shard"] == {
        "shape": "022",
        "infinity_degree": 28,
        "exhaustive_phase_zero_elevated_roles": [0, 2, 16],
    }
    assert row["all_statuses"] == "INFEASIBLE"
    assert all(shard["solver_status"] == "INFEASIBLE" for shard in row["shards"])
    assert "componentwise subtraction" in row["finite_field_sign_convention"]
    assert row["supersedes_original_raw_shards"] is True
    assert row["both_c_h_signs_excluded"] is True


def test_b16_profile_is_excluded_but_p19_endpoint_remains_open():
    row = p19_slack_twenty_b16_exclusion()
    assert row["p19_profiles_before"] == 5
    assert row["p19_profiles_after"] == 4
    assert row["remaining_slack_histogram"] == {20: 1, 24: 1, 28: 1, 32: 1}
    assert row["p19_second_all_finite_endpoint_closed"] is False
    assert row["L_status"] == "OPEN"
