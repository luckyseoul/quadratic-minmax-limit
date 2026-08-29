from e1_gmin_m4_prop15711 import p17_residue_zero_uniform_mean_exclusion


def test_uniform_mean_cell_bound_excludes_all_five_residue_zero_profiles():
    row = p17_residue_zero_uniform_mean_exclusion()
    assert row["proved_analytically"] is True
    assert row["uses_solver"] is False
    assert row["profile_count_before"] == 19
    assert row["profiles_excluded_here"] == 5
    assert row["profile_count_after"] == 14
    assert row["remaining_pair_slack_histogram"] == {
        96: 2,
        100: 3,
        104: 3,
        108: 2,
        112: 2,
        116: 1,
        128: 1,
    }
    assert row["remaining_residue_pair_histogram"] == {
        "u0=7,u1=0": 9,
        "u0=8,u1=0": 5,
    }


def test_uniform_mean_infinity_candidates_all_fail_the_fibre_bound():
    row = p17_residue_zero_uniform_mean_exclusion()
    assert row["uniform_directional_mean"] == 18
    assert row["all_finite_edges_forced_to_phase_one"] is True
    assert [
        (
            candidate["infinity_degree"],
            candidate["phase_one_b16_gauge"],
            candidate["phase_one_parallel_count"],
            candidate["cell_upper_bound_on_infinity_degree"],
        )
        for candidate in row["candidate_rows"]
    ] == [
        (6, 1, 7, 2),
        (24, 3, 5, 19),
        (42, 5, 3, 36),
        (60, 7, 1, 53),
    ]
    assert all(
        allocation["avoiding_rigid_b0_forces_every_phase_zero_quotient"] == 1
        and allocation["all_directional_means"] == 18
        for allocation in row["allocation_saturation_rows"]
    )
    assert row["p17_second_all_finite_endpoint_closed"] is False
    assert row["top_level_gates_changed"] is False
