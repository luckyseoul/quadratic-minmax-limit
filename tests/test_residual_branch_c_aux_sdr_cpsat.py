from scripts.residual_branch_c_aux_sdr_cpsat import (
    build_chart,
    generate_options,
    solve_scale,
)


def test_branch_c_affine_chart_has_exact_paley_type_partition() -> None:
    chart = build_chart(31)
    assert len(chart.hard_coordinates) == 16
    assert len(chart.opposite_coordinates) == 15
    assert set(chart.hard_coordinates).isdisjoint(chart.opposite_coordinates)
    assert set(chart.hard_coordinates) | set(chart.opposite_coordinates) == set(
        range(31)
    )
    assert all(chart.type_by_coordinate[value] == 1 for value in chart.hard_coordinates)
    assert all(
        chart.type_by_coordinate[value] == -1
        for value in chart.opposite_coordinates
    )


def test_constant_center_profile_has_exact_paired_sdr_witness() -> None:
    chart = build_chart(31)
    alphas = (1,) * len(chart.hard_coordinates)
    # c=18 is a deterministic exact witness in the canonical chart.  The
    # solver may choose a different matching, but its returned assignment is
    # independently replayed by the production function before it is exposed.
    result = solve_scale(
        chart,
        alphas,
        c=18,
        seconds=5.0,
        workers=1,
        random_seed=15766,
    )
    assert result["status"] in ("OPTIMAL", "FEASIBLE")
    witness = result["witness"]
    assert witness["exact_assignment_replay"]
    assert witness["hard_auxiliary_count"] == 14
    assert witness["opposite_auxiliary_count"] == 2
    assert len(set(witness["auxiliary_coordinates"])) == 16


def test_every_generated_option_is_locally_clean() -> None:
    chart = build_chart(31)
    alphas = tuple(range(1, 17))
    options = generate_options(chart, alphas, c=7)
    assert options
    for option in options:
        assert option.first_target < option.second_target
        assert option.first_auxiliary != option.second_auxiliary
        assert option.first_auxiliary not in (
            chart.hard_coordinates[option.first_target],
            chart.hard_coordinates[option.second_target],
        )
        assert option.second_auxiliary not in (
            chart.hard_coordinates[option.first_target],
            chart.hard_coordinates[option.second_target],
        )
