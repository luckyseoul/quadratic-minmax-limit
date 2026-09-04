from e1_gmin_m4_residual_min_degree_barrier import (
    antipodal_matching_aggregate_barrier,
    minimum_degree_chart_bound,
    outside_boundary_even_degree_bound,
    residual_chart_system,
)


def test_minimum_degree_forces_exact_isolated_chart_range() -> None:
    for p in (7, 11, 13, 17, 31):
        n = p * p + 1
        last_isolated = (n - 1) // 2
        assert minimum_degree_chart_bound(p, last_isolated)[
            "isolated_chart_forced"
        ] is True
        threshold = minimum_degree_chart_bound(p, n // 2)
        assert threshold["isolated_chart_forced"] is False
        assert threshold["minimum_degree_upper_bound"] == 1


def test_even_degree_outside_boundary_bound() -> None:
    row = outside_boundary_even_degree_bound(31, 200, 100)
    assert row["outside_even_degree_upper_bound"] == 0
    assert row["isolated_outside_vertex_forced"] is True
    full = outside_boundary_even_degree_bound(31, 481, 962)
    assert full["outside_vertex_count"] == 0
    assert full["outside_even_degree_upper_bound"] is None


def test_antipodal_matching_is_sharp_aggregate_barrier() -> None:
    for p in (7, 11, 13, 17, 19, 23, 31, 43):
        row = antipodal_matching_aggregate_barrier(p)
        assert row["proved"] is True
        assert row["minimum_degree_bound_is_sharp"] is True
        assert row["aggregate_system_passes"] is True
        assert row["all_phase_floors_pass"] is True
        assert row["common_difference_radon_energy"] == row["edge_count"] - 1
        assert row["pointwise_residual_box_passes"] is False
        assert row["is_residual_ii_counterexample"] is False
        assert row["residual_ii_closed"] is False


def test_exact_system_retains_pointwise_box() -> None:
    row = residual_chart_system()
    assert row["proved"] is True
    assert "S_G=2 implies f_e=+1" in row["residual_input"]
    assert "3<=T_H^eps<=Phi-2" in row["signed_shell_form"]
    assert "I+(p+1)P_L-eps_L*T-3p" in row["scaled_directional_mean"]
