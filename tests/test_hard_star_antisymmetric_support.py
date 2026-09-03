import pytest

from e1_gmin_m4_hard_star_antisymmetric_support import (
    antipodal_pair_total_ledger,
    balanced_ray_antisymmetric_support_barrier,
    equality_case_projective_pencil_lemma,
    exceptional_direction_support_bound,
    hard_star_antisymmetric_target,
    parallel_pair_total_ledger,
    symmetric_half_norm_barrier,
    theorem_record,
    two_zero_direction_pencil_counterexample,
)


def test_exact_real_and_combinatorial_support_bounds_are_separated():
    for p in (31, 43, 47, 59):
        r = (p - 3) // 4
        maximum_active = 2 * r + 2
        row = hard_star_antisymmetric_target(p, maximum_active)
        assert row["proved"]
        assert row["nonzero_cells_per_active_row"] == 2 * (p - 2)
        assert row["target_antisymmetric_norm_squared"] == maximum_active * 2 * (
            p - 2
        )
        assert row["real_bound_below_one_even_when_every_hard_row_is_active"]
        assert row["single_antipodal_orbit_floor"] == p - 2
        assert row["ternary_support_floor"] == 2 * (p - 2)


def test_exceptional_direction_incidence_refines_equality():
    row = exceptional_direction_support_bound(31, 9, 19)
    assert row["proved"]
    assert row["incidence_inequality"] == "A*(p-2)<=A*c-E"
    assert row["single_antipodal_orbit_floor"] == 32
    assert row["ternary_support_floor"] == 64
    assert row["equality_at_c_equals_p_minus_2_forces_E_zero"]


def test_scalar_orbit_and_parallel_ledgers_remain_feasible_at_both_endpoints():
    for p in (31, 43, 47):
        r = (p - 3) // 4
        for t in (2 * r * r - 4 * r - 2, 4 * r * r - 2 * r - 5):
            orbit = antipodal_pair_total_ledger(p, t)
            parallel = parallel_pair_total_ledger(p, t)
            assert orbit["proved"]
            assert parallel["proved"]
            assert not orbit["total_and_fixed_capacity_contradiction"]
            assert not parallel["scalar_parallel_pair_totals_contradict"]
            assert not parallel["common_edge_orbits_or_midpoints_constructed"]


def test_old_real_norm_bound_has_strictly_positive_slack():
    for p in (31, 43, 47, 59):
        row = symmetric_half_norm_barrier(p)
        assert row["proved"]
        assert row["twice_p_squared_times_margin"] > 0
        assert row["every_shifted_coefficient_positive"]
        assert not row["real_norm_and_symmetric_half_contradiction"]


def test_nine_active_equality_rows_force_one_projective_pencil():
    for p in (31, 43, 47, 59):
        row = equality_case_projective_pencil_lemma(p, 9)
        assert row["proved"]
        assert row["pair_resultant_degree"] == 8
        assert row["projective_root_count_threshold"] == 9
        assert row["conclusion_pairwise_projective_endpoint_intersection"]
        assert row["nonpencil_triangle_orbit_cap"] == 6
        assert row["single_orbit_count"] == p - 2
        assert row["conclusion_common_projective_source_vertex"]
        assert row["conclusion_active_phase_coherence"] == "j_L^2=L(P)^2"


def test_two_zero_rows_do_not_by_themselves_kill_the_pencil():
    for p in (31, 43, 47, 59):
        row = two_zero_direction_pencil_counterexample(p)
        assert row["proved"]
        assert row["single_nonfixed_edge_orbits"] == p - 2
        assert row["x_nonzero_fibre_sums_zero"]
        assert row["y_nonzero_fibre_sums_zero"]
        assert row["x_zero_fibre_sum"] == 1
        assert row["y_zero_fibre_sum"] == 1
        assert row["x_antisymmetric_edge_Radon_zero"]
        assert row["y_antisymmetric_edge_Radon_zero"]
        assert not row["active_hard_row_bijections_satisfied"]
        assert not row["residual_ii_counterexample"]


def test_balanced_barrier_and_theorem_record_keep_global_scope_open():
    row = balanced_ray_antisymmetric_support_barrier(31)
    assert row["proved"]
    assert row["single_orbit_floor"] == 29
    assert not row["total_edge_count_contradiction"]
    assert not row["parallel_pair_total_contradiction"]
    assert not row["real_norm_contradiction"]
    assert row["antisymmetric_ternary_lift_now_proved"]
    assert "Mobius" in row["former_antisymmetric_gate_superseded_by"]
    assert not row["one_common_simple_graph_constructed"]
    assert not row["residual_ii_closed"]

    record = theorem_record()
    assert record["proved"][
        "equality_with_at_least_nine_active_rows_forces_projective_pencil"
    ]
    assert not record["proved"]["two_zero_rows_alone_exclude_that_pencil"]
    assert record["proved"][
        "antisymmetric_ternary_lift_proved_by_subsequent_mobius_trade"
    ]
    assert not record["proved"]["coupled_symmetric_half_proved"]
    assert not record["proved"]["residual_ii_closed"]
    assert record["L_status"] == "OPEN"


def test_parameter_guards():
    with pytest.raises(ValueError):
        hard_star_antisymmetric_target(29, 1)
    with pytest.raises(ValueError):
        hard_star_antisymmetric_target(31, 17)
    with pytest.raises(ValueError):
        exceptional_direction_support_bound(31, 0, 0)
    with pytest.raises(ValueError):
        equality_case_projective_pencil_lemma(31, 8)
