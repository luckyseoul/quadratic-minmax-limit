import pytest

from e1_gmin_m4_prop15728 import (
    EVEN_B,
    p31_b2_direction_occupancy,
    p31_endpoint_block_paley_rows,
    p31_endpoint_floor_plus_two_obstruction,
    p31_endpoint_odd_fibre_ledger,
    p31_even_floor_table,
    p31_hard_type_b_split_cases,
    p31_phase_one_residue_ledger,
    p31_phase_zero_residue_parity,
    p31_type_phase_ledger,
    proposition_15728,
)


def test_endpoint_odd_fibre_sum_and_two_exact_type_budgets():
    odd = p31_endpoint_odd_fibre_ledger()
    assert odd["boundary_size"] == 32
    assert odd["endpoint_pair_slack_R"] == 10
    assert odd["possible_b_values"] == list(range(0, 31, 2))
    assert odd["sum_b"] == 72
    assert odd["proved"] is True

    phases = p31_type_phase_ledger()
    assert phases["edge_count"] == 125
    assert phases["parity_exponent"] == 61
    assert phases["directional_parity_sign"] == "-eps_d*c_H"
    assert phases["phase_one_type"] == "eps_d=c_H"
    assert phases["phase_zero_type"] == "eps_d=-c_H"
    assert phases["directions_per_type"] == 16
    assert phases["budget_per_type"] == 512
    assert phases["both_type_budgets_are_exact"] is True
    assert phases["proved"] is True


def test_p31_even_floor_table_has_only_two_phase_one_low_cells():
    row = p31_even_floor_table()
    assert tuple(row["even_b_values"]) == EVEN_B
    assert row["phase_zero_floors"] == {
        0: 0,
        2: 32,
        4: 56,
        6: 62,
        8: 62,
        10: 62,
        12: 62,
        14: 62,
        16: 62,
        18: 62,
        20: 62,
        22: 62,
        24: 62,
        26: 62,
        28: 56,
        30: 32,
    }
    assert row["phase_one_floors"] == {
        b: (30 if b in (2, 30) else 62) for b in EVEN_B
    }
    assert row["proved"] is True


def test_phase_one_common_residue_is_uniquely_thirty():
    row = p31_phase_one_residue_ledger()
    assert row["unique_u"] == 15
    assert row["common_residue"] == 30
    assert row["mean_multiset"] == {30: 15, 62: 1}
    assert row["baseline_direction_count"] == 15
    assert row["baseline_b_values"] == [2, 30]
    assert [item["u"] for item in row["residue_rows"] if item["feasible"]] == [15]
    assert row["floor_plus_two_cells"]["b=2"]["floor_plus_two_forbidden"] is True
    assert row["floor_plus_two_cells"]["b=30"]["floor_plus_two_forbidden"] is True
    assert row["proved"] is True


def test_endpoint_floor_plus_two_cells_force_an_impossible_tiny_lift():
    for b in (2, 30):
        row = p31_endpoint_floor_plus_two_obstruction(b)
        assert row["baseline_scaled_mean_2p_E_q0"] == 30
        assert row["proposed_scaled_mean_2p_E_A"] == 32
        assert row["lift_is_nonnegative_integral_quadratic"] is True
        assert row["induced_scaled_lift_mass_4p_E_C"] == 2
        assert row["prop_15_688_scaled_lift_floor"] == 28
        assert row["floor_plus_two_forbidden"] is True
        assert row["proved"] is True

    with pytest.raises(ValueError):
        p31_endpoint_floor_plus_two_obstruction(4)


def test_phase_zero_half_residue_is_even_because_infinity_degree_is_even():
    row = p31_phase_zero_residue_parity()
    assert row["type_residue_sum"] == "30+2*u0=2*I+6 (mod 32)"
    assert row["halved_congruence"] == "u0=I+4 (mod 16)"
    assert row["infinity_degree_I_even"] is True
    assert row["possible_phase_zero_u0"] == [0, 2, 4, 6, 8, 10, 12, 14]
    assert row["possible_phase_zero_residues"] == [0, 4, 8, 12, 16, 20, 24, 28]
    assert row["proved"] is True


def test_global_b_sum_leaves_at_least_fourteen_hard_type_b2_directions():
    row = p31_hard_type_b_split_cases()
    assert row["global_sum_b"] == 72
    assert row["necessary_case_count"] == 24
    assert row["at_most_one_hard_type_b30_direction"] is True
    assert row["minimum_hard_type_b2_directions"] == 14
    assert all(item["hard_type_sum_b"] + item["other_type_sum_b"] == 72 for item in row["cases"])
    assert all(item["other_type_sum_b"] >= 0 for item in row["cases"])
    assert row["proved"] is True


def test_b2_direction_fibre_accounting_is_exact():
    nonrich = p31_b2_direction_occupancy(0, 0)
    assert nonrich["occupancy_counts"] == {0: 14, 1: 2, 2: 15, 3: 0, 4: 0}
    assert nonrich["point_count"] == 32
    assert nonrich["fibre_count"] == 31
    assert nonrich["pair_count"] == 15
    assert nonrich["nonrich_profile"] is True
    assert nonrich["proved"] is True

    rich = p31_b2_direction_occupancy(2, 3)
    assert rich["occupancy_counts"] == {0: 19, 1: 0, 2: 7, 3: 2, 4: 3}
    assert rich["direction_slack"] == 8
    assert rich["pair_count"] == 31
    assert rich["proved"] is True

    for invalid in ((-1, 0), (0, -1), (3, 0), (2, 7)):
        with pytest.raises(ValueError):
            p31_b2_direction_occupancy(*invalid)


def test_each_block_row_forces_four_plus_y_nonrich_pairing_directions():
    rows = p31_endpoint_block_paley_rows()
    assert len(rows) == 6
    for y, row in enumerate(rows):
        assert row["four_secants_y"] == y
        assert row["trisecants_x"] == 10 - 2 * y
        assert row["rich_line_count"] == 10 - y
        assert row["minimum_hard_type_b2_directions"] == 14
        assert row["minimum_nonrich_hard_type_b2_directions"] == 4 + y
        assert row["nonrich_b2_occupancy"] == {0: 14, 1: 2, 2: 15, 3: 0, 4: 0}
        assert row["proved"] is True


def test_proposition_package_is_a_necessary_normal_form_not_a_close():
    row = proposition_15728()
    assert row["prop"] == "15.728"
    assert row["result_status"] == "proved necessary normal form"
    assert row["hypotheses"] == {
        "residual_affine_separator": True,
        "p": 31,
        "edge_count": 125,
        "all_finite_boundary_size": 32,
        "outside_pair_slack": 10,
        "prop_15_727_disjoint_block_normal_form": True,
    }
    assert row["conclusion"]["mean_multiset"] == {30: 15, 62: 1}
    assert row["conclusion"]["minimum_b2_directions"] == 14
    assert row["finite_configuration_search_used"] is False
    assert row["arc_classification_used"] is False
    assert row["p31_endpoint_excluded"] is False
    assert row["p_plus_one_shell_closed"] is False
    assert row["non_walsh_residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True
