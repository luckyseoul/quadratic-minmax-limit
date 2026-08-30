"""Prop. 15.724: full Miquelian-circle boundary exclusion."""
from __future__ import annotations

import pytest

from e1_gmin_m4_prop15724 import (
    full_circle_lift_contradiction,
    isolated_circle_chart,
    isolated_outside_vertex_ledger,
    phase_one_xnor_normal_form,
    theorem_full_circle_exclusion,
    universal_full_circle_certificate,
    xnor_coefficient_congruence_ledger,
    xnor_lift_certificate,
    zero_infinity_circle_arithmetic,
)


def test_an_outside_isolated_vertex_is_forced():
    for p in (17, 19, 23, 29, 101):
        row = isolated_outside_vertex_ledger(p)
        assert row["maximum_nonisolated_outside"] == 7 * p + 1
        assert row["guaranteed_isolated_outside_vertices"] == p * p - 8 * p - 1
        assert row["isolated_outside_vertex_exists"] is True


def test_isolated_chart_has_exact_aligned_zero_and_two_profile():
    for p in (17, 19, 23, 29, 31, 41):
        row = isolated_circle_chart(p)
        m = (p + 1) // 2
        assert row["transported_infinity_degree_I"] == 0
        assert row["phase_zero_type"] == {"directions": m, "b": 0}
        assert row["phase_one_type"] == {"directions": m, "b": 2}
        assert row["phase_one_floor"] == p - 1
        assert row["type_alignment_exact"] is True


def test_phase_one_xnor_type_has_one_elevation_only():
    for p in (17, 19, 23, 29, 101):
        row = phase_one_xnor_normal_form(p)
        m = (p + 1) // 2
        assert row["unique_residue_u"] == m - 1
        assert row["quotient_multiset"] == {"0": m - 1, "1": 1}
        assert [r["u"] for r in row["residue_rows"] if r["feasible"]] == [m - 1]
        assert row["floor_plus_two_xnor_lift_certificate"]["proved"] is True


def test_b2_xnor_lift_is_pointwise_valid_and_nonzero():
    for p in (17, 19, 23, 101):
        row = xnor_lift_certificate(p, 2)
        assert row["xnor_truth_table"] == {"00": 1, "01": 0, "10": 0, "11": 1}
        assert row["same_parity_as_b2_phase_one"] is True
        assert row["B_has_degree_at_most_two"] is True
        assert row["B_nonzero"] is True
        assert row["forced_4p_E_B"] == 2
        assert row["excluded"] is True


def test_xnor_coefficient_congruence_is_reestablished_in_this_chart():
    for p in (17, 19, 23, 101):
        row = xnor_coefficient_congruence_ledger(p)
        assert "XNOR/XOR" in row["prior_two_sign_targets"]
        assert row["sign_parameter"].startswith("tau in")
        assert row["divisibility"] == "q divides I+P_d-4"
        assert row["applies_to_b2_phase_one_here"] is True
        assert row["proved"] is True


def test_zero_infinity_arithmetic_forces_4_4_3():
    for p in (17, 19, 23, 29, 31, 41, 101):
        row = zero_infinity_circle_arithmetic(p)
        assert row["phase_zero_unique_u"] == 4
        assert row["x_plus_y"] == 7
        assert row["xnor_baseline_count_x"] == 4
        assert row["phase_zero_baseline_count_y"] == 3
        assert row["phase_zero_zero_quotient_directions_at_least"] == 4
        assert row["zero_quotient_scaled_mean"] == 8


def test_forced_mean_eight_violates_the_sharp_lift_floor():
    for p in (17, 19, 23, 29, 31, 41, 101):
        row = full_circle_lift_contradiction(p)
        assert row["forced_direction"]["forced_4p_E_B"] == 8
        assert row["prop_15_688_lower_bound"] == p - 3
        assert row["contradiction_gap"] == p - 11
        assert row["full_circle_excluded"] is True
        assert row["proved"] is True


def test_theorem_closes_only_R_zero_not_the_whole_shell():
    row = theorem_full_circle_exclusion()
    assert row["proved"] is True
    assert row["universal_certificate"]["proved"] is True
    assert row["theorem"]["full_Miquelian_circle_boundary"] == "EXCLUDED"
    assert row["theorem"]["outside_R_zero"] == "EXCLUDED"
    assert row["theorem"]["outside_R_two"] == "EXCLUDED_BY_15.722"
    assert row["theorem"]["outside_R_three"] == "EXCLUDED_BY_15.722"
    assert row["theorem"]["strict_outside_profiles_R_at_least_4"].startswith(
        "OPEN_AFTER"
    )
    assert row["theorem"]["whole_p_plus_one_shell"] == "OPEN"
    assert row["theorem"]["residual_ii"] is False
    assert row["theorem"]["limit_exists"] is False


def test_full_circle_scope_is_not_inferred_from_sample_primes():
    row = universal_full_circle_certificate()
    assert row["scope"] == "every odd prime p>=17"
    assert row["isolated_vertex_base_value_at_17"] > 0
    assert row["residue_base_check"] is True
    assert row["baseline_count_base_check"] is True
    assert row["proved"] is True


def test_full_circle_certificates_reject_odd_composite_moduli():
    for certificate in (
        isolated_outside_vertex_ledger,
        isolated_circle_chart,
        phase_one_xnor_normal_form,
        xnor_lift_certificate,
        xnor_coefficient_congruence_ledger,
        zero_infinity_circle_arithmetic,
        full_circle_lift_contradiction,
    ):
        with pytest.raises(ValueError, match="prime"):
            certificate(21)
