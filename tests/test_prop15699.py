from e1_gmin_m4_prop15699 import (
    p19_endpoint_boundary_unsat_certificate,
    p19_second_endpoint_exclusion,
)


def test_all_three_endpoint_profiles_have_completed_unsat_runs():
    row = p19_endpoint_boundary_unsat_certificate()
    assert row["excluded_slacks"] == [24, 28, 32]
    assert row["model"]["all_completed_unsat"] is True
    assert row["model"]["uses_edge_lift_variables"] is False
    assert row["model"]["uses_floor_relaxation"] is False
    assert row["sign_transfer"]["both_c_H_signs_excluded"] is True


def test_p19_second_endpoint_is_closed():
    row = p19_second_endpoint_exclusion()
    assert row["p19_profiles_before"] == 3
    assert row["p19_profiles_after"] == 0
    assert row["remaining_slack_histogram"] == {}
    assert row["p19_second_all_finite_endpoint_closed"] is True
    assert row["closes_residual_ii"] is False
    assert row["L_status"] == "OPEN"
