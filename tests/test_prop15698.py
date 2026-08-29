from e1_gmin_m4_prop15698 import (
    p19_allb2_boundary_unsat_certificate,
    p19_allb2_profile_exclusion,
)


def test_two_complete_unsat_runs_match_the_exact_boundary_model():
    row = p19_allb2_boundary_unsat_certificate()
    assert row["model"]["all_completed_unsat"] is True
    assert [run["threads"] for run in row["model"]["raw_runs"]] == [8, 16]
    assert row["sign_transfer"]["both_c_H_signs_excluded"] is True
    assert row["profile_excluded"] is True


def test_all_p19_slack_twenty_profiles_are_closed():
    row = p19_allb2_profile_exclusion()
    assert row["p19_profiles_before"] == 4
    assert row["p19_profiles_after"] == 3
    assert row["remaining_slack_histogram"] == {24: 1, 28: 1, 32: 1}
    assert row["all_p19_slack_twenty_profiles_closed"] is True
    assert row["p19_second_all_finite_endpoint_closed"] is False
    assert row["L_status"] == "OPEN"
