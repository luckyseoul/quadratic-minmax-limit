from src.e1_gmin_m4_prop15705 import p17_slack_sixteen_orbit_exclusion


def test_prop15705_orbit_census_is_honestly_limited_to_historical_targets():
    theorem = p17_slack_sixteen_orbit_exclusion()
    census = theorem["orbit_extension_census"]
    assert census["raw_four_point_extensions_with_core_secant_charge_at_most_four"] == 97122
    assert census["occupancy_valid_extensions"] == 47
    assert census["phase_labelled_target_hits"] == 0
    assert theorem["profile_count_before"] == 1228
    assert theorem["profiles_excluded_here"] == 13
    assert theorem["profile_count_after"] == 1215
    assert theorem["historical_orbiter_target_profile_count"] == 13
    assert theorem["corrected_zero_direction_slack_sixteen_profile_count"] == 87
    assert theorem["orbiter_uncovered_slack_sixteen_profile_count"] == 74
    assert theorem["remaining_slack_sixteen_profiles"] == 74
    assert theorem["slack_sixteen_block_closed_here"] is False
    assert theorem["historical_claim_of_final_slack_sixteen_closure"] is False
    assert theorem["remaining_slack_sixteen_status_here"] == "OPEN"
    assert theorem["proof_status"] == "PARTIAL"
    assert len(theorem["remaining_profile_indices"]) == 1215
