from src.e1_gmin_m4_prop15705 import p17_slack_sixteen_orbit_exclusion


def test_prop15705_orbit_census_and_accounting():
    theorem = p17_slack_sixteen_orbit_exclusion()
    census = theorem["orbit_extension_census"]
    assert census["raw_four_point_extensions_with_core_secant_charge_at_most_four"] == 97122
    assert census["occupancy_valid_extensions"] == 47
    assert census["phase_labelled_target_hits"] == 0
    assert theorem["profile_count_before"] == 654
    assert theorem["profiles_excluded_here"] == 13
    assert theorem["profile_count_after"] == 641
    assert theorem["remaining_slack_sixteen_profiles"] == 0
