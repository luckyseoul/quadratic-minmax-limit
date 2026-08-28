from fractions import Fraction

from e1_gmin_m4_prop15681 import (
    endpoint_residue_ledger,
    p29_arc_classification_ledger,
    p29_geometric_exclusion,
    p29_residue_zero_profiles,
    paired_cube_integral_quadratic_floor,
    pgl2_group_order,
    pgl2_subset_orbit_audit,
    theorem_record,
)


def test_paired_cube_integral_quadratic_floor():
    expected = {29: 14, 31: 16, 37: 18, 41: 20}
    for p, floor in expected.items():
        row = paired_cube_integral_quadratic_floor(p)
        assert row["rho"] == Fraction(1, p + 1)
        assert row["cube_average_floor"] == Fraction(1, 4)
        assert row["universal_scaled_mass_floor"] == floor
        assert row["universal_mass_floor"] == Fraction(floor, 4 * p)
        assert row["proved"] is True


def test_small_endpoint_positive_residues_are_all_removed():
    expected = {
        29: [0, 2, 3, 4, 5],
        31: [0, 2, 3, 4, 5, 6],
        37: [2, 3, 4, 5],
        41: [0, 2, 3, 4, 5, 6, 7],
    }
    for p, survivors in expected.items():
        row = endpoint_residue_ledger(p)
        assert [item["u0"] for item in row["pair_survivors"]] == survivors
        assert row["positive_residues_all_excluded"] is True
        assert all(item["therefore_b_zero"] for item in row["positive_residue_rows"])
        assert all(item["excluded"] for item in row["positive_residue_rows"])
        assert row["residue_zero_remains"] is (0 in survivors)


def test_p29_pair_slack_leaves_only_arc_and_one_triple_shapes():
    row = p29_residue_zero_profiles()
    assert row["near_arc_profile_count"] == 1
    assert row["arc_profile_count"] == 4
    assert row["arc_minimum_undetermined_directions"] == 4
    near = [item for item in row["profiles"] if item["pair_slack"] == 4]
    assert len(near) == 1
    assert near[0]["global_secant_distribution"] == {
        "0": 6,
        "11": 14,
        "12": 10,
    }
    arcs = [item for item in row["profiles"] if item["pair_slack"] == 0]
    assert {item["undetermined_directions"] for item in arcs} == {4, 5}


def test_pgl2_orbits_exhaust_classified_25_and_26_arcs():
    assert pgl2_group_order(29) == 29 * (29 * 29 - 1)
    four = pgl2_subset_orbit_audit(29, 4)
    five = pgl2_subset_orbit_audit(29, 5)
    assert four["subset_count"] == 27405
    assert four["orbit_count"] == 5
    assert five["subset_count"] == 142506
    assert five["orbit_count"] == 10
    classification = p29_arc_classification_ledger()
    assert classification["classified_projective_arc_classes"] == {25: 10, 26: 5}
    assert classification["all_25_and_26_arcs_conic_contained"] is True


def test_p29_geometric_contradiction_and_scope():
    geometry = p29_geometric_exclusion()
    assert geometry["arc_case"]["adjoin_two_size"] == 26
    assert geometry["near_arc_case"]["adjoin_two_size"] == 25
    assert geometry["excluded"] is True
    record = theorem_record()
    assert record["proved"] is True
    theorem = record["theorem"]
    assert theorem["p29_s24_next_all_finite_endpoint"] == "EXCLUDED"
    assert theorem["remaining_smaller_endpoints"] == [17, 19, 23, 31, 41]
    assert theorem["general_residual_ii"] is False
    assert theorem["R1"] is False
    assert theorem["limit_exists"] is False
