from e1_gmin_m4_prop15682 import (
    p31_complete_arc_classification,
    p31_geometric_exclusion,
    p31_residue_zero_profiles,
    theorem_record,
)


def test_exact_p31_profiles_have_only_arc_or_one_triple_slack():
    row = p31_residue_zero_profiles()
    assert len(row["profiles"]) == 14
    assert row["near_arc_profile_count"] == 3
    assert row["arc_profile_count"] == 11
    assert row["arc_minimum_undetermined_directions"] == 3
    assert {item["pair_slack"] for item in row["profiles"]} == {0, 4}
    near = [item for item in row["profiles"] if item["pair_slack"] == 4]
    arcs = [item for item in row["profiles"] if item["pair_slack"] == 0]
    assert min(item["undetermined_directions"] for item in near) == 5
    assert min(item["undetermined_directions"] for item in arcs) == 3
    for item in row["profiles"]:
        floor_pairs = sum(
            int(secants) * count
            for secants, count in item["global_secant_distribution"].items()
        )
        assert 2 * floor_pairs + item["pair_slack"] == 26 * 25


def test_complete_arc_classification_forces_conic_extensions():
    row = p31_complete_arc_classification()
    assert row["no_complete_arc_sizes"] == list(range(23, 32))
    assert row["largest_nonconic_complete_arc_size"] == 22
    assert row["consequence"] == (
        "every 27- and 28-arc in PG(2,31) is conic-contained"
    )
    assert row["proved_conditional_on_external_classification"] is True


def test_p31_three_infinity_point_contradiction():
    row = p31_geometric_exclusion()
    assert row["arc_case"]["adjoin_two_size"] == 28
    assert row["near_arc_case"]["adjoin_two_size"] == 27
    assert row["excluded"] is True


def test_theorem_record_closes_only_p31_endpoint():
    record = theorem_record()
    assert record["proved"] is True
    theorem = record["theorem"]
    assert theorem["p31_s26_next_all_finite_endpoint"] == "EXCLUDED"
    assert theorem["remaining_smaller_endpoints"] == [17, 19, 23, 41]
    assert theorem["p41_status"] == "ONLY_RESIDUE_ZERO_REMAINS"
    assert theorem["general_residual_ii"] is False
    assert theorem["R1"] is False
    assert theorem["limit_exists"] is False
