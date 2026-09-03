import copy

import pytest

import e1_gmin_m4_conic_odd_radon as conic
from e1_gmin_m4_equianharmonic_component_packing import (
    equianharmonic_component_packing_certificate,
    equianharmonic_exact_fiber_threshold,
    p3_odd_radon_centrality_component_upgrade,
    p43_equianharmonic_threshold_witness_certificate,
)


@pytest.fixture(scope="module")
def packing():
    return equianharmonic_component_packing_certificate()


def test_positive_excess_classification_and_uniformity(packing):
    assert packing["proved"] is True
    assert packing["component_identity"] == (
        "delta=K+2*AE-2+2*cycle_rank+caps"
    )
    assert packing["component_excess_identity"] == (
        "K-2*delta=4-K-4*AE-4*cycle_rank-2*caps"
    )
    tuples = packing["positive_excess_tuple_audit"]
    assert tuples["proved"] is True
    assert tuples["positive_integer_tuples_K_AE_cycle_rank_caps_excess"] == [
        (1, 0, 0, 0, 3),
        (1, 0, 0, 1, 1),
        (2, 0, 0, 0, 2),
        (3, 0, 0, 0, 1),
    ]
    assert packing["score_three_compact_no_go"]["proved"] is True
    assert (
        packing["score_three_compact_no_go"][
            "distinct_label_score_three_compact_exists"
        ]
        is False
    )
    assert packing["constant_target_signs_covered"] == [-1, 1]
    assert packing["orbit_support_collision_exception_primes"] == [
        2,
        3,
        7,
        13,
    ]
    assert packing[
        "discarded_branch_relation_change_exception_primes"
    ] == [2, 3, 5, 7, 13, 19]
    assert packing["all_exceptional_characteristics_below_31"] is True


def test_all_r1_reverse_positions_and_positive_outside_high(packing):
    r1 = packing["r1_exhaustive_audit"]
    assert r1["formal_branch_count"] == 192
    assert r1["attempts_by_reverse_position"] == {
        "positive": 64,
        "negative_ac": 64,
        "negative_bc": 64,
    }
    assert r1["consistent_branches_by_reverse_position"] == {
        "positive": 8,
        "negative_ac": 0,
        "negative_bc": 0,
    }
    assert r1["valid_branches_by_reverse_position"] == {
        "positive": 8,
        "negative_ac": 0,
        "negative_bc": 0,
    }
    assert r1["inconsistent_relation_exception_primes"] == [2, 3, 7]
    assert packing["r1_triangle_flips_equal_exhaustive_supports"] is True
    assert packing["r1_supports_equal_o2_supports"] is True

    high = packing["positive_outside_high_audit"]
    assert high["formal_multiplier_pairs"] == 4
    assert high["equal_multiplier_repeated_label_families"] == 2
    assert high["unequal_multiplier_reverse_target_families"] == 2
    assert high["self_antipodal_ordered_candidates"] == 2
    assert high["distinct_self_antipodal_supports"] == 1
    assert high["support"] == packing["cap"]["support"]


def test_symbolic_support_intersections_and_weighted_packing(packing):
    assert packing["o2"]["branch_counts"] == {
        "rank_two": 176,
        "rank_one": 8,
        "inconsistent": 72,
    }
    assert packing["o2"]["valid_rank_two_assignments"] == 48
    assert packing["hh"]["valid_rank_two_assignments"] == 8
    assert packing["o2_pairwise_intersections"] == [2] * 6
    assert packing["hh_o2_intersections"] == [3] * 4
    assert packing["hh_cap_intersection"] == 0
    assert packing["cap_o2_intersections"] == [0] * 4
    assert packing["weighted_disjoint_packing_maximum"] == 3
    assert ["HH", "cap"] in packing[
        "weighted_disjoint_packing_maximizers"
    ]


def test_exact_fiber_threshold_and_sign_symmetry():
    below31 = equianharmonic_exact_fiber_threshold(31, 6)
    boundary31 = equianharmonic_exact_fiber_threshold(31, 7)
    below43 = equianharmonic_exact_fiber_threshold(43, 8)
    boundary43 = equianharmonic_exact_fiber_threshold(43, 9)
    assert below31["minimum_compact_atom_count"] == 7
    assert below31["exact_fiber_excluded"] is True
    assert boundary31["exact_fiber_excluded"] is False
    assert below43["minimum_compact_atom_count"] == 9
    assert below43["exact_fiber_excluded"] is True
    assert boundary43["exact_fiber_excluded"] is False
    assert boundary43["constant_target_signs_covered"] == [-1, 1]
    assert boundary43["threshold_is_necessary_not_sufficient"] is True
    assert boundary43["residual_ii_closed"] is False
    for bad in ((29, 1), (31, -1), (31, 8), (43, True)):
        with pytest.raises(ValueError):
            equianharmonic_exact_fiber_threshold(*bad)


def test_p43_threshold_witness_replays_every_edge_and_odd_channel():
    witness31 = conic.p31_equianharmonic_witness_certificate()
    threshold31 = equianharmonic_exact_fiber_threshold(31, witness31["b"])
    assert witness31["proved"] is True
    assert witness31["edge_orbit_replay_exact"] is True
    assert witness31["all_odd_channels_zero"] is True
    assert witness31["b"] == threshold31["minimum_compact_atom_count"] == 7

    witness = p43_equianharmonic_threshold_witness_certificate()
    assert witness["proved"] is True
    assert (witness["p"], witness["b"], witness["k"], witness["q"]) == (
        43,
        9,
        13,
        36,
    )
    assert len(witness["ae_atoms"]) == 9
    assert len(witness["compact_atoms"]) == 9
    assert witness["target_support"] == 41
    assert witness["target_l1"] == 41
    assert witness["edge_orbit_replay_exact"] is True
    assert witness["odd_channel_count"] == 210
    assert witness["all_odd_channels_zero"] is True
    assert witness["degree_six"] == [37, 19, 8]
    assert witness["degree_eight"] == [18, 17, 10, 32]
    assert witness["degree_six_and_eight_both_zero"] is False
    assert witness["attains_equianharmonic_compact_threshold"] is True
    assert witness["residual_ii_closed"] is False


def test_dependency_gated_zero_odd_centrality_and_scope(monkeypatch):
    p31 = p3_odd_radon_centrality_component_upgrade(31, 6)
    p43 = p3_odd_radon_centrality_component_upgrade(43, 8)
    p47 = p3_odd_radon_centrality_component_upgrade(47, 11)
    assert p31["proved"] is True
    assert p31["compact_count_hypothesis"] == "3*b<=2*r+4"
    assert p43["proved"] is True
    assert p47["proved"] is True
    assert p47["compact_count_hypothesis"] == "0<=b<=r"
    assert p47["equianharmonic_branch_exists"] is False
    assert p47["assumes_zero_odd_global_forms"] is True
    assert p47["nonzero_global_forms_ruled_out"] is False
    assert p47["global_common_edge_lift_constructed"] is False
    assert p47["Boolean_lift_constructed"] is False
    assert p47["residual_ii_closed"] is False
    with pytest.raises(ValueError):
        p3_odd_radon_centrality_component_upgrade(31, 7)

    original = conic.theorem_record

    def broken_theorem_record():
        record = copy.deepcopy(original())
        record["proved"][
            "constant_branch_forces_q_cubed_equals_one"
        ] = False
        return record

    monkeypatch.setattr(conic, "theorem_record", broken_theorem_record)
    with pytest.raises(ArithmeticError):
        p3_odd_radon_centrality_component_upgrade(31, 6)
