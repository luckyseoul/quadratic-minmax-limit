import json
from pathlib import Path

from e1_gmin_m4_prop15745 import (
    collision_one_coordinate_bound,
    difference_radon_partition_ledger,
    moment_root_rank_audit,
    p13_t4_u0_ledger,
    proposition_15745,
    row_energy_certificate,
    translated_cut_vectors,
)


ROOT = Path(__file__).resolve().parents[1]


def test_u0_normalization_profiles_and_moments_are_exact():
    row = p13_t4_u0_ledger()
    assert row["proved"] is True
    assert row["hard_quotient_identity"] == "k_L>=1 and sum_L k_L=11"
    assert row["hard_excess_sum"] == 4
    assert row["exact_mean_14_floor_compatible_b_values"] == [2, 12]
    assert row["b2_candidate_excess_above_XNOR"] == 2
    assert row["b2_candidate_excluded_by_integral_lift_floor"] is True
    assert row["forced_exact_literal_b"] == 12
    assert row["exact_parallel_candidates"] == [5]
    assert row["hard_sign_times_global_T"] == 17
    assert row["general_hard_parallel_count"] == "P_L=4+k_L=5+e_L"
    assert row["hard_parallel_edge_total"] == 39
    assert row["opposite_parallel_edge_total"] == 22
    assert row["opposite_parallel_profile"] == [3] * 6 + [4]
    assert [
        entry["exact_literal_star_count"]
        for entry in row["hard_excess_partitions"]
    ] == [6, 5, 5, 4, 3]
    assert [
        entry["forced_global_moment_degrees"]
        for entry in row["hard_excess_partitions"]
    ] == [[2, 4], [2, 4], [2, 4], [2], [2]]


def test_full_cut_catalog_root_ranks_and_radon_bases_are_exact():
    assert len(translated_cut_vectors()) == 74
    roots = moment_root_rank_audit()
    assert roots["degree_2"] == {
        "homogeneous_degree": 2,
        "root_count": 3,
        "root_count_exceeds_degree": True,
        "subsets_checked": 364,
        "common_evaluation_rank": 3,
    }
    assert roots["degree_4"] == {
        "homogeneous_degree": 4,
        "root_count": 5,
        "root_count_exceeds_degree": True,
        "subsets_checked": 2002,
        "common_evaluation_rank": 5,
    }
    radon = difference_radon_partition_ledger()
    assert radon["proved"] is True
    assert {
        key: value["nonexact_parseval_base"]
        for key, value in radon["partition_ledgers"].items()
    } == {
        "4": 625,
        "3+1": 661,
        "2+2": 665,
        "2+1+1": 693,
        "1+1+1+1": 721,
    }


def test_broad_row_maxima_and_independent_threshold_replays_are_exact():
    expected = {
        "opposite_q3_m24": (False, None, None),
        "hard_e1_m2": (True, 31, [0, 3, 1, 4, 1, 2]),
        "hard_e2_m2": (True, 96, [5, 0, -3, -2, 3, 7]),
        "opposite_q3_m2": (True, 76, [-5, -3, -2, -5, -2, -3]),
        "opposite_q4_m2": (True, 111, [-4, -6, -1, -3, -7, 0]),
        "hard_e2_c1_m2": (True, 66, [4, 1, -3, 0, 2, 6]),
    }
    for key, (feasible, maximum, maximizer) in expected.items():
        row = row_energy_certificate(key)
        assert row["proved"] is True
        assert row["feasible"] is feasible
        assert row["sharp_energy_maximum"] == maximum
        assert row["explicit_maximizer"] == maximizer
        assert row["optimization_model"]["status"] == (
            "OPTIMAL" if feasible else "INFEASIBLE"
        )
        assert row["independent_replay_model"]["status"] == "INFEASIBLE"
        assert len(row["optimization_model"]["model_proto_sha256"]) == 64
        assert len(row["independent_replay_model"]["model_proto_sha256"]) == 64


def test_collision_one_tightening_closes_the_last_u0_partition():
    collision = collision_one_coordinate_bound()
    assert collision["proved"] is True
    assert collision["minimum_collision_contribution"] == 1
    assert collision["unique_sorted_equality_profile"] == [1, 1, 1, 1, 1, 2]
    assert collision["hard_sign_classes_per_nonzero_bin"] == 6
    assert collision["opposite_sign_classes_per_nonzero_bin"] == 7
    assert collision["hard_normalized_coordinate_interval"] == [-7, 6]
    bucket = collision["finite_field_bucket_sign_audit"]
    assert bucket["chi_minus_one"] == 1
    assert bucket["every_nonzero_F13_scalar_is_a_square_in_F13_squared"] is True
    assert bucket["hard_normalized_nonzero_buckets_checked"] == 84

    result = proposition_15745()
    assert result["proved"] is True
    assert result["p13_t4_u0_closed"] is True
    assert result["p13_k_eq_60_closed"] is False
    assert result["not_addressed_by_prop_15745"] == [3, 4, 6]
    assert result["prior_prop_15744_u3_closed"] is True
    assert result["remaining_p13_t4_residues"] == [4, 6]
    assert result["partition_close"]["three_exact_stars_1+1+1+1"]["closed"]
    assert result["partition_close"]["four_exact_stars_2+1+1"]["closed"]


def test_checked_in_evidence_matches_live_certificate():
    expected = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15745.json").read_text()
    )
    assert expected == proposition_15745()
